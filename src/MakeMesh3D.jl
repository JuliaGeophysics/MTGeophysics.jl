using ArchGDAL
using Printf

const _TOPO_CACHE = Dict{String, Any}()

# modem's covariance reader skips exactly 16 header lines - keep this 16 lines long
const _MESH_COV_HEADER = """
+-----------------------------------------------------------------------------+
| This file defines model covariance for a recursive autoregression scheme.   |
| The model space may be divided into distinct areas using integer masks.     |
| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |
| air, ocean and the rest of the model is turned off automatically. You can   |
| also define exceptions to override smoothing between any two model areas.   |
| To turn off smoothing set it to zero. This header is 16 lines long.         |
| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
| 2. Smoothing in the X direction (NzEarth real values)                       |
| 3. Smoothing in the Y direction (NzEarth real values)                       |
| 4. Vertical smoothing (1 real value)                                        |
| 5. Number of times the smoothing should be applied (1 integer >= 0)         |
| 6. Number of exceptions (1 integer >= 0)                                    |
| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 2 & 3) |
| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
+-----------------------------------------------------------------------------+
"""

#----- EM skin depth from resistivity and period -----
_skin_depth(ρ, T) = 503.0 * sqrt(ρ * T)

_snap(v, r) = collect(r)[argmin(abs.(collect(r) .- v))]

#----- find the typical station spacing -----
function nearest_neighbour_spacing(sx, sy)
    n = length(sx)
    n < 2 && return (median = 1000.0, min = 1000.0)
    dmin = fill(Inf, n)
    @inbounds for i in 1:n, j in 1:n
        i == j && continue
        d = hypot(sx[i] - sx[j], sy[i] - sy[j])
        d < dmin[i] && (dmin[i] = d)
    end
    dmin = filter(isfinite, dmin)
    return (median = median(dmin), min = minimum(dmin))
end

#----- count how many sites fall in each core cell -----
function site_occupancy(sx, sy, x_edges, y_edges, ix_core, iy_core)
    counts = Dict{Tuple{Int, Int}, Int}()
    nx = length(x_edges) - 1
    ny = length(y_edges) - 1
    for k in eachindex(sx)
        ix = clamp(searchsortedlast(x_edges, sx[k]), 1, nx)
        iy = clamp(searchsortedlast(y_edges, sy[k]), 1, ny)
        counts[(ix, iy)] = get(counts, (ix, iy), 0) + 1
    end
    maxper = isempty(counts) ? 0 : maximum(values(counts))
    incore = count(p -> first(p)[1] in ix_core && first(p)[2] in iy_core, collect(counts))
    ncore = length(ix_core) * length(iy_core)
    occ = ncore == 0 ? 0.0 : incore / ncore
    return (maxper = maxper, occupied = occ)
end

function _data_origin_latlon(d)
    lat0 = (length(d.origin) >= 1 && isfinite(d.origin[1])) ? Float64(d.origin[1]) : mean(d.loc[:, 1])
    lon0 = (length(d.origin) >= 2 && isfinite(d.origin[2])) ? Float64(d.origin[2]) : mean(d.loc[:, 2])
    return lat0, lon0
end

function _local_tm_to_target_context(d, target_crs::AbstractString)
    lat0, lon0 = _data_origin_latlon(d)
    src = "+proj=tmerc +lat_0=$(float(lat0)) +lon_0=$(float(lon0)) +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs"
    local_tm_to_target = Proj.Transformation(src, strip(target_crs); always_xy = true)
    wgs84_to_target = Proj.Transformation("EPSG:4326", strip(target_crs); always_xy = true)

    pred_x = Float64[]
    pred_y = Float64[]
    actual_x = Float64[]
    actual_y = Float64[]
    for i in 1:d.ns
        p = local_tm_to_target((Float64(d.y[i]) + 500000.0, Float64(d.x[i])))
        a = wgs84_to_target((Float64(d.loc[i, 2]), Float64(d.loc[i, 1])))
        push!(pred_x, Float64(p[1]))
        push!(pred_y, Float64(p[2]))
        push!(actual_x, Float64(a[1]))
        push!(actual_y, Float64(a[2]))
    end

    return (
        trans = local_tm_to_target,
        shift_x = median(pred_x .- actual_x),
        shift_y = median(pred_y .- actual_y),
        lat0 = lat0,
        lon0 = lon0,
    )
end

function _local_xy_to_target_xy(ctx, x_local::Real, y_local::Real)
    p = ctx.trans((Float64(y_local) + 500000.0, Float64(x_local)))
    return Float64(p[1]) - ctx.shift_x, Float64(p[2]) - ctx.shift_y
end

function _mesh_target_bbox(m, ctx)
    xt = Float64[]
    yt = Float64[]
    for x in (m.x_edges[1], m.x_edges[end]), y in (m.y_edges[1], m.y_edges[end])
        px, py = _local_xy_to_target_xy(ctx, x, y)
        push!(xt, px)
        push!(yt, py)
    end
    return (minimum(xt), maximum(xt), minimum(yt), maximum(yt))
end

function _data_target_bbox(d, ctx)
    xt = Float64[]
    yt = Float64[]
    for i in 1:d.ns
        px, py = _local_xy_to_target_xy(ctx, d.x[i], d.y[i])
        push!(xt, px)
        push!(yt, py)
    end
    return (minimum(xt), maximum(xt), minimum(yt), maximum(yt))
end

function _load_topography_geotiff(path::AbstractString)
    ArchGDAL.read(path) do dataset
        band = ArchGDAL.getband(dataset, 1)
        zgrid = Float64.(ArchGDAL.read(band))
        nodata = ArchGDAL.getnodatavalue(band)
        gt = ArchGDAL.getgeotransform(dataset)
        x0, dx, rx, y0, ry, dy = Float64.(gt)

        (abs(rx) > 1e-9 || abs(ry) > 1e-9) && error("Rotated GeoTIFF geotransforms are not supported by this mesh builder.")

        ny, nx = size(zgrid)
        xs = x0 .+ ((0:(nx - 1)) .+ 0.5) .* dx
        ys = y0 .+ ((0:(ny - 1)) .+ 0.5) .* dy

        if xs[end] < xs[1]
            xs = reverse(xs)
            zgrid = reverse(zgrid; dims = 2)
        end
        if ys[end] < ys[1]
            ys = reverse(ys)
            zgrid = reverse(zgrid; dims = 1)
        end

        if nodata !== nothing
            nd = Float64(nodata)
            zgrid[zgrid .== nd] .= NaN
        end
        # Some DEM products use extreme finite sentinels instead of explicit NaN.
        zgrid[abs.(zgrid) .> 1e30] .= NaN
        zgrid[.!isfinite.(zgrid)] .= NaN

        valid_rows = vec(any(isfinite.(zgrid), dims = 2))
        valid_cols = vec(any(isfinite.(zgrid), dims = 1))
        any(valid_rows) || error("Topography GeoTIFF contains no finite elevation samples.")
        any(valid_cols) || error("Topography GeoTIFF contains no finite elevation samples.")
        valid_bbox = (
            minimum(xs[valid_cols]), maximum(xs[valid_cols]),
            minimum(ys[valid_rows]), maximum(ys[valid_rows]),
        )

        return (
            x = collect(xs),
            y = collect(ys),
            z = zgrid,
            bbox = (minimum(xs), maximum(xs), minimum(ys), maximum(ys)),
            valid_bbox = valid_bbox,
        )
    end
end

function _ensure_topography(path::AbstractString, bbox)
    haskey(_TOPO_CACHE, path) && return _TOPO_CACHE[path]
    isfile(path) || error("Topography file not found: $path")
    return _TOPO_CACHE[path] = _load_topography_geotiff(path)
end

function _sample_topography(topo, xq::Real, yq::Real)
    xs = topo.x
    ys = topo.y
    z = topo.z
    nx = length(xs)
    ny = length(ys)

    # Padding cells can extend beyond the DEM. Clamp them to the raster edge so
    # the topography over the padding decays to the nearest available elevation.
    xlo, xhi, ylo, yhi = topo.valid_bbox
    xqc = clamp(Float64(xq), xlo, xhi)
    yqc = clamp(Float64(yq), ylo, yhi)

    ix = xqc == xs[end] ? nx - 1 : clamp(searchsortedlast(xs, xqc), 1, nx - 1)
    iy = yqc == ys[end] ? ny - 1 : clamp(searchsortedlast(ys, yqc), 1, ny - 1)

    x1, x2 = xs[ix], xs[ix + 1]
    y1, y2 = ys[iy], ys[iy + 1]
    tx = x2 == x1 ? 0.0 : (xqc - x1) / (x2 - x1)
    ty = y2 == y1 ? 0.0 : (yqc - y1) / (y2 - y1)

    z11 = z[iy, ix]
    z12 = z[iy, ix + 1]
    z21 = z[iy + 1, ix]
    z22 = z[iy + 1, ix + 1]
    if all(isfinite, (z11, z12, z21, z22))
        return (1 - tx) * (1 - ty) * z11 + tx * (1 - ty) * z12 + (1 - tx) * ty * z21 + tx * ty * z22
    end

    candidates = [(z11, hypot(xqc - x1, yqc - y1)), (z12, hypot(xqc - x2, yqc - y1)),
                  (z21, hypot(xqc - x1, yqc - y2)), (z22, hypot(xqc - x2, yqc - y2))]
    finite_candidates = filter(c -> isfinite(c[1]), candidates)
    if !isempty(finite_candidates)
        return first(finite_candidates[argmin(last.(finite_candidates))])
    end

    # If the immediate bilinear stencil falls in a no-data moat, expand locally
    # until the nearest finite DEM sample is found.
    for radius in 1:25
        ix1 = max(1, ix - radius)
        ix2 = min(nx, ix + radius + 1)
        iy1 = max(1, iy - radius)
        iy2 = min(ny, iy + radius + 1)
        best_val = NaN
        best_dist = Inf
        for jy in iy1:iy2, jx in ix1:ix2
            val = z[jy, jx]
            isfinite(val) || continue
            dist = hypot(xqc - xs[jx], yqc - ys[jy])
            if dist < best_dist
                best_dist = dist
                best_val = val
            end
        end
        isfinite(best_val) && return best_val
    end
    return NaN
end

function _sample_surface_local(m, d, ctx, topo_path::AbstractString)
    bbox = _mesh_target_bbox(m, ctx)
    topo = _ensure_topography(topo_path, bbox)

    datum_candidates = Float64[]
    finite_station_elev = Float64[]
    for i in eachindex(d.x)
        tx, ty = _local_xy_to_target_xy(ctx, d.x[i], d.y[i])
        elev_abs = _sample_topography(topo, tx, ty)
        isfinite(elev_abs) || continue
        push!(finite_station_elev, elev_abs)
        push!(datum_candidates, d.z[i] + elev_abs)
    end
    isempty(datum_candidates) && error("Could not estimate the local vertical datum from the topography file.")
    datum_elev = median(datum_candidates)
    fallback_elev = median(finite_station_elev)

    surface_local = Matrix{Float64}(undef, m.nx, m.ny)
    fallback_count = 0
    for ix in 1:m.nx, iy in 1:m.ny
        tx, ty = _local_xy_to_target_xy(ctx, m.x_centers[ix], m.y_centers[iy])
        elev_abs = _sample_topography(topo, tx, ty)
        if !isfinite(elev_abs)
            elev_abs = fallback_elev
            fallback_count += 1
        end
        surface_local[ix, iy] = datum_elev - elev_abs
    end

    fallback_count == 0 || @warn("Used median DEM elevation fallback for padded mesh columns that still sampled nodata.", fallback_count = fallback_count)

    return surface_local, datum_elev
end

function _derive_topography_data(d, ctx, topo_path::AbstractString)
    bbox = _data_target_bbox(d, ctx)
    topo = _ensure_topography(topo_path, bbox)

    station_elev_abs = Vector{Float64}(undef, d.ns)
    for i in 1:d.ns
        tx, ty = _local_xy_to_target_xy(ctx, d.x[i], d.y[i])
        elev_abs = _sample_topography(topo, tx, ty)
        isfinite(elev_abs) || error("Could not derive a finite DEM elevation for site $(d.site[i]).")
        station_elev_abs[i] = elev_abs
    end

    # Use the highest sampled station elevation as the local vertical datum so the
    # derived site Z values depend only on the DEM, not on the input data file.
    datum_elev = maximum(station_elev_abs)
    topo_z = datum_elev .- station_elev_abs
    delta_z = topo_z .- Float64.(d.z)

    d_topo = deepcopy(d)
    d_topo.z = collect(Float64.(topo_z))
    d_topo.loc[:, 3] .= topo_z

    return d_topo, (
        ns = d.ns,
        datum_elev = datum_elev,
        zmin = minimum(topo_z),
        zmax = maximum(topo_z),
        nchanged = count(abs.(delta_z) .> 1e-6),
        median_shift = median(abs.(delta_z)),
        max_shift = maximum(abs.(delta_z)),
    )
end

#----- build the mesh from sites, padding and depth settings -----
function design_mesh(sx, sy, T; ρ_bg, dx_core, dy_core, nx_pad, ny_pad,
                     pad_factor, z_first, z_factor, depth_mult,
                     nx_core = nothing, ny_core = nothing,
                     surface_local_range = nothing, n_air = 0,
                     air_factor = 1.35, air_first = nothing,
                     air_first_div = 2.0)
    dxc = max(1.0, dx_core)
    dyc = max(1.0, dy_core)
    nx_pad = max(0, nx_pad)
    ny_pad = max(0, ny_pad)

    cx_mid = (minimum(sx) + maximum(sx)) / 2
    cy_mid = (minimum(sy) + maximum(sy)) / 2
    nx_core = nx_core === nothing ? max(1, ceil(Int, max(maximum(sx) - minimum(sx), dxc) / dxc)) : max(1, Int(nx_core))
    ny_core = ny_core === nothing ? max(1, ceil(Int, max(maximum(sy) - minimum(sy), dyc) / dyc)) : max(1, Int(ny_core))
    x0 = cx_mid - nx_core * dxc / 2
    y0 = cy_mid - ny_core * dyc / 2

    padx = [dxc * pad_factor^i for i in 1:nx_pad]
    pady = [dyc * pad_factor^i for i in 1:ny_pad]
    dx = vcat(reverse(padx), fill(dxc, nx_core), padx)
    dy = vcat(reverse(pady), fill(dyc, ny_core), pady)

    origin_x = x0 - sum(padx)
    origin_y = y0 - sum(pady)

    Tmin, Tmax = minimum(T), maximum(T)
    δmin = _skin_depth(ρ_bg, Tmin)
    δmax = _skin_depth(ρ_bg, Tmax)
    z_target = depth_mult * δmax

    ground_dz = Float64[]
    z = 0.0
    t = z_first
    while z < z_target
        push!(ground_dz, t)
        z += t
        t *= z_factor
    end

    air_dz = Float64[]
    origin_z = 0.0
    if surface_local_range !== nothing && Int(n_air) > 0
        surface_min = Float64(surface_local_range[1])
        air_seed = air_first === nothing ? max(5.0, z_first / air_first_div) : max(1.0, Float64(air_first))
        air_dz = reverse([air_seed * air_factor^(i - 1) for i in 1:Int(n_air)])
        origin_z = surface_min - sum(air_dz)
    end

    dz = vcat(air_dz, ground_dz)
    x_edges = origin_x .+ vcat(0.0, cumsum(dx))
    y_edges = origin_y .+ vcat(0.0, cumsum(dy))
    z_edges = origin_z .+ vcat(0.0, cumsum(dz))
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    y_centers = (y_edges[1:end-1] .+ y_edges[2:end]) ./ 2
    z_centers = (z_edges[1:end-1] .+ z_edges[2:end]) ./ 2

    ix_core = (nx_pad + 1):(nx_pad + nx_core)
    iy_core = (ny_pad + 1):(ny_pad + ny_core)
    diag = site_occupancy(sx, sy, x_edges, y_edges, ix_core, iy_core)

    core_x0 = x_edges[nx_pad + 1]
    core_x1 = x_edges[nx_pad + nx_core + 1]
    core_y0 = y_edges[ny_pad + 1]
    core_y1 = y_edges[ny_pad + ny_core + 1]
    pad_extent = min(sum(padx), sum(pady))
    nz_air = length(air_dz)
    nz_earth = length(ground_dz)
    surface_min = surface_local_range === nothing ? 0.0 : Float64(surface_local_range[1])
    surface_max = surface_local_range === nothing ? 0.0 : Float64(surface_local_range[2])

    return (
        dx = dx, dy = dy, dz = dz, origin = [origin_x, origin_y, origin_z], ρ_bg = ρ_bg,
        nx = length(dx), ny = length(dy), nz = length(dz), nz_air = nz_air, nz_earth = nz_earth,
        nx_core = nx_core, ny_core = ny_core, dx_core = dxc, dy_core = dyc,
        x_edges = x_edges, y_edges = y_edges, z_edges = z_edges,
        x_centers = x_centers, y_centers = y_centers, z_centers = z_centers,
        x_edges_km = x_edges ./ 1000, y_edges_km = y_edges ./ 1000, z_edges_km = z_edges ./ 1000,
        core_x0_km = core_x0 / 1000, core_x1_km = core_x1 / 1000,
        core_y0_km = core_y0 / 1000, core_y1_km = core_y1 / 1000,
        full_x0_km = x_edges[1] / 1000, full_x1_km = x_edges[end] / 1000,
        full_y0_km = y_edges[1] / 1000, full_y1_km = y_edges[end] / 1000,
        δmin_km = δmin / 1000, δmax_km = δmax / 1000,
        depth_km = sum(ground_dz) / 1000, pad_extent_km = pad_extent / 1000,
        surface_min_km = surface_min / 1000, surface_max_km = surface_max / 1000,
        pad_ok = pad_extent >= δmax,
        maxper = diag.maxper, occupied = diag.occupied,
    )
end

function _build_resistivity_model(m, surface_local)
    A = fill(Float64(m.ρ_bg), m.nx, m.ny, m.nz)
    surface_local === nothing && return A
    @inbounds for ix in 1:m.nx, iy in 1:m.ny, iz in 1:m.nz
        if m.z_centers[iz] < surface_local[ix, iy]
            A[ix, iy, iz] = NaN
        end
    end
    return A
end

function _build_covariance_mask(m, surface_local)
    mask = fill(Int8(1), m.nx, m.ny, m.nz_earth)
    surface_local === nothing && return mask
    @inbounds for ke in 1:m.nz_earth
        iz = m.nz_air + ke
        zc = m.z_centers[iz]
        for ix in 1:m.nx, iy in 1:m.ny
            if zc < surface_local[ix, iy]
                mask[ix, iy, ke] = 0
            end
        end
    end
    return mask
end

function _write_start_model_ws(path::AbstractString, m; rotation::Real = 0.0)
    nz_air = m.nz_air
    nz_earth = m.nz - nz_air
    nz_earth > 0 || error("Model must contain at least one earth layer.")

    Aw = copy(Float64.(m.A[:, :, (nz_air + 1):end]))
    Aw[isnan.(Aw)] .= 1e17
    Aw .= log.(Aw)

    origin_z = m.origin[3] + sum(m.dz[1:nz_air])

    open(path, "w") do io
        println(io, "# Written by MTGeophysics.jl MakeMesh3D WS-compatible earth-only model")
        println(io, "$(m.nx) $(m.ny) $(nz_earth) 0 LOGE")

        for v in m.dx; print(io, "$v "); end; println(io)
        for v in m.dy; print(io, "$v "); end; println(io)
        for v in m.dz[(nz_air + 1):end]; print(io, "$v "); end; println(io)

        for k in 1:nz_earth
            println(io)
            for j in 1:m.ny
                for i in m.nx:-1:1
                    @printf(io, "%15.5E", Aw[i, j, k])
                end
                println(io)
            end
        end

        println(io, "$(m.origin[1]) $(m.origin[2]) $(origin_z)")
        println(io, "$rotation")
    end

    println("Model written to: $path")
    return path
end

function _print_real_line(io, values)
    for v in values
        @printf(io, " %.3f ", Float64(v))
    end
    println(io)
end

function write_covariance_file(path::AbstractString, mask::Array{<:Integer, 3};
                               cov_x::Real, cov_y::Real, cov_z::Real,
                               n_apply::Integer = 2,
                               exceptions::Vector{Tuple{Int, Int, Float64}} = Tuple{Int, Int, Float64}[])
    nx, ny, nz = size(mask)
    open(path, "w") do io
        print(io, _MESH_COV_HEADER)
        println(io)
        @printf(io, " %d %12d %12d \n\n", nx, ny, nz)
        _print_real_line(io, fill(Float64(cov_x), nz))
        _print_real_line(io, fill(Float64(cov_y), nz))
        _print_real_line(io, (Float64(cov_z),))
        println(io)
        @printf(io, " %d \n\n", Int(n_apply))
        @printf(io, " %d \n", length(exceptions))
        for (a, b, v) in exceptions
            @printf(io, " %d %d %.3f\n", a, b, v)
        end
        println(io)
        for k in 1:nz
            @printf(io, "\n %d %12d \n", k, k)
            for ix in 1:nx
                for iy in 1:ny
                    @printf(io, " %d ", Int(mask[ix, iy, k]))
                end
                println(io)
            end
        end
    end
    return path
end

function _earth_only_z_edges(m)
    z0 = m.origin[3] + sum(m.dz[1:m.nz_air])
    return z0 .+ vcat(0.0, cumsum(m.dz[(m.nz_air + 1):end]))
end

function _site_column_constraints(m, d)
    cols = Dict{Tuple{Int, Int}, Tuple{String, Float64}}()
    for isite in 1:d.ns
        ix = clamp(searchsortedlast(m.x_edges, d.x[isite]), 1, m.nx)
        iy = clamp(searchsortedlast(m.y_edges, d.y[isite]), 1, m.ny)
        key = (ix, iy)
        site_z = Float64(d.z[isite])
        current = get(cols, key, ("", Inf))
        if site_z < current[2]
            cols[key] = (String(d.site[isite]), site_z)
        end
    end
    return cols
end

function _apply_site_surface_constraints(m, d, surface_local)
    surface_local === nothing && return nothing, (
        nadjusted = 0,
        adjustments = NamedTuple[],
        impossible = NamedTuple[],
        text = "Covariance follows topography only.",
    )

    earth_z_edges = _earth_only_z_edges(m)
    earth_z_centers = Float64.(m.z_centers[(m.nz_air + 1):end])
    adjusted = copy(surface_local)
    adjustments = NamedTuple[]
    impossible = NamedTuple[]

    for ((ix, iy), (site_name, site_z)) in _site_column_constraints(m, d)
        max_first = searchsortedlast(earth_z_edges, site_z)
        if max_first < 1
            push!(impossible, (
                site = site_name,
                ix = ix,
                iy = iy,
                site_z = site_z,
                top_edge = Float64(earth_z_edges[1]),
                clearance = site_z - Float64(earth_z_edges[1]),
            ))
            continue
        end

        desired_first = min(max_first, m.nz_earth)
        old_surface = Float64(adjusted[ix, iy])
        new_surface = min(old_surface, earth_z_centers[desired_first])
        adjusted[ix, iy] = new_surface

        if new_surface < old_surface - 1e-9
            push!(adjustments, (
                site = site_name,
                ix = ix,
                iy = iy,
                site_z = site_z,
                surface_from = old_surface,
                surface_to = new_surface,
                lift = old_surface - new_surface,
            ))
        end
    end

    sort!(adjustments; by = row -> -row.lift)
    sort!(impossible; by = row -> row.clearance)
    text = isempty(adjustments) ?
        "Covariance follows topography only." :
        @sprintf("Covariance auto-adjusted in %d station column(s) to honour site elevations.", length(adjustments))

    return adjusted, (
        nadjusted = length(adjustments),
        adjustments = adjustments,
        impossible = impossible,
        text = text,
    )
end

function _site_air_report(m, d)
    earth_z_edges = _earth_only_z_edges(m)
    bad = NamedTuple[]

    for isite in 1:d.ns
        ix = clamp(searchsortedlast(m.x_edges, d.x[isite]), 1, m.nx)
        iy = clamp(searchsortedlast(m.y_edges, d.y[isite]), 1, m.ny)
        first_earth = findfirst(!=(0), vec(m.cov_mask[ix, iy, :]))
        first_earth === nothing && error("No earth cells remain in covariance column for site $(d.site[isite]) at mesh column ($ix, $iy).")

        top_edge = earth_z_edges[first_earth]
        clearance = Float64(d.z[isite]) - Float64(top_edge)
        if clearance < 0.0
            push!(bad, (
                site = String(d.site[isite]),
                ix = ix,
                iy = iy,
                site_z = Float64(d.z[isite]),
                top_edge = Float64(top_edge),
                clearance = clearance,
            ))
        end
    end

    sort!(bad; by = row -> row.clearance)
    nbad = length(bad)
    ok = nbad == 0
    sample_rows = bad[1:min(end, 5)]
    sample = isempty(sample_rows) ? "" : join([
        @sprintf("%s needs %.1f m", row.site, -row.clearance) for row in sample_rows
    ], ", ")
    max_raise = isempty(bad) ? 0.0 : -bad[1].clearance
    summary = ok ?
        @sprintf("✓ All %d sites lie at or below the first non-air covariance cell.", d.ns) :
        @sprintf("• %d site(s) lie in air relative to the written covariance mask. Raise the local model surface at those XY columns until the top of the first non-air cell is at or above each site Z. Example adjustments: %s.",
            nbad, sample)
    fix = ok ? "" :
        "Fix the model, not the data: remove leading air-mask 0 values in the flagged covariance columns and raise the matching topography/model surface until every flagged site has nonnegative clearance (site_z - top_of_first_earth >= 0)."

    return (
        ok = ok,
        nbad = nbad,
        rows = bad,
        sample = sample,
        max_raise = max_raise,
        text = summary,
        fix = fix,
    )
end

function _print_site_surface_report(m)
    if m.site_surface.nadjusted > 0
        println(m.site_surface.text)
        for row in m.site_surface.adjustments[1:min(end, 10)]
            println(@sprintf(
                "  %s: xcell=%d ycell=%d surface %.3f -> %.3f (lift %.3f m)",
                row.site, row.ix, row.iy, row.surface_from, row.surface_to, row.lift))
        end
        if m.site_surface.nadjusted > 10
            println(@sprintf("  ... %d more adjusted column(s) omitted.", m.site_surface.nadjusted - 10))
        end
    end

    if !isempty(m.site_surface.impossible)
        @warn(@sprintf("%d site(s) still sit above the shallowest earth cell; revise the vertical mesh or topography datum.", length(m.site_surface.impossible)))
        for row in m.site_surface.impossible[1:min(end, 10)]
            println(@sprintf(
                "  %s: xcell=%d ycell=%d site_z=%.3f shallowest_earth=%.3f clearance=%.3f",
                row.site, row.ix, row.iy, row.site_z, row.top_edge, row.clearance))
        end
        if length(m.site_surface.impossible) > 10
            println(@sprintf("  ... %d more impossible site(s) omitted.", length(m.site_surface.impossible) - 10))
        end
    end
    return nothing
end

function _print_site_air_report(m)
    m.site_air.ok && return nothing

    @warn(@sprintf("%d site(s) are in air relative to the written covariance mask.", m.site_air.nbad))
    println(m.site_air.fix)
    println(@sprintf("Worst-case upward surface adjustment needed: %.1f m.", m.site_air.max_raise))
    for row in m.site_air.rows[1:min(end, 10)]
        println(@sprintf(
            "  %s: xcell=%d ycell=%d site_z=%.3f top_of_first_earth=%.3f clearance=%.3f",
            row.site, row.ix, row.iy, row.site_z, row.top_edge, row.clearance))
    end
    if m.site_air.nbad > 10
        println(@sprintf("  ... %d more site(s) omitted.", m.site_air.nbad - 10))
    end
    return nothing
end

function build_mesh_bundle(d, sx, sy, Tobs;
                           ρ_bg, dx_core, dy_core, nx_core, ny_core,
                           nx_pad, ny_pad, pad_factor, z_first, z_factor, depth_mult,
                           cov_value, topo_ctx = nothing,
                           topo_path::AbstractString = "",
                           n_air::Int = 6, air_factor::Real = 1.35,
                           air_first_div::Real = 2.0)
    use_topo = !isempty(strip(topo_path))
    base = design_mesh(sx, sy, Tobs;
        ρ_bg = Float64(ρ_bg),
        dx_core = Float64(dx_core),
        dy_core = Float64(dy_core),
        nx_core = Int(nx_core),
        ny_core = Int(ny_core),
        nx_pad = Int(nx_pad),
        ny_pad = Int(ny_pad),
        pad_factor = Float64(pad_factor),
        z_first = Float64(z_first),
        z_factor = Float64(z_factor),
        depth_mult = Float64(depth_mult))

    surface_local = nothing
    datum_elev = NaN
    if use_topo
        topo_ctx === nothing && error("Topography mode was requested without a valid projection context.")
        surface_local, datum_elev = _sample_surface_local(base, d, topo_ctx, topo_path)
        base = design_mesh(sx, sy, Tobs;
            ρ_bg = Float64(ρ_bg),
            dx_core = Float64(dx_core),
            dy_core = Float64(dy_core),
            nx_core = Int(nx_core),
            ny_core = Int(ny_core),
            nx_pad = Int(nx_pad),
            ny_pad = Int(ny_pad),
            pad_factor = Float64(pad_factor),
            z_first = Float64(z_first),
            z_factor = Float64(z_factor),
            depth_mult = Float64(depth_mult),
            surface_local_range = (minimum(surface_local), maximum(surface_local)),
            n_air = n_air,
            air_factor = Float64(air_factor),
            air_first = max(5.0, Float64(z_first) / air_first_div))
    end

    site_surface = (
        nadjusted = 0,
        adjustments = NamedTuple[],
        impossible = NamedTuple[],
        text = "Covariance follows topography only.",
    )
    if surface_local !== nothing
        surface_local, site_surface = _apply_site_surface_constraints(base, d, surface_local)
    end

    A = _build_resistivity_model(base, surface_local)
    cov_mask = _build_covariance_mask(base, surface_local)
    site_air = _site_air_report((; base..., A = A, cov_mask = cov_mask), d)
    return (; base..., A = A, cov_mask = cov_mask, cov_value = Float64(cov_value),
              topo_active = surface_local !== nothing,
              topo_name = surface_local === nothing ? "flat" : basename(topo_path),
              datum_elev = datum_elev,
              site_surface = site_surface,
              site_air = site_air)
end

function _save_mesh_outputs(m, out_model::AbstractString, out_cov::AbstractString;
                            cov_apply::Integer = 2)
    mkpath(dirname(abspath(out_model)))
    mkpath(dirname(abspath(out_cov)))
    _write_start_model_ws(out_model, m; rotation = 0.0)
    write_covariance_file(out_cov, m.cov_mask;
        cov_x = m.cov_value, cov_y = m.cov_value, cov_z = m.cov_value,
        n_apply = cov_apply)
    _print_site_surface_report(m)
    _print_site_air_report(m)
    return nothing
end

function MakeMesh3D(data_file::AbstractString;
    out_model::AbstractString = joinpath(dirname(abspath(data_file)), "mesh_start_model.rho"),
    out_covariance::AbstractString = joinpath(dirname(out_model), "C3.dat"),
    topo_file::AbstractString = "",
    topo_crs::AbstractString = "EPSG:26918",
    cov_value::Real = 0.3,
    cov_apply::Int = 2,
    mode::Symbol = :gui,
    cell_width_frac::Real = 0.5,
    n_pad::Int = 12,
    pad_factor::Real = 1.5,
    first_layer_div::Real = 4.0,
    vertical_factor::Real = 1.2,
    depth_mult::Real = 3.0,
    air_layers::Int = 6,
    air_factor::Real = 1.35,
    air_first_div::Real = 2.0,
    colormap = :Spectral,
    resistivity_range = (0.0, 4.0),
    site_color = :black,
    site_size_full = 6,
    site_size_core = 5,
    grid_color = (:grey, 0.7),
    grid_linewidth = 1.0,
    fig_size = (1500, 1000),
    interactive::Bool = true)

    mode in (:gui, :nogui) || error("mode must be :gui or :nogui, got $mode")
    isempty(data_file) && error("data_file is required")
    isfile(data_file) || error("Data file not found: $data_file")

    use_topo = !isempty(strip(topo_file))
    use_topo && (isfile(topo_file) || error("Topography file not found: $topo_file"))

    d_raw = load_data_modem(data_file)
    d_raw.name = data_file
    d_raw.origin = collect(Float64.(d_raw.origin))
    d_raw.loc = Float64.(d_raw.loc)
    topo_ctx = _local_tm_to_target_context(d_raw, topo_crs)

    if use_topo
        d, d_topo_diag = _derive_topography_data(d_raw, topo_ctx, topo_file)
    else
        d = d_raw
        d_topo_diag = nothing
    end
    sx = collect(Float64.(d.x))
    sy = collect(Float64.(d.y))
    Tobs = collect(Float64.(d.T))

    spacing = nearest_neighbour_spacing(sx, sy)
    ρ_off = filter(x -> isfinite(x) && x > 0, vec(d.ρ[:, [2, 3], :]))
    ρ_bg_data = isempty(ρ_off) ? 100.0 : median(ρ_off)
    ρ_bg0 = _snap(ρ_bg_data, 10:10:20000)
    dx_core0 = _snap(spacing.median * cell_width_frac, 50:25:10000)
    dy_core0 = dx_core0
    span_x0 = max(maximum(sx) - minimum(sx), 1.0)
    span_y0 = max(maximum(sy) - minimum(sy), 1.0)
    nx_core0 = max(1, ceil(Int, span_x0 / dx_core0))
    ny_core0 = max(1, ceil(Int, span_y0 / dy_core0))
    z_first0 = _snap(_skin_depth(Float64(ρ_bg0), minimum(Tobs)) / first_layer_div, 5:5:2000)

    @info @sprintf("Loaded %d sites, %d periods (%.3g–%.3g s); median spacing %.0f m; suggested background ρ = %.0f Ω·m; topography %s",
                   length(sx), length(Tobs), minimum(Tobs), maximum(Tobs), spacing.median, ρ_bg_data,
                   use_topo ? basename(topo_file) : "off")

    if mode === :nogui
        m = build_mesh_bundle(d, sx, sy, Tobs;
            ρ_bg = Float64(ρ_bg0), dx_core = Float64(dx_core0), dy_core = Float64(dy_core0),
            nx_core = nx_core0, ny_core = ny_core0,
            nx_pad = n_pad, ny_pad = n_pad,
            pad_factor = pad_factor, z_first = Float64(z_first0), z_factor = vertical_factor,
            depth_mult = depth_mult, cov_value = cov_value, topo_ctx = topo_ctx,
            topo_path = topo_file, n_air = air_layers, air_factor = air_factor,
            air_first_div = air_first_div)
        @printf("grid %d×%d×%d (%d cells); core %d×%d @ %.0f×%.0f m; air %d; max sites/cell %d; occupied %.0f%%; pad %.0f km/side; depth %.0f km; cov %.2f\n",
            m.nx, m.ny, m.nz, m.nx * m.ny * m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core,
            m.nz_air, m.maxper, 100 * m.occupied, m.pad_extent_km, m.depth_km, m.cov_value)
        _save_mesh_outputs(m, out_model, out_covariance; cov_apply = cov_apply)
        @printf("wrote %s\n", out_model)
        @printf("wrote %s\n", out_covariance)
        return m
    end

    isdefined(MTGeophysics, :_make_mesh3D_gui) ||
        error("GLMakie is not available; use mode = :nogui on headless systems")

    return _make_mesh3D_gui(;
        data_file = data_file,
        out_model = out_model,
        out_covariance = out_covariance,
        topo_file = topo_file,
        cov_value = cov_value,
        cov_apply = cov_apply,
        n_pad = n_pad,
        pad_factor = pad_factor,
        vertical_factor = vertical_factor,
        depth_mult = depth_mult,
        air_layers = air_layers,
        air_factor = air_factor,
        air_first_div = air_first_div,
        colormap = colormap,
        resistivity_range = resistivity_range,
        site_color = site_color,
        site_size_full = site_size_full,
        site_size_core = site_size_core,
        grid_color = grid_color,
        grid_linewidth = grid_linewidth,
        fig_size = fig_size,
        interactive = interactive,
        d = d, sx = sx, sy = sy, Tobs = Tobs,
        spacing = spacing, ρ_bg_data = ρ_bg_data,
        ρ_bg0 = ρ_bg0, dx_core0 = dx_core0, dy_core0 = dy_core0,
        nx_core0 = nx_core0, ny_core0 = ny_core0, z_first0 = z_first0,
        topo_ctx = topo_ctx)
end

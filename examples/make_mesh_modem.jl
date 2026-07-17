# Interactive ModEM mesh builder with DEM-derived station elevations.
# Inputs (optional CLI args):
#   [1] data_file (ModEM .dat)
#   [2] out_model (.rho)
#   [3] topo_geotiff (optional EPSG:26918 topography)
#   [4] out_covariance (optional, defaults to C3.dat beside the model)
#   [5] covariance value (optional, defaults to 0.3)
#   [6] mode (optional: gui|nogui)
# This variant also writes a new `*_topo.dat` copy whose station elevations are
# derived from the GeoTIFF/topography workflow instead of the input site Z values.
# Set ENV["MTGEO_MESH_MODE"] = "nogui" to build+write the mesh headlessly.

using MTGeophysics
using GLMakie
using ArchGDAL
using Proj
using Statistics
using Printf
using Dates

GLMakie.activate!()

#----- user controls (edit here) -----
const MODE            = length(ARGS) >= 6 ? lowercase(ARGS[6]) : lowercase(get(ENV, "MTGEO_MESH_MODE", "gui"))
const DATA_PATH       = length(ARGS) >= 1 ? ARGS[1] : normpath(@__DIR__, "MT3DINV4", "data7O10Dg.dat")
const OUT_PATH        = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(DATA_PATH), "mesh_start_model.rho")
const DEFAULT_TOPO_PATH = normpath(@__DIR__, "MT3DINV4", "ArcticDEM_30m_EPSG26918.tif")
const TOPO_PATH       = length(ARGS) >= 3 ? ARGS[3] : (isfile(DEFAULT_TOPO_PATH) ? DEFAULT_TOPO_PATH : "")
const OUT_COV_PATH    = length(ARGS) >= 4 ? ARGS[4] : joinpath(dirname(OUT_PATH), "C3.dat")
const OUT_DATA_PATH   = let base = basename(DATA_PATH)
    outbase = if occursin(r"(?i)_topo\.dat$", base)
        base
    elseif occursin(r"(?i)\.dat$", base)
        replace(base, r"(?i)\.dat$" => "_topo.dat")
    else
        base * "_topo.dat"
    end
    joinpath(dirname(DATA_PATH), outbase)
end
const COV_VALUE0      = length(ARGS) >= 5 ? something(tryparse(Float64, ARGS[5]), 0.3) : 0.3
const CELL_WIDTH_FRAC = 0.5
const N_PAD           = 12
const PAD_FACTOR      = 1.5
const FIRST_LAYER_DIV = 4.0
const VERTICAL_FACTOR = 1.2
const DEPTH_MULT      = 3.0
const TOPO_CRS        = "EPSG:26918"
const AIR_LAYERS      = 6
const AIR_FACTOR      = 1.35
const AIR_FIRST_DIV   = 2.0
const COV_APPLY       = 2
const CMAP            = :Spectral
const CRANGE          = (0.0, 4.0)
const SITE_COLOR      = :black
const SITE_SIZE_FULL  = 6
const SITE_SIZE_CORE  = 5
const GRID_COLOR      = (:grey, 0.7)
const GRID_WIDTH      = 1.0
const FIG_SIZE        = (1500, 1000)
const USE_TOPO        = !isempty(strip(TOPO_PATH))

const TOPO_CACHE = Ref{Any}(nothing)

USE_TOPO || error("make_mesh_modem2.jl requires a topography GeoTIFF because station elevations are derived from the DEM.")

# modem's covariance reader skips exactly 16 header lines - keep this 16 lines long
const COV_HEADER = """
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
skin_depth(ρ, T) = 503.0 * sqrt(ρ * T)

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
    topo = TOPO_CACHE[]
    if topo === nothing
        isfile(path) || error("Topography file not found: $path")
        TOPO_CACHE[] = _load_topography_geotiff(path)
    end
    return TOPO_CACHE[]
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

function _sample_surface_local(m, d, ctx)
    bbox = _mesh_target_bbox(m, ctx)
    topo = _ensure_topography(TOPO_PATH, bbox)

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

function _derive_topography_data(d, ctx)
    bbox = _data_target_bbox(d, ctx)
    topo = _ensure_topography(TOPO_PATH, bbox)

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
                     air_factor = AIR_FACTOR, air_first = nothing)
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
    δmin = skin_depth(ρ_bg, Tmin)
    δmax = skin_depth(ρ_bg, Tmax)
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
        air_seed = air_first === nothing ? max(5.0, z_first / AIR_FIRST_DIV) : max(1.0, Float64(air_first))
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
        println(io, "# Written by MTGeophysics.jl make_mesh_modem WS-compatible earth-only model")
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
                               n_apply::Integer = COV_APPLY,
                               exceptions::Vector{Tuple{Int, Int, Float64}} = Tuple{Int, Int, Float64}[])
    nx, ny, nz = size(mask)
    open(path, "w") do io
        print(io, COV_HEADER)
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

function _has_finite_complex(A)
    return any(isfinite.(real.(A))) || any(isfinite.(imag.(A)))
end

function _print_topography_data_report(diag)
    println(@sprintf(
        "Topography-derived data: %s | datum %.3f m | site z range %.3f to %.3f m | updated %d/%d site elevations | median |Δz| %.3f m | max |Δz| %.3f m",
        basename(OUT_DATA_PATH), diag.datum_elev, diag.zmin, diag.zmax,
        diag.nchanged, diag.ns, diag.median_shift, diag.max_shift))
    return nothing
end

function _format_like_token(value::Real, template::AbstractString)
    t = strip(template)
    if occursin('E', t) || occursin('e', t)
        mantissa = split(replace(t, 'e' => 'E'), 'E')[1]
        decimals = occursin('.', mantissa) ? length(split(mantissa, '.')[2]) : 0
        return @sprintf("%.*E", decimals, Float64(value))
    end
    if occursin('.', t)
        decimals = length(split(t, '.')[2])
        return @sprintf("%.*f", decimals, Float64(value))
    end
    return string(round(Int, value))
end

function _is_data_line(line::AbstractString)
    t = strip(line)
    return !isempty(t) && !startswith(t, "#") && !startswith(t, ">") && occursin(r"^[\s\+\-\.0-9]", t)
end

function _write_topography_data_copy(input_path::AbstractString, output_path::AbstractString, d_topo::Data)
    site_z = Dict{String, Float64}(d_topo.site[i] => Float64(d_topo.z[i]) for i in 1:d_topo.ns)
    open(input_path, "r") do io_in
        open(output_path, "w") do io_out
            for line in eachline(io_in)
                t = strip(line)
                if _is_data_line(t)
                    parts = split(line)
                    if length(parts) >= 11
                        site = parts[2]
                        haskey(site_z, site) || error("Site $site from $input_path was not found in the topo-corrected data object.")
                        parts[7] = _format_like_token(site_z[site], parts[7])
                        println(io_out, join(parts, " "))
                        continue
                    end
                end
                println(io_out, line)
            end
        end
    end
    println("ModEM data written to: $output_path")
    return output_path
end

function build_mesh_bundle(d, sx, sy, Tobs;
                           ρ_bg, dx_core, dy_core, nx_core, ny_core,
                           nx_pad, ny_pad, pad_factor, z_first, z_factor, depth_mult,
                           cov_value, topo_ctx = nothing)
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
    if USE_TOPO
        topo_ctx === nothing && error("Topography mode was requested without a valid projection context.")
        surface_local, datum_elev = _sample_surface_local(base, d, topo_ctx)
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
            n_air = AIR_LAYERS,
            air_factor = AIR_FACTOR,
            air_first = max(5.0, Float64(z_first) / AIR_FIRST_DIV))
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
              topo_name = surface_local === nothing ? "flat" : basename(TOPO_PATH),
              datum_elev = datum_elev,
              site_surface = site_surface,
              site_air = site_air)
end

function save_outputs(m, d_out, d_topo_diag)
    _write_start_model_ws(OUT_PATH, m; rotation = 0.0)
    write_covariance_file(OUT_COV_PATH, m.cov_mask; cov_x = m.cov_value, cov_y = m.cov_value, cov_z = m.cov_value)
    _write_topography_data_copy(DATA_PATH, OUT_DATA_PATH, d_out)
    _print_topography_data_report(d_topo_diag)
    _print_site_surface_report(m)
    _print_site_air_report(m)
    return nothing
end

#----- snap a value to the nearest in a range -----
snap(v, r) = collect(r)[argmin(abs.(collect(r) .- v))]

d_raw = load_data_modem(DATA_PATH)
d_raw.name = DATA_PATH
d_raw.origin = collect(Float64.(d_raw.origin))
d_raw.loc = Float64.(d_raw.loc)
topo_ctx = _local_tm_to_target_context(d_raw, TOPO_CRS)
d, d_topo_diag = _derive_topography_data(d_raw, topo_ctx)
sx = collect(Float64.(d.x))
sy = collect(Float64.(d.y))
Tobs = collect(Float64.(d.T))

spacing = nearest_neighbour_spacing(sx, sy)
ρ_off = filter(x -> isfinite(x) && x > 0, vec(d.ρ[:, [2, 3], :]))
ρ_bg_data = isempty(ρ_off) ? 100.0 : median(ρ_off)
ρ_bg0 = snap(ρ_bg_data, 10:10:20000)
dx_core0 = snap(spacing.median * CELL_WIDTH_FRAC, 50:25:10000)
dy_core0 = dx_core0
span_x0 = max(maximum(sx) - minimum(sx), 1.0)
span_y0 = max(maximum(sy) - minimum(sy), 1.0)
nx_core0 = max(1, ceil(Int, span_x0 / dx_core0))
ny_core0 = max(1, ceil(Int, span_y0 / dy_core0))
z_first0 = snap(skin_depth(Float64(ρ_bg0), minimum(Tobs)) / FIRST_LAYER_DIV, 5:5:2000)

@info @sprintf("Loaded %d sites, %d periods (%.3g–%.3g s); median spacing %.0f m; suggested background ρ = %.0f Ω·m; topography %s; topo-data output %s",
               length(sx), length(Tobs), minimum(Tobs), maximum(Tobs), spacing.median, ρ_bg_data,
               basename(TOPO_PATH), basename(OUT_DATA_PATH))

if MODE == "nogui"
    m = build_mesh_bundle(d, sx, sy, Tobs;
        ρ_bg = Float64(ρ_bg0), dx_core = Float64(dx_core0), dy_core = Float64(dy_core0),
        nx_core = nx_core0, ny_core = ny_core0,
        nx_pad = N_PAD, ny_pad = N_PAD,
        pad_factor = PAD_FACTOR, z_first = Float64(z_first0), z_factor = VERTICAL_FACTOR,
        depth_mult = DEPTH_MULT, cov_value = COV_VALUE0, topo_ctx = topo_ctx)
    @printf("grid %d×%d×%d (%d cells); core %d×%d @ %.0f×%.0f m; air %d; max sites/cell %d; occupied %.0f%%; pad %.0f km/side; depth %.0f km; cov %.2f (DEM-derived surface + DEM-derived site z)\n",
        m.nx, m.ny, m.nz, m.nx * m.ny * m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core,
        m.nz_air, m.maxper, 100 * m.occupied, m.pad_extent_km, m.depth_km, m.cov_value)
    save_outputs(m, d, d_topo_diag)
    @printf("wrote %s\n", OUT_PATH)
    @printf("wrote %s\n", OUT_COV_PATH)
    @printf("wrote %s\n", OUT_DATA_PATH)
    exit(0)
end

#----- rectangle outline for the core box -----
corner_points(yl, yh, xl, xh) = [Point2f(yl, xl), Point2f(yh, xl), Point2f(yh, xh),
                                 Point2f(yl, xh), Point2f(yl, xl)]

#----- check the mesh and advise what to fix -----
function suggestions(m)
    msgs = String[]
    m.maxper > 1 && push!(msgs,
        "• More than one site falls in some cells — decrease 'Core cell N–S (m)' / 'Core cell E–W (m)'.")
    m.pad_ok || push!(msgs,
        "• Lateral padding is too small — increase 'Pad cells N–S' / 'Pad cells E–W' or 'Pad factor'; " *
        "the boundary should sit at least one skin depth at the longest period (δ(Tmax)) from the core.")
    m.depth_km < m.δmax_km && push!(msgs,
        "• Model is too shallow — increase 'Depth × δ(Tmax)' so the mesh spans the depth of investigation.")
    !isempty(m.site_surface.impossible) && push!(msgs,
        "• Some sites still sit above the shallowest available earth cell — reduce the top air thickness or revise the topography/datum so the earth model starts above those site elevations.")
    m.site_air.ok || push!(msgs, m.site_air.text * " " * m.site_air.fix)
    ok = isempty(msgs)
    text = ok ? "✓ Mesh looks well-sized for ModEM inversion." : join(msgs, "    ")
    return (text = text, ok = ok)
end

const PARAM_HELP = Dict(
    "Core cell N–S (m)" =>
        "Core cell width in the N–S (X) direction, in metres. Use about half the " *
        "station spacing so no cell holds more than one site and structure between " *
        "stations is resolved. Smaller = finer resolution but more cells and slower.",
    "Core cell E–W (m)" =>
        "Core cell width in the E–W (Y) direction, in metres. Use about half the " *
        "station spacing so no cell holds more than one site and structure between " *
        "stations is resolved. Smaller = finer resolution but more cells and slower.",
    "Core cells N–S" =>
        "Number of core cells in the N–S (X) direction. Independent of 'Core cell N–S (m)': " *
        "the core spans count × cell width, centred on the stations.",
    "Core cells E–W" =>
        "Number of core cells in the E–W (Y) direction. Independent of 'Core cell E–W (m)': " *
        "the core spans count × cell width, centred on the stations.",
    "ρ background (Ω·m)" =>
        "Half-space resistivity that fills the starting model. It also sets the skin " *
        "depths that size the vertical mesh and padding. The default is suggested from " *
        "the median off-diagonal apparent resistivity of the data.",
    "Pad factor" =>
        "Geometric growth ratio of the padding cells (~1.4–1.6). A larger ratio expands " *
        "the grid quickly with few cells so the boundaries sit far from the sites.",
    "Pad cells N–S" =>
        "Number of padding cells added on each side in N–S (X). More cells push the " *
        "boundary farther out; it should reach ≥ one skin depth at the longest period " *
        "δ(Tmax) so boundary conditions do not affect the core.",
    "Pad cells E–W" =>
        "Number of padding cells added on each side in E–W (Y). More cells push the " *
        "boundary farther out; it should reach ≥ one skin depth at the longest period " *
        "δ(Tmax) so boundary conditions do not affect the core.",
    "First layer (m)" =>
        "Thickness of the top vertical cell. Set from the shortest-period skin depth " *
        "(≈ δ(Tmin) / 4) so the shallowest structure the data sees is resolved.",
    "Vertical factor" =>
        "Geometric growth of layer thickness with depth (~1.1–1.3). Larger = fewer " *
        "layers but coarser resolution at depth.",
    "Depth × δ(Tmax)" =>
        "Model bottom depth as a multiple of the longest-period skin depth δ(Tmax). " *
        "Use ≥ 1 so the model spans the depth of investigation; ~1.5 is typical.",
    "Covariance" =>
        "Recursive autoregression smoothing value written into the ModEM covariance file. " *
        "Use 0 to turn smoothing off, ~0.1 for weak smoothing, and ~0.3 for moderate smoothing.",
)

#----- draw the map view of the mesh -----
function draw2d!(ax, m)
    heatmap!(ax, m.y_edges_km, m.x_edges_km, fill(log10(m.ρ_bg), m.ny, m.nx);
             colormap = CMAP, colorrange = CRANGE)
    vlines!(ax, m.y_edges_km; color = GRID_COLOR, linewidth = GRID_WIDTH)
    hlines!(ax, m.x_edges_km; color = GRID_COLOR, linewidth = GRID_WIDTH)
    lines!(ax, corner_points(m.core_y0_km, m.core_y1_km, m.core_x0_km, m.core_x1_km);
           color = :black, linewidth = 2.5)
    scatter!(ax, sy ./ 1000, sx ./ 1000; color = SITE_COLOR, marker = :circle,
             markersize = site_size)
    return ax
end

fig = Figure(size = FIG_SIZE, figure_padding = (40, 14, 14, 14))

status = Observable("Outputs → $(basename(OUT_PATH)), $(basename(OUT_COV_PATH)), $(basename(OUT_DATA_PATH))")
help_obs = Observable("Click  ?  next to a parameter for an explanation.")
sug_obs = Observable("")
sug_color = Observable(RGBf(0.15, 0.5, 0.25))
show_full = Observable(true)
site_size = Observable(SITE_SIZE_FULL)

Colorbar(fig[1, 3]; colormap = CMAP, colorrange = CRANGE,
    label = "log₁₀ ρ (Ω·m)", width = 16)

left = GridLayout(fig[1, 1]; valign = :top)

infotext = join([
    @sprintf("datafile  >  %s", basename(DATA_PATH)),
    @sprintf("sites     >  %d", length(sx)),
    @sprintf("periods   >  %d   (%.3g – %.3g s)", length(Tobs), minimum(Tobs), maximum(Tobs)),
    @sprintf("spacing   >  %.1f km", spacing.median / 1000),
    @sprintf("rho(data) >  %.0f ohm-m", ρ_bg_data),
    @sprintf("topo      >  %s", basename(TOPO_PATH)),
    @sprintf("elev(src) >  DEM-derived, input site z ignored"),
    @sprintf("data(out) >  %s", basename(OUT_DATA_PATH)),
    @sprintf("cov(mode) >  DEM-derived surface + DEM-derived site z"),
    @sprintf("cov(out)  >  %s", basename(OUT_COV_PATH)),
], "\n")
Label(left[1, 1:3], infotext; font = "Consolas", fontsize = 13, halign = :left,
    justification = :left, color = :gray25, tellwidth = false)

next_row = Ref(1)

function add_param!(name, default; isint = false)
    r = (next_row[] += 1)
    Label(left[r, 1], name; halign = :right, fontsize = 13)
    s = isint ? string(Int(round(default))) : string(default)
    tb = Textbox(left[r, 2]; stored_string = s, validator = Float64, width = 84)
    ib = Button(left[r, 3]; label = "?", fontsize = 13, width = 26)
    on(ib.clicks) do _
        help_obs[] = get(PARAM_HELP, name, "")
    end
    return tb
end

tb_dx   = add_param!("Core cell N–S (m)", dx_core0; isint = true)
tb_dy   = add_param!("Core cell E–W (m)", dy_core0; isint = true)
tb_nxc  = add_param!("Core cells N–S", nx_core0; isint = true)
tb_nyc  = add_param!("Core cells E–W", ny_core0; isint = true)
tb_rho  = add_param!("ρ background (Ω·m)", ρ_bg0; isint = true)
tb_npx  = add_param!("Pad cells N–S", N_PAD; isint = true)
tb_npy  = add_param!("Pad cells E–W", N_PAD; isint = true)
tb_pf   = add_param!("Pad factor", PAD_FACTOR)
tb_zf   = add_param!("First layer (m)", z_first0; isint = true)
tb_zfac = add_param!("Vertical factor", VERTICAL_FACTOR)
tb_dm   = add_param!("Depth × δ(Tmax)", DEPTH_MULT)
tb_cov  = add_param!("Covariance", COV_VALUE0)

updatebtn = Button(left[(next_row[] += 1), 1:3]; label = "Update mesh", fontsize = 14)

viewgrid = GridLayout(left[(next_row[] += 1), 1:3])
corefullbtn = Button(viewgrid[1, 1]; label = "Show core", fontsize = 14, tellwidth = false)
resetbtn    = Button(viewgrid[1, 2]; label = "Reset zoom", fontsize = 14, tellwidth = false)
colsize!(viewgrid, 1, Relative(0.5)); colsize!(viewgrid, 2, Relative(0.5))
colgap!(viewgrid, 8)

actiongrid = GridLayout(left[(next_row[] += 1), 1:3])
savebtn   = Button(actiongrid[1, 1]; label = "Save model + cov + data", fontsize = 14, tellwidth = false)
exportbtn = Button(actiongrid[1, 2]; label = "Export figure", fontsize = 14, tellwidth = false)
colsize!(actiongrid, 1, Relative(0.5)); colsize!(actiongrid, 2, Relative(0.5))
colgap!(actiongrid, 8)

rowgap!(left, 7)
colgap!(left, 8)
colsize!(left, 2, Fixed(86))
colsize!(left, 3, Fixed(26))

Label(fig[2, 1:3], sug_obs; fontsize = 13, halign = :left, justification = :left,
    color = sug_color, word_wrap = true, tellwidth = false)
Label(fig[3, 1:3], help_obs; fontsize = 12.5, halign = :left, justification = :left,
    color = :gray30, word_wrap = true, tellwidth = false)
Label(fig[4, 1:3], status; fontsize = 11, halign = :left, color = :gray45,
    word_wrap = true, tellwidth = false)

colsize!(fig.layout, 1, Fixed(290))
colsize!(fig.layout, 3, Fixed(70))
rowsize!(fig.layout, 1, Relative(0.82))
for r in 2:4
    rowsize!(fig.layout, r, Auto())
end
rowgap!(fig.layout, 6)

tbget(tb, dflt) = begin
    s = tb.displayed_string[]
    (s === nothing || isempty(strip(s))) && return dflt
    v = tryparse(Float64, s)
    v === nothing ? dflt : v
end

function build_mesh()
    build_mesh_bundle(d, sx, sy, Tobs;
        ρ_bg = tbget(tb_rho, Float64(ρ_bg0)),
        dx_core = tbget(tb_dx, Float64(dx_core0)),
        dy_core = tbget(tb_dy, Float64(dy_core0)),
        nx_core = Int(round(tbget(tb_nxc, Float64(nx_core0)))),
        ny_core = Int(round(tbget(tb_nyc, Float64(ny_core0)))),
        nx_pad = Int(round(tbget(tb_npx, Float64(N_PAD)))),
        ny_pad = Int(round(tbget(tb_npy, Float64(N_PAD)))),
        pad_factor = tbget(tb_pf, PAD_FACTOR),
        z_first = tbget(tb_zf, Float64(z_first0)),
        z_factor = tbget(tb_zfac, VERTICAL_FACTOR),
        depth_mult = tbget(tb_dm, DEPTH_MULT),
        cov_value = tbget(tb_cov, COV_VALUE0),
        topo_ctx = topo_ctx)
end

mesh = Observable(build_mesh())
current_ax = Ref{Any}(nothing)

function apply_limits!(ax, m)
    if show_full[]
        xlims!(ax, m.full_y0_km, m.full_y1_km)
        ylims!(ax, m.full_x0_km, m.full_x1_km)
    else
        padx = 0.05 * (m.core_y1_km - m.core_y0_km)
        pady = 0.05 * (m.core_x1_km - m.core_x0_km)
        xlims!(ax, m.core_y0_km - padx, m.core_y1_km + padx)
        ylims!(ax, m.core_x0_km - pady, m.core_x1_km + pady)
    end
    return nothing
end

function refresh!(m)
    current_ax[] === nothing || delete!(current_ax[])
    ax = Axis(fig[1, 2]; xlabel = "Y East (km)", ylabel = "X North (km)",
        aspect = DataAspect(), halign = :right,
        xgridvisible = false, ygridvisible = false,
        xlabelfont = :bold, ylabelfont = :bold,
        xticklabelcolor = :gray55, yticklabelcolor = :gray55)
    draw2d!(ax, m)
    current_ax[] = ax
    apply_limits!(ax, m)
    s = suggestions(m)
    topo_info = m.topo_active ? @sprintf("topo %s · datum %.1f m · air %d", m.topo_name, m.datum_elev, m.nz_air) : "topo flat"
    site_info = m.site_air.ok ?
        (m.site_surface.nadjusted == 0 ? "all sites below surface" : @sprintf("all sites below surface (%d auto-adjusted column(s))", m.site_surface.nadjusted)) :
        @sprintf("%d site(s) in air", m.site_air.nbad)
    gridinfo = @sprintf("%d × %d × %d cells · core %d × %d (%.0f × %.0f m) · max sites/cell %d · depth %.0f km · cov %.2f (DEM-derived surface + DEM-derived site z) · %s · %s",
        m.nx, m.ny, m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core, m.maxper, m.depth_km, m.cov_value, topo_info, site_info)
    sug_obs[] = gridinfo * "\n" * s.text
    sug_color[] = s.ok ? RGBf(0.15, 0.5, 0.25) : RGBf(0.75, 0.2, 0.15)
    return nothing
end

on(mesh) do m
    refresh!(m)
end
refresh!(mesh[])

on(updatebtn.clicks) do _
    mesh[] = build_mesh()
end

for tb in (tb_dx, tb_dy, tb_nxc, tb_nyc, tb_rho, tb_npx, tb_npy, tb_pf, tb_zf, tb_zfac, tb_dm, tb_cov)
    on(tb.stored_string) do _
        mesh[] = build_mesh()
    end
end

on(corefullbtn.clicks) do _
    show_full[] = !show_full[]
    corefullbtn.label[] = show_full[] ? "Show core" : "Show full"
    site_size[] = show_full[] ? SITE_SIZE_FULL : SITE_SIZE_CORE
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(resetbtn.clicks) do _
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(savebtn.clicks) do _
    m = mesh[]
    save_outputs(m, d, d_topo_diag)
    air_status = m.site_air.ok ?
        (m.site_surface.nadjusted == 0 ? "all sites below surface" : @sprintf("all sites below surface after %d auto-adjusted column(s)", m.site_surface.nadjusted)) :
        @sprintf("WARNING: %d site(s) in air; see guidance above", m.site_air.nbad)
    status[] = @sprintf("Saved %d×%d×%d model → %s | covariance %.2f (DEM-derived surface + DEM-derived site z) → %s | topo-data → %s | %s", m.nx, m.ny, m.nz, OUT_PATH, m.cov_value, OUT_COV_PATH, OUT_DATA_PATH, air_status)
    @info status[]
end

on(exportbtn.clicks) do _
    m = mesh[]
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    efig = Figure(size = (1200, 1000))
    eax = Axis(efig[1, 1]; xlabel = "Y East (km)", ylabel = "X North (km)",
        aspect = DataAspect(), halign = :right, xgridvisible = false, ygridvisible = false)
    draw2d!(eax, m)
    apply_limits!(eax, m)
    Colorbar(efig[1, 2]; colormap = CMAP, colorrange = CRANGE,
        label = "log₁₀ ρ (Ω·m)", width = 16)
    fn = joinpath(dirname(OUT_PATH), "mesh_map_$ts.png")
    save(fn, efig; px_per_unit = 3)
    status[] = "Saved figure → $fn"
    @info status[]
end

screen = display(fig)
if !isinteractive()
    wait(screen)
end

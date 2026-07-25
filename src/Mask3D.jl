# bathymetry/topography extraction and water/air masks for WS3D models
# water cells are detected from log10 resistivity below a threshold, air cells
# from the NaN tag the loader applies above 10^15 Ω·m; both files store the base
# of the tagged column run at physical cell centers in the model frame (x, y and
# z carry the model origin), so they can be re-applied to any grid on the same
# coordinate frame and plotted directly on a model section

using Dates
using Printf

# base of the contiguous run of tagged cells at the top of each column, as depth
# below the model top and as z in the model frame
function _tagged_column_base(m::WS3DModel, istagged)
    zbot = cumsum(m.dz)
    xs = Float64[]; ys = Float64[]; ds = Float64[]; zs = Float64[]
    n_stray = 0
    for j in 1:m.ny, i in 1:m.nx
        klast = 0
        for k in 1:m.nz
            istagged(m.A[i, j, k]) || break
            klast = k
        end
        n_stray += count(k -> istagged(m.A[i, j, k]), klast+1:m.nz)
        if klast > 0
            push!(xs, m.cx[i]); push!(ys, m.cy[j])
            push!(ds, zbot[klast]); push!(zs, m.z[klast+1])
        end
    end
    return (x=xs, y=ys, depth=ds, z=zs), n_stray
end

function _write_surface_table(path::AbstractString, cols, title::AbstractString,
                              base_label::AbstractString, meta::AbstractString,
                              origin::AbstractVector{<:Real})
    open(path, "w") do io
        println(io, "# MTGeophysics ", title, " — ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "# ", meta, "  n_columns=", length(cols.depth))
        @printf(io, "# origin_x_m=%.3f  origin_y_m=%.3f  origin_z_m=%.3f\n",
                origin[1], origin[2], origin[3])
        println(io, "# x_center_m  y_center_m  ", base_label, "  base_z_m")
        hasz = hasproperty(cols, :z)
        for t in eachindex(cols.depth)
            @printf(io, "%15.3f %15.3f %12.3f %12.3f\n",
                    cols.x[t], cols.y[t], cols.depth[t], hasz ? cols.z[t] : NaN)
        end
    end
    return path
end

function _read_surface_table(path::AbstractString)
    xs = Float64[]; ys = Float64[]; ds = Float64[]; zs = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        p = split(s)
        push!(xs, parse(Float64, p[1]))
        push!(ys, parse(Float64, p[2]))
        push!(ds, parse(Float64, p[3]))
        push!(zs, length(p) >= 4 ? parse(Float64, p[4]) : NaN)
    end
    return (x=xs, y=ys, depth=ds, z=zs)
end

# depth is measured from the model top while m.cz carries origin[3], so the cell
# test uses cumsum(dz); comparing centers leaves half a cell of slack, which the
# 3-decimal file format needs
function _mask_above_base(m::WS3DModel, cols)
    mask = falses(m.nx, m.ny, m.nz)
    zcen = cumsum(m.dz) .- m.dz ./ 2
    for t in eachindex(cols.depth)
        i = argmin(abs.(m.cx .- cols.x[t]))
        j = argmin(abs.(m.cy .- cols.y[t]))
        for k in 1:m.nz
            zcen[k] < cols.depth[t] || break
            mask[i, j, k] = true
        end
    end
    return mask
end

"""
    extract_bathymetry(m; water_log10) -> (x, y, depth)

Water-bottom depth (m, positive down) per wet column of a WS3D model. A column
is wet where its surface cells have `log10 ρ < water_log10`; depth is the bottom
edge of the deepest contiguous water cell. Non-contiguous water below the
seafloor triggers a warning and is ignored. Returns coordinate vectors of the
wet cell-center columns and their depths.
"""
function extract_bathymetry(m::WS3DModel; water_log10::Float64)
    cols, n_stray = _tagged_column_base(m, v -> v < water_log10)
    n_stray > 0 && @warn "extract_bathymetry: $n_stray water-like cells below the seafloor ignored"
    return cols
end

"""
    write_bathymetry(path, bathy; water_log10=NaN, origin=[0,0,0])

Write a bathymetry table (`x_center_m  y_center_m  water_depth_m  base_z_m`, one
row per wet column) as produced by [`extract_bathymetry`](@ref). `x`/`y` are cell
centers and `base_z_m` the seafloor in the model frame, both carrying `origin`,
so the surface plots directly on a model section.
"""
write_bathymetry(path::AbstractString, bathy; water_log10::Float64=NaN,
                 origin::AbstractVector{<:Real}=[0.0, 0.0, 0.0]) =
    _write_surface_table(path, bathy, "bathymetry", "water_depth_m",
                         string("water_log10=", water_log10), origin)

"""
    read_bathymetry(path) -> (x, y, depth, z)

Read a bathymetry table written by [`write_bathymetry`](@ref). `z` is `NaN` for
tables written before the column existed.
"""
read_bathymetry(path::AbstractString) = _read_surface_table(path)

"""
    extract_topography(m) -> (x, y, depth, z)

Ground-surface depth (m, positive down) per column of a WS3D model that carries
topographic air. Air is the `NaN` tag applied by [`load_ws3d_model`](@ref) above
10^15 Ω·m; depth is the bottom edge of the contiguous air run at the top of the
column, `z` the same interface in the model frame. Air below the ground surface
triggers a warning and is ignored.
"""
function extract_topography(m::WS3DModel)
    cols, n_stray = _tagged_column_base(m, isnan)
    n_stray > 0 && @warn "extract_topography: $n_stray air cells below the ground surface ignored"
    return cols
end

"""
    write_topography(path, topo; origin=[0,0,0])

Write a topography table (`x_center_m  y_center_m  ground_depth_m  base_z_m`,
one row per column with air) as produced by [`extract_topography`](@ref). Same
frame convention as [`write_bathymetry`](@ref).
"""
write_topography(path::AbstractString, topo;
                 origin::AbstractVector{<:Real}=[0.0, 0.0, 0.0]) =
    _write_surface_table(path, topo, "topography", "ground_depth_m",
                         "air_tag=NaN(log10rho>15)", origin)

"""
    read_topography(path) -> (x, y, depth, z)

Read a topography table written by [`write_topography`](@ref).
"""
read_topography(path::AbstractString) = _read_surface_table(path)

"""
    water_mask_from_bathymetry(m, bathy) -> BitArray{3}

Full-grid mask of water cells: for each wet column (matched to the nearest
cell-center column of `m`), cells whose centers lie above the water-bottom
depth are water.
"""
water_mask_from_bathymetry(m::WS3DModel, bathy) = _mask_above_base(m, bathy)

"""
    air_mask_from_topography(m, topo) -> BitArray{3}

Full-grid mask of air cells: for each column with air (matched to the nearest
cell-center column of `m`), cells whose centers lie above the ground surface.
"""
air_mask_from_topography(m::WS3DModel, topo) = _mask_above_base(m, topo)

"""
    water_mask_from_model(m; water_log10) -> BitArray{3}

Full-grid mask of cells with `log10 ρ < water_log10`.
"""
water_mask_from_model(m::WS3DModel; water_log10::Float64) = BitArray(m.A .< water_log10)

"""
    air_mask_from_model(m) -> BitArray{3}

Full-grid mask of topographic air, i.e. the cells the loader tagged `NaN`.
"""
air_mask_from_model(m::WS3DModel) = BitArray(isnan.(m.A))

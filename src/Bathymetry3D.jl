# bathymetry extraction and water masks for WS3D models
# water cells are detected from log10 resistivity below a threshold; the
# bathymetry file stores water-bottom depth per wet column at physical cell
# centers, so it can be re-applied to any grid on the same coordinate frame

using Dates
using Printf

"""
    extract_bathymetry(m; water_log10) -> (x, y, depth)

Water-bottom depth (m, positive down) per wet column of a WS3D model. A column
is wet where its surface cells have `log10 ρ < water_log10`; depth is the bottom
edge of the deepest contiguous water cell. Non-contiguous water below the
seafloor triggers a warning and is ignored. Returns coordinate vectors of the
wet cell-center columns and their depths.
"""
function extract_bathymetry(m::WS3DModel; water_log10::Float64)
    zbot = cumsum(m.dz)
    xs = Float64[]; ys = Float64[]; ds = Float64[]
    n_stray = 0
    for j in 1:m.ny, i in 1:m.nx
        klast = 0
        for k in 1:m.nz
            m.A[i, j, k] < water_log10 || break
            klast = k
        end
        n_stray += count(k -> m.A[i, j, k] < water_log10, klast+1:m.nz)
        if klast > 0
            push!(xs, m.cx[i]); push!(ys, m.cy[j]); push!(ds, zbot[klast])
        end
    end
    n_stray > 0 && @warn "extract_bathymetry: $n_stray water-like cells below the seafloor ignored"
    return (x=xs, y=ys, depth=ds)
end

"""
    write_bathymetry(path, bathy; water_log10=NaN)

Write a bathymetry table (`x_center_m  y_center_m  water_depth_m`, one row per
wet column) as produced by [`extract_bathymetry`](@ref).
"""
function write_bathymetry(path::AbstractString, bathy; water_log10::Float64=NaN)
    open(path, "w") do io
        println(io, "# MTGeophysics bathymetry — ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "# water_log10=", water_log10, "  n_wet=", length(bathy.depth))
        println(io, "# x_center_m  y_center_m  water_depth_m")
        for t in eachindex(bathy.depth)
            @printf(io, "%15.3f %15.3f %12.3f\n", bathy.x[t], bathy.y[t], bathy.depth[t])
        end
    end
    return path
end

"""
    read_bathymetry(path) -> (x, y, depth)

Read a bathymetry table written by [`write_bathymetry`](@ref).
"""
function read_bathymetry(path::AbstractString)
    xs = Float64[]; ys = Float64[]; ds = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        p = split(s)
        push!(xs, parse(Float64, p[1]))
        push!(ys, parse(Float64, p[2]))
        push!(ds, parse(Float64, p[3]))
    end
    return (x=xs, y=ys, depth=ds)
end

"""
    water_mask_from_bathymetry(m, bathy) -> BitArray{3}

Full-grid mask of water cells: for each wet column (matched to the nearest
cell-center column of `m`), cells whose centers lie above the water-bottom
depth are water.
"""
function water_mask_from_bathymetry(m::WS3DModel, bathy)
    mask = falses(m.nx, m.ny, m.nz)
    zc = abs.(m.cz)
    for t in eachindex(bathy.depth)
        i = argmin(abs.(m.cx .- bathy.x[t]))
        j = argmin(abs.(m.cy .- bathy.y[t]))
        for k in 1:m.nz
            zc[k] < bathy.depth[t] || break
            mask[i, j, k] = true
        end
    end
    return mask
end

"""
    water_mask_from_model(m; water_log10) -> BitArray{3}

Full-grid mask of cells with `log10 ρ < water_log10`.
"""
water_mask_from_model(m::WS3DModel; water_log10::Float64) = BitArray(m.A .< water_log10)

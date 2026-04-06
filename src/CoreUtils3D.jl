# Headless 3D model core/padding detection utilities.
# Author: @pankajkmishra
# These functions identify unpadded core regions and depth ranges without
# requiring any visualization backend. Safe for clusters and CI.

"""
    edges_from_centers(c::AbstractVector)

Compute cell-edge coordinates from cell centers by midpoint extrapolation.
"""
edges_from_centers(c::AbstractVector) = begin
    n = length(c)
    n < 2 && return [c[1] - 0.5; c[1] + 0.5]
    mids = (c[1:end-1] .+ c[2:end]) ./ 2
    first_edge = c[1] - (c[2] - c[1]) / 2
    last_edge  = c[end] + (c[end] - c[end-1]) / 2
    vcat(first_edge, mids, last_edge)
end

"""
    core_indices(c::AbstractVector; tol::Real=0.2)

Return the index range of the largest contiguous block of cells whose widths
are within `tol * w_ref` of the interior reference width (median over the
central 60% of the mesh).
"""
function core_indices(c::AbstractVector; tol::Real = 0.2)
    e = edges_from_centers(c)
    w = abs.(diff(e))
    n = length(w)
    n <= 4 && return 1:n
    s = max(1, round(Int, 0.2n))
    t = min(n, round(Int, 0.8n))
    w_ref = median(@view w[(s+1):t])
    is_core = abs.(w .- w_ref) .<= tol * w_ref
    best_i, best_len = 1, 0
    i = 1
    while i <= n
        if is_core[i]
            j = i
            while j <= n && is_core[j]; j += 1; end
            if (j - i) > best_len
                best_i, best_len = i, (j - i)
            end
            i = j
        else
            i += 1
        end
    end
    best_len > 0 ? (best_i:(best_i + best_len - 1)) : (1:n)
end

"""
    z_indices_for_max_depth(zc::AbstractVector, max_depth::Real)

Return `1:k` where cumulative depth best matches `max_depth`.
"""
function z_indices_for_max_depth(zc::AbstractVector, max_depth::Real)
    e = edges_from_centers(zc)
    dz = abs.(diff(e))
    cum = cumsum(dz)
    if max_depth <= cum[1]
        return 1:1
    elseif max_depth >= cum[end]
        return 1:length(zc)
    else
        k = findfirst(>=(max_depth), cum)::Int
        if k > 1 && (max_depth - cum[k-1]) < (cum[k] - max_depth)
            return 1:(k-1)
        else
            return 1:k
        end
    end
end

"""
    lateral_core_ranges(m; tol=0.2)

Return `(ix, iy)` core index ranges for any model with `.cx` and `.cy` fields.
"""
function lateral_core_ranges(m; tol::Real = 0.2)
    ix = core_indices(m.cx; tol=tol)
    iy = core_indices(m.cy; tol=tol)
    return ix, iy
end

"""
    core_view(m; tol=0.2)

Return `(Rview, ix, iy)` where `Rview = @view m.A[ix, iy, :]`.
"""
function core_view(m; tol::Real = 0.2)
    ix, iy = lateral_core_ranges(m; tol=tol)
    return (@view m.A[ix, iy, :]), ix, iy
end

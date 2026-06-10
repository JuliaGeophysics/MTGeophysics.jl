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
    _mode_width(w::AbstractVector)

Most frequently occurring cell width in `w`. Widths are rounded to 6 significant
digits before tallying, since cell-edge reconstruction blends only the 1–2 cells
straddling a core/padding boundary while every interior core cell shares an
identical width. Ties are broken toward the smaller width (finer cells).
"""
function _mode_width(w::AbstractVector)
    counts = Dict{Float64,Int}()
    @inbounds for v in w
        key = round(float(v); sigdigits = 6)
        counts[key] = get(counts, key, 0) + 1
    end
    best_w, best_c = first(keys(counts)), 0
    for (k, c) in counts
        if c > best_c || (c == best_c && k < best_w)
            best_w, best_c = k, c
        end
    end
    return best_w
end

"""
    core_indices(c::AbstractVector; tol::Real=0.2)

Return the index range of the largest contiguous block of cells whose widths
are within `tol * w_ref` of the reference core width `w_ref`.

The core (survey) region is meshed with uniform fine cells, so `w_ref` is taken
as the *most common* cell width (mode). This recovers the fine-cell plateau no
matter where it sits or how small a fraction of the axis it occupies — unlike a
central-window median, which is polluted by padding when the core is narrow.
"""
function core_indices(c::AbstractVector; tol::Real = 0.2)
    e = edges_from_centers(c)
    w = abs.(diff(e))
    n = length(w)
    n <= 4 && return 1:n
    w_ref = _mode_width(w)
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

# Shared plotting utilities for model viewers.
# Author: @pankajkmishra
# This file contains grid/index/depth/color helper functions reused by 2D and 3D visualization scripts.
# Keeping these helpers here avoids duplicated plotting logic in examples.

function edges_from_centers(c::AbstractVector)
    n = length(c)
    n < 2 && return [c[1] - 0.5; c[1] + 0.5]
    mids = (c[1:end-1] .+ c[2:end]) ./ 2
    first_edge = c[1] - (c[2] - c[1]) / 2
    last_edge  = c[end] + (c[end] - c[end-1]) / 2
    vcat(first_edge, mids, last_edge)
end

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

function compute_colorrange(R::AbstractArray;
        resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}, AbstractVector{<:Real}} = nothing)
    if isnothing(resistivity_range)
        vals = R[isfinite.(R)]
        isempty(vals) && (vals = [0.0, 1.0])
        qlo, qhi = quantile(vec(vals), (0.02, 0.98))
        cmin, cmax = min(qlo, qhi), max(qlo, qhi)
        if cmin == cmax
            ϵ = max(1e-12, 1e-6 * abs(cmin))
            cmin -= ϵ; cmax += ϵ
        end
    else
        lo, hi = Float64(resistivity_range[1]), Float64(resistivity_range[2])
        (!isfinite(lo) || !isfinite(hi)) && error("resistivity_range must contain finite numbers")
        lo > hi && ((lo, hi) = (hi, lo))
        if lo == hi
            ϵ = max(1e-12, 1e-6 * abs(lo))
            lo -= ϵ; hi += ϵ
        end
        cmin, cmax = lo, hi
    end
    return cmin, cmax
end

function prepare_model_arrays(M;
        log10scale::Bool = true,
        withPadding::Bool = true,
        max_depth::Union{Nothing, Real} = nothing,
        pad_tol::Real = 0.2)
    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : copy(M.A)

    ix = withPadding ? (1:length(x_all)) : core_indices(x_all; tol = pad_tol)
    iy = withPadding ? (1:length(y_all)) : core_indices(y_all; tol = pad_tol)
    kz = isnothing(max_depth) ? (1:length(z_all)) : z_indices_for_max_depth(z_all, float(max_depth))

    x = x_all[ix]
    y = y_all[iy]
    z = z_all[kz]
    R = A_all[ix, iy, kz]

    return x, y, z, R, ix, iy, kz
end

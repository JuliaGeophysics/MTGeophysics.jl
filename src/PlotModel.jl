# Shared plotting utilities for model viewers.
# Author: @pankajkmishra
# This file contains color/model-array helper functions reused by 2D and 3D visualization scripts.
# Core grid/index/depth utilities (edges_from_centers, core_indices, z_indices_for_max_depth)
# are defined in CoreUtils3D.jl and always available, even without GLMakie.

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

using GLMakie
using Statistics

# --- helpers ---------------------------------------------------------------

edges_from_centers(c::AbstractVector) = begin
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

# --- viewer ----------------------------------------------------------------

function gl_modem_viewer(
    M;
    log10scale::Bool = true,
    cmap = :viridis,
    figsize = (1200, 900),
    withPadding::Bool = true,
    max_depth::Union{Nothing, Real} = nothing,
    pad_tol::Real = 0.2,
    resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}, AbstractVector{<:Real}} = nothing # NEW
)
    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : M.A

    if withPadding
        ix = 1:length(x_all)
        iy = 1:length(y_all)
    else
        ix = core_indices(x_all; tol = pad_tol)
        iy = core_indices(y_all; tol = pad_tol)
    end

    x = x_all[ix]
    y = y_all[iy]
    R = @view A_all[ix, iy, :]

    if isnothing(max_depth)
        kz = 1:length(z_all)
    else
        kz = z_indices_for_max_depth(z_all, float(max_depth))
    end
    z = z_all[kz]
    R = @view R[:, :, kz]

    # ---- Color scaling -----------------------------------------------------
    if isnothing(resistivity_range)
        vals = R[isfinite.(R)]
        if isempty(vals)
            vals = [0.0, 1.0]
        end
        qlo, qhi = quantile(vec(vals), (0.02, 0.98))
        cmin, cmax = min(qlo, qhi), max(qlo, qhi)
        if cmin == cmax
            ϵ = max(1e-12, 1e-6 * abs(cmin))
            cmin -= ϵ; cmax += ϵ
        end
    else
        lo, hi = resistivity_range[1], resistivity_range[2]
        if !isfinite(lo) || !isfinite(hi)
            error("resistivity_range must contain finite numbers")
        end
        if lo > hi
            lo, hi = hi, lo
        end
        if lo == hi
            ϵ = max(1e-12, 1e-6 * abs(lo))
            lo -= ϵ; hi += ϵ
        end
        cmin, cmax = lo, hi
    end
    # -----------------------------------------------------------------------

    fig = Figure(resolution = figsize)
    ax  = LScene(fig[1, 1], show_axis = false)

    plt = volumeslices!(ax, x, y, z, R; colormap = cmap, colorrange = (cmin, cmax), bbox_visible = true)


    controls = fig[2, 1] = GridLayout()

    sl_yz = Slider(controls[1, 1], range = 1:length(x), startvalue = round(Int, length(x) / 2))
    sl_xz = Slider(controls[2, 1], range = 1:length(y), startvalue = round(Int, length(y) / 2))
    sl_xy = Slider(controls[3, 1], range = 1:length(z), startvalue = round(Int, length(z) / 2))

    Label(controls[1, 3], "YZ")
    Label(controls[2, 3], "XZ")
    Label(controls[3, 3], "XY")

    tog_yz = Toggle(controls[1, 2], active = true)
    tog_xz = Toggle(controls[2, 2], active = true)
    tog_xy = Toggle(controls[3, 2], active = true)


    on(sl_yz.value) do v; plt[:update_yz][](v) end
    on(sl_xz.value) do v; plt[:update_xz][](v) end
    on(sl_xy.value) do v; plt[:update_xy][](v) end

    h_yz = plt[:heatmap_yz][]
    h_xz = plt[:heatmap_xz][]
    h_xy = plt[:heatmap_xy][]
    on(tog_yz.active) do a; h_yz.visible = a end
    on(tog_xz.active) do a; h_xz.visible = a end
    on(tog_xy.active) do a; h_xy.visible = a end

    Colorbar(fig[1, 2], h_xy, label = log10scale ? "log10 ρ (Ω·m)" : "ρ (Ω·m)")

    return fig, (ax = ax, plt = plt,
                 sliders = (sl_yz, sl_xz, sl_xy),
                 toggles = (tog_yz, tog_xz, tog_xy),
                 indices = (ix = ix, iy = iy, kz = kz),
                 colorrange = (cmin = cmin, cmax = cmax))
end

# Examples:
# fig, parts = gl_modem_viewer(M; resistivity_range = (1.0, 4.0))                # log10 scale by default
# fig, parts = gl_modem_viewer(M; log10scale=false, resistivity_range=(10, 1e4)) # linear Ω·m
# fig, parts = gl_modem_viewer(M; withPadding=false, max_depth=50000.0, resistivity_range=[1.5, 3.5])
# display(fig)

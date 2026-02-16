# Example: interactive 3D volume-slice model viewer.
# Author: @pankajkmishra
# This script visualizes XY/XZ/YZ slices with controls for depth, padding, and export.
# Use it to explore overall 3D model geometry and resistivity trends.

using Pkg
Pkg.activate(dirname(@__DIR__))

using GLMakie
using Statistics
using Dates

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))

model_file = joinpath(@__DIR__, "I_NLCG_140.rho")

log10_scale = true

colormap = Reverse(:turbo)

resistivity_range = (1.0, 4.0)

max_depth = 50000.0

show_padding = false

pad_tolerance = 0.2

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

function modem_3d_viewer(
    M;
    log10scale::Bool = true,
    cmap = Reverse(:turbo),
    figsize = (1200, 800),
    withPadding::Bool = true,
    max_depth::Union{Nothing, Real} = nothing,
    pad_tol::Real = 0.2,
    resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}} = nothing
)

    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : M.A

    ix_full = 1:length(x_all)
    iy_full = 1:length(y_all)
    ix_core = core_indices(x_all; tol = pad_tol)
    iy_core = core_indices(y_all; tol = pad_tol)

    if withPadding
        ix = ix_full
        iy = iy_full
    else
        ix = ix_core
        iy = iy_core
    end

    x = x_all[ix]
    y = y_all[iy]
    R = @view A_all[ix, iy, :]

    if isnothing(max_depth)
        kz = 1:length(z_all)
    else
        kz = z_indices_for_max_depth(z_all, float(max_depth))
    end

    z = -z_all[kz]
    R = A_all[ix, iy, kz]

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
        if lo > hi
            lo, hi = hi, lo
        end
        if lo == hi
            ϵ = max(1e-12, 1e-6 * abs(lo))
            lo -= ϵ; hi += ϵ
        end
        cmin, cmax = lo, hi
    end

    current_colormap = Observable(cmap)
    show_full_model = Observable(withPadding)
    current_x = Observable(copy(x))
    current_y = Observable(copy(y))
    current_z = Observable(copy(z))
    current_R = Observable(copy(R))

    current_plt = Ref{Any}(nothing)
    current_heatmaps = Ref{NamedTuple}((yz=nothing, xz=nothing, xy=nothing))

    fig = Figure(size = figsize)

    left_panel = fig[1, 1] = GridLayout(valign = :top)
    controls = left_panel[1, 1] = GridLayout()

    ax = LScene(fig[1, 2], show_axis = false)

    cb_label = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"

    colsize!(fig.layout, 1, Fixed(180))
    colsize!(fig.layout, 2, Auto())

    plt = volumeslices!(ax, current_x[], current_y[], current_z[], current_R[]; 
                        colormap = current_colormap[], 
                        colorrange = (cmin, cmax), 
                        bbox_visible = true)
    current_plt[] = plt
    current_heatmaps[] = (yz=plt[:heatmap_yz][], xz=plt[:heatmap_xz][], xy=plt[:heatmap_xy][])

    cam3d!(ax.scene, projectiontype = Makie.Perspective)

    cb = Colorbar(fig[1, 3], current_heatmaps[].xy, label = cb_label, width = 15)
    colsize!(fig.layout, 3, Fixed(60))

    row = 1

    Label(controls[row, 1], "X:", halign = :right, fontsize = 12)
    sl_yz = Slider(controls[row, 2], range = 1:length(x), startvalue = round(Int, length(x) / 2), width = 120)
    tog_yz = Toggle(controls[row, 3], active = true)
    row += 1
    lbl_yz = Label(controls[row, 1:3], "$(round(x[sl_yz.value[]], digits=0)) m", fontsize = 11, color = :gray30)
    row += 1

    Label(controls[row, 1], "Y:", halign = :right, fontsize = 12)
    sl_xz = Slider(controls[row, 2], range = 1:length(y), startvalue = round(Int, length(y) / 2), width = 120)
    tog_xz = Toggle(controls[row, 3], active = true)
    row += 1
    lbl_xz = Label(controls[row, 1:3], "$(round(y[sl_xz.value[]], digits=0)) m", fontsize = 11, color = :gray30)
    row += 1

    Label(controls[row, 1], "Z:", halign = :right, fontsize = 12)
    sl_xy = Slider(controls[row, 2], range = 1:length(z), startvalue = round(Int, length(z) / 2), width = 120)
    tog_xy = Toggle(controls[row, 3], active = true)
    row += 1
    lbl_xy = Label(controls[row, 1:3], "Depth: $(round(-z[sl_xy.value[]], digits=0)) m", fontsize = 11, color = :gray30)
    row += 1

    btn_toggle = Button(controls[row, 1:3], label = show_full_model[] ? "Show Core" : "Show Full", fontsize = 11)
    row += 1

    btn_reset = Button(controls[row, 1:3], label = "Reset View", fontsize = 11)
    row += 1

    btn_export = Button(controls[row, 1:3], label = "Export", fontsize = 11)
    row += 1

    rowgap!(controls, 2)
    colgap!(controls, 2)

    on(sl_yz.value) do v
        if current_plt[] !== nothing
            current_plt[][:update_yz][](v)
        end
        lbl_yz.text[] = "$(round(current_x[][v], digits=0)) m"
    end
    on(sl_xz.value) do v
        if current_plt[] !== nothing
            current_plt[][:update_xz][](v)
        end
        lbl_xz.text[] = "$(round(current_y[][v], digits=0)) m"
    end
    on(sl_xy.value) do v
        if current_plt[] !== nothing
            current_plt[][:update_xy][](v)
        end
        lbl_xy.text[] = "Depth: $(round(-current_z[][v], digits=0)) m"
    end

    on(tog_yz.active) do a
        if current_heatmaps[].yz !== nothing
            current_heatmaps[].yz.visible = a
        end
    end
    on(tog_xz.active) do a
        if current_heatmaps[].xz !== nothing
            current_heatmaps[].xz.visible = a
        end
    end
    on(tog_xy.active) do a
        if current_heatmaps[].xy !== nothing
            current_heatmaps[].xy.visible = a
        end
    end

    function rebuild_plot!(new_x, new_y, new_z, new_R, new_cmap)

        for plot in copy(ax.scene.plots)
            delete!(ax.scene, plot)
        end

        new_plt = volumeslices!(ax, new_x, new_y, new_z, new_R; 
                                colormap = new_cmap, 
                                colorrange = (cmin, cmax), 
                                bbox_visible = true)

        current_plt[] = new_plt
        current_heatmaps[] = (yz=new_plt[:heatmap_yz][], xz=new_plt[:heatmap_xz][], xy=new_plt[:heatmap_xy][])

        sl_yz.range[] = 1:length(new_x)
        sl_xz.range[] = 1:length(new_y)
        sl_xy.range[] = 1:length(new_z)

        set_close_to!(sl_yz, round(Int, length(new_x) / 2))
        set_close_to!(sl_xz, round(Int, length(new_y) / 2))
        set_close_to!(sl_xy, round(Int, length(new_z) / 2))

        current_heatmaps[].yz.visible = tog_yz.active[]
        current_heatmaps[].xz.visible = tog_xz.active[]
        current_heatmaps[].xy.visible = tog_xy.active[]

        cb.colormap[] = new_cmap

        mid_x = round(Int, length(new_x) / 2)
        mid_y = round(Int, length(new_y) / 2)
        mid_z = round(Int, length(new_z) / 2)
        lbl_yz.text[] = "$(round(new_x[mid_x], digits=0)) m"
        lbl_xz.text[] = "$(round(new_y[mid_y], digits=0)) m"
        lbl_xy.text[] = "Depth: $(round(-new_z[mid_z], digits=0)) m"

        return new_plt
    end

    on(btn_toggle.clicks) do _
        show_full_model[] = !show_full_model[]
        btn_toggle.label[] = show_full_model[] ? "Show Core" : "Show Full"

        if show_full_model[]
            ix_new = ix_full
            iy_new = iy_full
        else
            ix_new = ix_core
            iy_new = iy_core
        end

        new_x = x_all[ix_new]
        new_y = y_all[iy_new]
        new_z = -z_all[kz]
        new_R = A_all[ix_new, iy_new, kz]

        current_x[] = new_x
        current_y[] = new_y
        current_R[] = new_R

        rebuild_plot!(new_x, new_y, new_z, new_R, current_colormap[])
    end

    on(btn_reset.clicks) do _

        update_cam!(ax.scene, Vec3f(1, 1, 0.5), Vec3f(0, 0, 0), Vec3f(0, 0, 1))
    end

    function export_figure()

        ix_pos = sl_yz.value[]
        iy_pos = sl_xz.value[]
        iz_pos = sl_xy.value[]

        cur_x = current_x[]
        cur_y = current_y[]
        cur_z = current_z[]

        x_val = cur_x[ix_pos]
        y_val = cur_y[iy_pos]
        z_val = cur_z[iz_pos]

        export_fig = Figure(size = (1200, 900), fontsize = 14)

        depth_str = if abs(z_val) < 1000
            "$(round(z_val, digits=1)) m"
        else
            "$(round(z_val/1000, digits=2)) km"
        end
        Label(export_fig[1, 1:2], "3D View | XY at $depth_str, X=$(round(x_val, digits=0))m, Y=$(round(y_val, digits=0))m", 
              fontsize = 18, font = :bold)

        export_ax = LScene(export_fig[2, 1], show_axis = true)

        export_plt = volumeslices!(export_ax, cur_x, cur_y, cur_z, current_R[]; 
                                   colormap = current_colormap[], 
                                   colorrange = (cmin, cmax), 
                                   bbox_visible = true)

        export_plt[:update_yz][](ix_pos)
        export_plt[:update_xz][](iy_pos)
        export_plt[:update_xy][](iz_pos)

        export_plt[:heatmap_yz][].visible = tog_yz.active[]
        export_plt[:heatmap_xz][].visible = tog_xz.active[]
        export_plt[:heatmap_xy][].visible = tog_xy.active[]

        cb_lbl = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
        Colorbar(export_fig[2, 2], export_plt[:heatmap_xy][], label = cb_lbl, labelsize = 14)

        colsize!(export_fig.layout, 2, Relative(0.05))

        timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
        view_type = show_full_model[] ? "full" : "core"
        filename = "modem_3d_$(view_type)_$(timestamp).png"

        save(filename, export_fig, px_per_unit = 3)

        println("Figure exported: $filename")
        println("  Resolution: 3600 x 2700 pixels")
        println("  Slices: X=$(round(x_val, digits=1))m, Y=$(round(y_val, digits=1))m, Depth=$depth_str")

        return filename
    end

    on(btn_export.clicks) do _
        export_figure()
    end

    return fig, (
        ax = ax,
        sliders = (yz = sl_yz, xz = sl_xz, xy = sl_xy),
        toggles = (yz = tog_yz, xz = tog_xz, xy = tog_xy),
        colormap = current_colormap,
        data = (x = current_x, y = current_y, z = current_z, R = current_R),
        colorrange = (cmin = cmin, cmax = cmax)
    )
end

function main()
    println("Loading ModEM model from: $model_file")
    M = load_model_modem(model_file)

    println("Model dimensions: $(size(M.A))")
    println("  X cells: $(length(M.cx))")
    println("  Y cells: $(length(M.cy))")
    println("  Z cells: $(length(M.cz))")

    println("\nLaunching 3D viewer...")
    println("  Log10 scale: $log10_scale")
    println("  Colormap: $colormap")
    println("  Resistivity range: $resistivity_range")
    println("  Show padding: $show_padding")
    if !isnothing(max_depth)
        println("  Max depth: $max_depth m")
    end

    fig, parts = modem_3d_viewer(M;
        log10scale = log10_scale,
        cmap = colormap,
        withPadding = show_padding,
        max_depth = max_depth,
        pad_tol = pad_tolerance,
        resistivity_range = resistivity_range
    )

    println("\n3D Viewer Controls:")
    println("  - Drag: Rotate view")
    println("  - Scroll: Zoom in/out")
    println("  - Right-drag: Pan")
    println("  - Sliders: Move slice planes")
    println("  - Toggles: Show/hide planes")
    println("  - 'Show Core/Full': Toggle padding cells")
    println("  - 'Reset View': Reset camera angle")
    println("  - 'Export Figure': Save high-quality PNG")
    println("  - Colormap dropdown: Change color scheme")

    screen = display(fig)
    println("\nViewer is open. Close the window to exit.")
    wait(screen)

    return fig, parts
end

fig, parts = main()

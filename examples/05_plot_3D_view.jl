# Example: interactive 3D volume-slice model viewer.
# Author: @pankajkmishra
# This script visualizes XY/XZ/YZ slices with controls for depth, padding, and export.
# Use it to explore overall 3D model geometry and resistivity trends.

using Pkg
#Pkg.activate(dirname(@__DIR__))

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

default_view_direction = Vec3f(-0.90f0, 0.95f0, 0.75f0)
default_view_scale = 1.15f0

function modem_3d_viewer(
    M;
    log10scale::Bool = true,
    cmap = Reverse(:turbo),
    figsize = (1800, 920),
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

    ax = LScene(fig[1, 1], show_axis = false)
    controls = fig[2, 1] = GridLayout()
    axis_bar = fig[3, 1] = GridLayout()

    cb_label = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"

    rowsize!(fig.layout, 1, Relative(0.82))
    rowsize!(fig.layout, 2, Auto(56))
    rowsize!(fig.layout, 3, Auto(28))
    colsize!(fig.layout, 1, Relative(0.96))

    plt = volumeslices!(ax, current_x[], current_y[], current_z[], current_R[]; 
                        colormap = current_colormap[], 
                        colorrange = (cmin, cmax), 
                        bbox_visible = true)
    current_plt[] = plt
    current_heatmaps[] = (yz=plt[:heatmap_yz][], xz=plt[:heatmap_xz][], xy=plt[:heatmap_xy][])

    cam3d!(ax.scene, projectiontype = Makie.Perspective)

    cb = Colorbar(fig[1, 2], current_heatmaps[].xy, label = cb_label, width = 16)
    colsize!(fig.layout, 2, Relative(0.04))

    function fit_camera!(xv, yv, zv)
        xmin, xmax = extrema(xv)
        ymin, ymax = extrema(yv)
        zmin, zmax = extrema(zv)
        center = Vec3f((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
        spanx = Float32(xmax - xmin)
        spany = Float32(ymax - ymin)
        spanz = Float32(zmax - zmin)
        horiz = max(spanx, spany)
        zlift = max(spanz, 0.45f0 * horiz)
        dir_norm = sqrt(default_view_direction[1]^2 + default_view_direction[2]^2 + default_view_direction[3]^2)
        dirx = default_view_direction[1] / dir_norm
        diry = default_view_direction[2] / dir_norm
        dirz = default_view_direction[3] / dir_norm
        eye = center + Vec3f(default_view_scale * dirx * horiz,
                             default_view_scale * diry * horiz,
                             default_view_scale * dirz * zlift)
        update_cam!(ax.scene, eye, center, Vec3f(0, 0, 1))
    end

    fit_camera!(current_x[], current_y[], current_z[])

    Label(controls[1, 1], "X:", halign = :right, fontsize = 12)
    btn_prev_yz = Button(controls[1, 2], label = "Prev", fontsize = 10)
    sl_yz = Slider(controls[1, 3], range = 1:length(x), startvalue = round(Int, length(x) / 2), width = 180)
    btn_next_yz = Button(controls[1, 4], label = "Next", fontsize = 10)
    tog_yz = Toggle(controls[1, 5], active = true)
    lbl_yz = Label(controls[1, 6], "$(round(x[sl_yz.value[]], digits=0)) m", fontsize = 10, color = :gray35)

    Label(controls[1, 7], "Y:", halign = :right, fontsize = 12)
    btn_prev_xz = Button(controls[1, 8], label = "Prev", fontsize = 10)
    sl_xz = Slider(controls[1, 9], range = 1:length(y), startvalue = round(Int, length(y) / 2), width = 180)
    btn_next_xz = Button(controls[1, 10], label = "Next", fontsize = 10)
    tog_xz = Toggle(controls[1, 11], active = true)
    lbl_xz = Label(controls[1, 12], "$(round(y[sl_xz.value[]], digits=0)) m", fontsize = 10, color = :gray35)

    Label(controls[1, 13], "Z:", halign = :right, fontsize = 12)
    btn_prev_xy = Button(controls[1, 14], label = "Prev", fontsize = 10)
    sl_xy = Slider(controls[1, 15], range = 1:length(z), startvalue = round(Int, length(z) / 2), width = 180)
    btn_next_xy = Button(controls[1, 16], label = "Next", fontsize = 10)
    tog_xy = Toggle(controls[1, 17], active = true)
    lbl_xy = Label(controls[1, 18], "Depth: $(round(-z[sl_xy.value[]], digits=0)) m", fontsize = 10, color = :gray35)

    btn_toggle = Button(controls[1, 19], label = show_full_model[] ? "Show Core" : "Show Full", fontsize = 10)
    btn_reset = Button(controls[1, 20], label = "Reset View", fontsize = 10)
    btn_export = Button(controls[1, 21], label = "Export", fontsize = 10)

    axis_info = Observable("")
    Label(axis_bar[1, 1], axis_info, fontsize = 11, color = :gray30, halign = :left, tellwidth = false)

    rowgap!(controls, 2)
    colgap!(controls, 6)

    function update_axis_info!()
        xv = current_x[]
        yv = current_y[]
        zv = current_z[]
        sx = xv[sl_yz.value[]]
        sy = yv[sl_xz.value[]]
        sz = -zv[sl_xy.value[]]
        axis_info[] = "X: [$(round(minimum(xv), digits=0)), $(round(maximum(xv), digits=0))] m   Y: [$(round(minimum(yv), digits=0)), $(round(maximum(yv), digits=0))] m   Depth: [0, $(round(maximum(-zv), digits=0))] m   |   Slice: X=$(round(sx, digits=0)) m, Y=$(round(sy, digits=0)) m, Depth=$(round(sz, digits=0)) m"
    end

    function apply_plane_visibility!()
        yz_on = tog_yz.active[]
        xz_on = tog_xz.active[]
        xy_on = tog_xy.active[]
        if !yz_on && !xz_on && !xy_on
            tog_xy.active[] = true
            yz_on = false
            xz_on = false
            xy_on = true
        end
        if current_heatmaps[].yz !== nothing
            current_heatmaps[].yz.visible[] = yz_on
        end
        if current_heatmaps[].xz !== nothing
            current_heatmaps[].xz.visible[] = xz_on
        end
        if current_heatmaps[].xy !== nothing
            current_heatmaps[].xy.visible[] = xy_on
        end
    end

    on(sl_yz.value) do v
        if current_plt[] !== nothing
            current_plt[][:update_yz][](v)
        end
        lbl_yz.text[] = "$(round(current_x[][v], digits=0)) m"
        update_axis_info!()
    end
    on(sl_xz.value) do v
        if current_plt[] !== nothing
            current_plt[][:update_xz][](v)
        end
        lbl_xz.text[] = "$(round(current_y[][v], digits=0)) m"
        update_axis_info!()
    end
    on(sl_xy.value) do v
        if current_plt[] !== nothing
            current_plt[][:update_xy][](v)
        end
        lbl_xy.text[] = "Depth: $(round(-current_z[][v], digits=0)) m"
        update_axis_info!()
    end

    on(btn_prev_yz.clicks) do _
        set_close_to!(sl_yz, max(1, sl_yz.value[] - 1))
    end
    on(btn_next_yz.clicks) do _
        set_close_to!(sl_yz, min(length(current_x[]), sl_yz.value[] + 1))
    end
    on(btn_prev_xz.clicks) do _
        set_close_to!(sl_xz, max(1, sl_xz.value[] - 1))
    end
    on(btn_next_xz.clicks) do _
        set_close_to!(sl_xz, min(length(current_y[]), sl_xz.value[] + 1))
    end
    on(btn_prev_xy.clicks) do _
        set_close_to!(sl_xy, max(1, sl_xy.value[] - 1))
    end
    on(btn_next_xy.clicks) do _
        set_close_to!(sl_xy, min(length(current_z[]), sl_xy.value[] + 1))
    end

    on(tog_yz.active) do a
        apply_plane_visibility!()
    end
    on(tog_xz.active) do a
        apply_plane_visibility!()
    end
    on(tog_xy.active) do a
        apply_plane_visibility!()
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

        apply_plane_visibility!()

        cb.colormap[] = new_cmap

        mid_x = round(Int, length(new_x) / 2)
        mid_y = round(Int, length(new_y) / 2)
        mid_z = round(Int, length(new_z) / 2)
        lbl_yz.text[] = "$(round(new_x[mid_x], digits=0)) m"
        lbl_xz.text[] = "$(round(new_y[mid_y], digits=0)) m"
        lbl_xy.text[] = "Depth: $(round(-new_z[mid_z], digits=0)) m"
        update_axis_info!()

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
        current_z[] = new_z
        current_R[] = new_R

        rebuild_plot!(new_x, new_y, new_z, new_R, current_colormap[])
        fit_camera!(current_x[], current_y[], current_z[])
    end

    on(btn_reset.clicks) do _
        fit_camera!(current_x[], current_y[], current_z[])
    end

    update_axis_info!()

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

        xedges_export = edges_from_centers(cur_x)
        yedges_export = edges_from_centers(cur_y)
        zedges_export = edges_from_centers(cur_z)
        xlims_export = extrema(xedges_export)
        ylims_export = extrema(yedges_export)
        zlims_export = extrema(zedges_export)

        export_ax = Axis3(export_fig[2, 1];
            xlabel = "x",
            ylabel = "y",
            zlabel = "z",
            limits = (xlims_export, ylims_export, zlims_export),
            xautolimitmargin = (0.0f0, 0.0f0),
            yautolimitmargin = (0.0f0, 0.0f0),
            zautolimitmargin = (0.0f0, 0.0f0),
            aspect = :data,
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xticksize = 5,
            yticksize = 5,
            zticksize = 5,
            xticklabelsize = 13,
            yticklabelsize = 13,
            zticklabelsize = 13,
            xticklabelfont = :bold,
            yticklabelfont = :bold,
            zticklabelfont = :bold,
            xlabelsize = 14,
            ylabelsize = 14,
            zlabelsize = 14,
            xlabelfont = :bold,
            ylabelfont = :bold,
            zlabelfont = :bold,
            xticklabelcolor = :gray70,
            yticklabelcolor = :gray70,
            zticklabelcolor = :gray70,
            xtickcolor = :gray70,
            ytickcolor = :gray70,
            ztickcolor = :gray70,
            xlabelcolor = :gray65,
            ylabelcolor = :gray65,
            zlabelcolor = :gray65
        )

        export_plt = volumeslices!(export_ax, cur_x, cur_y, cur_z, current_R[]; 
                                   colormap = current_colormap[], 
                                   colorrange = (cmin, cmax), 
                       bbox_visible = false)

        export_plt[:update_yz][](ix_pos)
        export_plt[:update_xz][](iy_pos)
        export_plt[:update_xy][](iz_pos)

        export_plt[:heatmap_yz][].visible[] = tog_yz.active[]
        export_plt[:heatmap_xz][].visible[] = tog_xz.active[]
        export_plt[:heatmap_xy][].visible[] = tog_xy.active[]

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
    println("  - Prev/Next + sliders: Step or move slice planes")
    println("  - Toggles: Show/hide planes")
    println("  - 'Show Core/Full': Toggle padding cells")
    println("  - 'Reset View': Reset camera angle")
    println("  - 'Export Figure': Save high-quality PNG")

    screen = display(fig)
    println("\nViewer is open. Close the window to exit.")
    wait(screen)

    return fig, parts
end

fig, parts = main()

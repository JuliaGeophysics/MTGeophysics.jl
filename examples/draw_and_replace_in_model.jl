# Example: polygon-based manual model editing.
# Author: @pankajkmishra
# This script lets you draw zones on slices and replace resistivity within selected depth intervals.
# Use it for targeted interactive edits with visual feedback.

using Pkg

using GLMakie
using Statistics
using Dates

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "CoreUtils3D.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))

model_file = joinpath(@__DIR__, "Cascadia", "cascad_half_inverse.ws")

replacement_resistivity = 10000.0
layers_above = 2
layers_below = 2
transition_layers = 1
apply_to_all_depths = false
depth_range = (0.0, 50000.0)

log10_scale = true
colormap = Reverse(:turbo)
with_padding = true
max_depth = nothing
resistivity_range = (0.0, 4.0)
show_grid = true
grid_color = :black
grid_linewidth = 0.5 
grid_alpha = 0.3
pad_tol = 0.5

function point_in_polygon(px, py, poly_x, poly_y)
    n = length(poly_x)
    inside = false
    j = n
    for i in 1:n
        xi, yi = poly_x[i], poly_y[i]
        xj, yj = poly_x[j], poly_y[j]
        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function get_cells_in_polygon(poly_x, poly_y, cx, cy)
    mask = falses(length(cx), length(cy))
    for (i, x) in enumerate(cx)
        for (j, y) in enumerate(cy)
            if point_in_polygon(x, y, poly_x, poly_y)
                mask[i, j] = true
            end
        end
    end
    return mask
end

function apply_zone_modification!(A, mask, center_layer, layers_above, layers_below, 
                                   transition_layers, replacement_value, nz)
    count = 0

    core_start = max(1, center_layer - layers_above)
    core_end = min(nz, center_layer + layers_below)

    for k in core_start:core_end
        for j in 1:size(A, 2)
            for i in 1:size(A, 1)
                if mask[i, j]
                    A[i, j, k] = replacement_value
                    count += 1
                end
            end
        end
    end

    if transition_layers > 0
        trans_above_start = max(1, core_start - transition_layers)
        for k in trans_above_start:(core_start-1)
            weight = (k - trans_above_start + 1) / (transition_layers + 1)
            for j in 1:size(A, 2)
                for i in 1:size(A, 1)
                    if mask[i, j]
                        A[i, j, k] = (1 - weight) * A[i, j, k] + weight * replacement_value
                        count += 1
                    end
                end
            end
        end

        trans_below_end = min(nz, core_end + transition_layers)
        for k in (core_end+1):trans_below_end
            weight = 1 - (k - core_end) / (transition_layers + 1)
            for j in 1:size(A, 2)
                for i in 1:size(A, 1)
                    if mask[i, j]
                        A[i, j, k] = (1 - weight) * A[i, j, k] + weight * replacement_value
                        count += 1
                    end
                end
            end
        end
    end

    return count, core_start, core_end
end

function apply_zone_all_depths!(A, mask, replacement_value)
    count = 0
    for k in 1:size(A, 3)
        for j in 1:size(A, 2)
            for i in 1:size(A, 1)
                if mask[i, j]
                    A[i, j, k] = replacement_value
                    count += 1
                end
            end
        end
    end
    return count
end

function get_depth_layer_indices(layer_depths, depth_range)
    d_min, d_max = depth_range
    indices = Int[]
    for (k, depth) in enumerate(layer_depths)
        prev_depth = k == 1 ? 0.0 : layer_depths[k-1]
        if prev_depth < d_max && depth > d_min
            push!(indices, k)
        end
    end
    return indices
end

function zone_editor(M_input; model_name="model", M=M_input, replacement_resistivity=1000.0,
                     apply_to_all_depths=false, depth_range=(0.0, 50000.0),
                     log10scale=true, cmap=Reverse(:turbo), figsize=(1200, 950),
                     withPadding=true, max_depth=nothing, pad_tol=0.2,
                     resistivity_range=nothing, show_grid=true, grid_color=:black,
                     grid_linewidth=0.5, grid_alpha=0.3)

    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_working = copy(M.A)
    A_display = log10scale ? log10.(A_working) : copy(A_working)

    ix_full = 1:length(x_all)
    iy_full = 1:length(y_all)
    ix_core = core_indices(x_all; tol=pad_tol)
    iy_core = core_indices(y_all; tol=pad_tol)

    kz = isnothing(max_depth) ? (1:length(z_all)) : z_indices_for_max_depth(z_all, float(max_depth))
    z = z_all[kz]

    x_full = x_all[ix_full]
    y_full = y_all[iy_full]

    x_edges_full = edges_from_centers(x_full)
    y_edges_full = edges_from_centers(y_full)
    z_edges = edges_from_centers(z)
    layer_depths = cumsum(diff(z_edges))

    if isnothing(resistivity_range)
        vals = A_display[isfinite.(A_display)]
        qlo, qhi = quantile(vec(vals), (0.02, 0.98))
        cmin, cmax = qlo, qhi
    else
        cmin, cmax = resistivity_range
    end

    current_x_edges = Observable(x_edges_full)
    current_y_edges = Observable(y_edges_full)

    fig = Figure(size=figsize)

    current_layer = Observable(1)
    slice_data = Observable(A_display[:, :, 1]')

    title_str = Observable("Zone Editor - Layer 1 / $(length(z))")
    Label(fig[0, 1:2], title_str, fontsize=18, font=:bold)

    ax = Axis(fig[1, 1], xlabel="Y (m)", ylabel="X (m)", aspect=DataAspect())

    current_colormap = Observable(cmap)
    hm = heatmap!(ax, current_y_edges[], current_x_edges[], slice_data,
                  colormap=current_colormap, colorrange=(cmin, cmax))

    cb_label = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
    Colorbar(fig[1, 2], hm, label=cb_label)

    polygon_points = Observable(Point2f[])
    polygon_closed = Observable(false)

    polygon_lines = @lift length($polygon_points) >= 2 ? $polygon_points : Point2f[]

    scatter!(ax, polygon_points, color=:red, markersize=10)
    lines!(ax, polygon_lines, color=:red, linewidth=2)

    closed_line = Observable(Point2f[])
    lines!(ax, closed_line, color=:red, linewidth=2, linestyle=:dash)

    slider_grid = fig[2, 1:2] = GridLayout()
    btn_prev = Button(slider_grid[1, 1], label="<< Prev")
    Label(slider_grid[1, 2], "Depth Layer:", fontsize=14)
    sl = Slider(slider_grid[1, 3], range=1:length(z), startvalue=1, width=400)
    layer_label = Observable("1 / $(length(z))")
    Label(slider_grid[1, 4], layer_label, fontsize=14)
    btn_next = Button(slider_grid[1, 5], label="Next >>")

    button_grid = fig[3, 1:2] = GridLayout()
    btn_clear = Button(button_grid[1, 1], label="Clear Polygon")
    btn_apply = Button(button_grid[1, 2], label="Apply to Zone")
    btn_undo = Button(button_grid[1, 3], label="Undo Last")
    btn_reset = Button(button_grid[1, 4], label="Reset Model")
    btn_export = Button(button_grid[1, 5], label="Export Figure")

    options_grid = fig[4, 1:2] = GridLayout()
    Label(options_grid[1, 1], "ρ (Ω·m):", fontsize=12)
    resistivity_tb = Textbox(options_grid[1, 2], placeholder="10000", 
                             stored_string=string(Int(replacement_resistivity)), width=80)

    Label(options_grid[1, 3], "Above:", fontsize=12)
    above_tb = Textbox(options_grid[1, 4], stored_string=string(layers_above), width=40)

    Label(options_grid[1, 5], "Below:", fontsize=12)
    below_tb = Textbox(options_grid[1, 6], stored_string=string(layers_below), width=40)

    Label(options_grid[1, 7], "Trans:", fontsize=12)
    trans_tb = Textbox(options_grid[1, 8], stored_string=string(transition_layers), width=40)

    Label(options_grid[1, 9], "All depths:", fontsize=12)
    depth_toggle = Toggle(options_grid[1, 10], active=apply_to_all_depths)

    info_grid = fig[5, 1:2] = GridLayout()
    status_text = Observable("Click to draw polygon vertices. Right-click to close polygon.")
    Label(info_grid[1, 1], status_text, fontsize=12, color=:blue)

    history_stack = Vector{Array{Float64,3}}()
    push!(history_stack, copy(A_working))

    function update_display()
        A_display = log10scale ? log10.(A_working) : copy(A_working)
        layer_idx = current_layer[]
        slice_data[] = A_display[:, :, layer_idx]'
    end

    function update_slice(layer_idx)
        layer_idx = clamp(layer_idx, 1, length(z))
        current_layer[] = layer_idx
        update_display()
        layer_label[] = "$layer_idx / $(length(z))"

        depth_top = layer_idx == 1 ? 0.0 : layer_depths[layer_idx-1]
        depth_bottom = layer_depths[layer_idx]
        title_str[] = "Zone Editor - Layer $layer_idx / $(length(z)) | Depth: $(round(depth_top, digits=1)) - $(round(depth_bottom, digits=1)) m"
    end

    on(sl.value) do val
        update_slice(val)
    end

    on(btn_prev.clicks) do _
        set_close_to!(sl, max(1, sl.value[] - 1))
    end

    on(btn_next.clicks) do _
        set_close_to!(sl, min(length(z), sl.value[] + 1))
    end

    y_min, y_max = extrema(y_all)
    x_min, x_max = extrema(x_all)

    deregister_interaction!(ax, :rectanglezoom)

    on(events(ax).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            if !polygon_closed[]
                pos = mouseposition(ax)
                y_click, x_click = pos[1], pos[2]

                if y_click < y_min || y_click > y_max || x_click < x_min || x_click > x_max
                    status_text[] = "Click outside model bounds - ignored."
                    return
                end

                pts = copy(polygon_points[])
                push!(pts, Point2f(y_click, x_click))
                polygon_points[] = pts
                status_text[] = "$(length(pts)) vertices. Right-click to close, or continue adding."
            end
        elseif event.button == Mouse.right && event.action == Mouse.press
            pts = polygon_points[]
            if length(pts) >= 3
                polygon_closed[] = true
                closed_line[] = [pts[end], pts[1]]
                status_text[] = "Polygon closed with $(length(pts)) vertices. Click 'Apply to Zone' to modify."
            else
                status_text[] = "Need at least 3 vertices to close polygon."
            end
        end
    end

    on(btn_clear.clicks) do _
        polygon_points[] = Point2f[]
        closed_line[] = Point2f[]
        polygon_closed[] = false
        status_text[] = "Polygon cleared. Click to draw new vertices."
    end

    on(btn_apply.clicks) do _
        pts = polygon_points[]
        if length(pts) < 3
            status_text[] = "Need at least 3 vertices to define a zone."
            return
        end

        push!(history_stack, copy(A_working))

        poly_y = [p[1] for p in pts]
        poly_x = [p[2] for p in pts]

        println("Polygon vertices (Y, X):")
        for (i, p) in enumerate(pts)
            println("  $i: Y=$(p[1]), X=$(p[2])")
        end
        println("X range: $(minimum(x_all)) to $(maximum(x_all))")
        println("Y range: $(minimum(y_all)) to $(maximum(y_all))")

        mask = get_cells_in_polygon(poly_x, poly_y, x_all, y_all)
        println("Mask cells selected: $(sum(mask))")

        repl_val = try
            parse(Float64, resistivity_tb.stored_string[])
        catch
            replacement_resistivity
        end

        n_above = try
            parse(Int, above_tb.stored_string[])
        catch
            layers_above
        end

        n_below = try
            parse(Int, below_tb.stored_string[])
        catch
            layers_below
        end

        n_trans = try
            parse(Int, trans_tb.stored_string[])
        catch
            transition_layers
        end

        n_cells_2d = sum(mask)

        if depth_toggle.active[]
            count = apply_zone_all_depths!(A_working, mask, repl_val)
            status_text[] = "Applied to ALL depths: $(n_cells_2d) cells × $(length(z)) layers = $(count) cells. ρ = $(repl_val) Ω·m"
        else
            count, core_start, core_end = apply_zone_modification!(
                A_working, mask, current_layer[], n_above, n_below, n_trans, repl_val, length(z))
            status_text[] = "Applied: $(n_cells_2d) cells, layers $(core_start)-$(core_end) + $(n_trans) trans. ρ = $(repl_val) Ω·m"
        end

        update_display()

        polygon_points[] = Point2f[]
        closed_line[] = Point2f[]
        polygon_closed[] = false
    end

    on(btn_undo.clicks) do _
        if length(history_stack) > 1
            pop!(history_stack)
            A_working .= history_stack[end]
            update_display()
            status_text[] = "Undo successful. $(length(history_stack)-1) changes remaining."
        else
            status_text[] = "Nothing to undo."
        end
    end

    on(btn_reset.clicks) do _
        A_working .= history_stack[1]
        history_stack = [copy(A_working)]
        update_display()
        polygon_points[] = Point2f[]
        closed_line[] = Point2f[]
        polygon_closed[] = false
        status_text[] = "Model reset to original."
    end

    on(btn_export.clicks) do _
        layer_idx = current_layer[]
        filename = "$(model_name)_zone_layer$(layer_idx).png"

        export_fig = Figure(size=(900, 900), fontsize=16)
        depth_top = layer_idx == 1 ? 0.0 : layer_depths[layer_idx-1]
        depth_bottom = layer_depths[layer_idx]
        Label(export_fig[1, 1:2], "Layer $layer_idx | Depth: $(round(depth_top, digits=1)) - $(round(depth_bottom, digits=1)) m", 
              fontsize=20, font=:bold)

        export_ax = Axis(export_fig[2, 1], xlabel="Y (m)", ylabel="X (m)", aspect=AxisAspect(1))

        A_exp = log10scale ? log10.(A_working) : copy(A_working)
        heatmap!(export_ax, y_edges_full, x_edges_full, A_exp[:, :, layer_idx]',
                 colormap=current_colormap[], colorrange=(cmin, cmax))

        cb_lbl = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
        Colorbar(export_fig[2, 2], colormap=current_colormap[], colorrange=(cmin, cmax), 
                 label=cb_lbl, labelsize=16)
        colsize!(export_fig.layout, 2, Relative(0.05))

        save(filename, export_fig, px_per_unit=3)
        status_text[] = "Exported: $filename"
    end

    save_grid = fig[6, 1:2] = GridLayout()
    btn_save = Button(save_grid[1, 1], label = "Save Model", fontsize = 14)
    save_label = Observable("")
    Label(save_grid[1, 2], save_label, fontsize = 12, color = :green)

    on(btn_save.clicks) do _
        repl_val = try
            parse(Float64, resistivity_tb.stored_string[])
        catch
            replacement_resistivity
        end

        n_above = try
            parse(Int, above_tb.stored_string[])
        catch
            layers_above
        end

        n_below = try
            parse(Int, below_tb.stored_string[])
        catch
            layers_below
        end

        n_trans = try
            parse(Int, trans_tb.stored_string[])
        catch
            transition_layers
        end

        mode_tag = depth_toggle.active[] ? "allz" : "layer$(current_layer[])"
        rho_tag = "rho$(Int(round(repl_val)))"
        above_tag = "ab$(n_above)"
        below_tag = "bl$(n_below)"
        trans_tag = "tr$(n_trans)"

        outname = splitext(M.name)[1] * "_edited_$(rho_tag)_$(above_tag)_$(below_tag)_$(trans_tag)_$(mode_tag).rho"
        write_model_modem(outname, M.dx, M.dy, M.dz, A_working, M.origin)
        save_label[] = "Saved: $(basename(outname))"
        status_text[] = "Model saved to: $outname"
    end

    return fig, (ax=ax, A_working=A_working, history=history_stack, 
                 current_layer=current_layer, layer_depths=layer_depths)
end

function main()
    if !isfile(model_file)
        println("ERROR: Model file not found: $model_file")
        return nothing, nothing
    end

    println("Loading ModEM model: $model_file")
    M = load_model_modem(model_file)

    println("Model loaded:")
    println("  Grid size: $(M.nx) × $(M.ny) × $(M.nz)")

    model_name = splitext(basename(model_file))[1]

    println("\nCreating zone editor...")

    fig, parts = zone_editor(M;
        model_name=model_name,
        M=M,
        replacement_resistivity=replacement_resistivity,
        apply_to_all_depths=apply_to_all_depths,
        depth_range=depth_range,
        log10scale=log10_scale,
        cmap=colormap,
        withPadding=with_padding,
        max_depth=max_depth,
        pad_tol=pad_tol,
        resistivity_range=resistivity_range,
        show_grid=show_grid,
        grid_color=grid_color,
        grid_linewidth=grid_linewidth,
        grid_alpha=grid_alpha
    )

    println("\nZone Editor ready!")
    println("  - Left-click to add polygon vertices")
    println("  - Right-click to close polygon")
    println("  - Click 'Apply to Zone' to replace resistivity")
    println("  - Toggle to apply to current layer or all depths")

    screen = display(fig)
    println("\nClose the figure window to exit...")
    wait(screen)

    return fig, parts
end

fig, parts = main()

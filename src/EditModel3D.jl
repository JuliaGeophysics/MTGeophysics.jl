function _edit_depth_slice_viewer(
    M;
    model_name::String = "model",
    log10scale::Bool = true,
    cmap = :turbo,
    figsize = (1100, 950),
    withPadding::Bool = true,
    max_depth::Union{Nothing, Real} = nothing,
    pad_tol::Real = 0.2,
    resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}} = nothing,
    show_grid::Bool = true,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.3,
    axis_xlabel::String = "Y (m)",
    axis_ylabel::String = "X (m)",
    axis_aspect::Union{Nothing, Real} = nothing
)
    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : copy(M.A)

    ix_full = 1:length(x_all)
    iy_full = 1:length(y_all)
    ix_core, iy_core = lateral_core_ranges(M; tol = pad_tol)

    if isnothing(max_depth)
        kz = 1:length(z_all)
    else
        kz = z_indices_for_max_depth(z_all, float(max_depth))
    end
    z = z_all[kz]

    x_full = x_all[ix_full]
    y_full = y_all[iy_full]
    R_full = A_all[ix_full, iy_full, kz]

    x_core = x_all[ix_core]
    y_core = y_all[iy_core]
    R_core = A_all[ix_core, iy_core, kz]

    x_edges_full = edges_from_centers(x_full)
    y_edges_full = edges_from_centers(y_full)
    x_edges_core = edges_from_centers(x_core)
    y_edges_core = edges_from_centers(y_core)
    z_edges = edges_from_centers(z)

    layer_depths = cumsum(diff(z_edges))

    if isnothing(resistivity_range)
        vals = R_full[isfinite.(R_full)]
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
        cmin, cmax = resistivity_range[1], resistivity_range[2]
        if cmin > cmax
            cmin, cmax = cmax, cmin
        end
        if cmin == cmax
            ϵ = max(1e-12, 1e-6 * abs(cmin))
            cmin -= ϵ; cmax += ϵ
        end
    end

    show_full_model = Observable(withPadding)

    current_x_edges = Observable(withPadding ? x_edges_full : x_edges_core)
    current_y_edges = Observable(withPadding ? y_edges_full : y_edges_core)
    current_R = withPadding ? R_full : R_core

    fig = Figure(size = figsize)

    title_str = Observable("Depth Layer 1 / $(length(z)) | Depth: 0 - $(round(layer_depths[1], digits=1)) m")
    Label(fig[0, 1:2], title_str, fontsize = 18, font = :bold)

    ax_aspect = isnothing(axis_aspect) ? DataAspect() : AxisAspect(Float64(axis_aspect))

    ax = Axis(fig[1, 1],
        xlabel = axis_xlabel,
        ylabel = axis_ylabel,
        aspect = ax_aspect,
        title = ""
    )

    current_layer = Observable(1)
    slice_data = Observable(current_R[:, :, 1]')

    current_colormap = Observable(cmap)

    hm = heatmap!(ax, current_y_edges[], current_x_edges[], slice_data,
                  colormap = current_colormap,
                  colorrange = (cmin, cmax))

    grid_plots = Ref{Vector{Any}}([])

    function draw_grid!(ax, x_edges, y_edges)
        for p in grid_plots[]
            delete!(ax, p)
        end
        grid_plots[] = []

        if show_grid
            for xe in x_edges
                p = lines!(ax, [y_edges[1], y_edges[end]], [xe, xe],
                    color = (grid_color, grid_alpha),
                    linewidth = grid_linewidth)
                push!(grid_plots[], p)
            end
            for ye in y_edges
                p = lines!(ax, [ye, ye], [x_edges[1], x_edges[end]],
                    color = (grid_color, grid_alpha),
                    linewidth = grid_linewidth)
                push!(grid_plots[], p)
            end
        end
    end

    draw_grid!(ax, current_x_edges[], current_y_edges[])

    cb_label = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
    Colorbar(fig[1, 2], hm, label = cb_label)

    slider_grid = fig[2, 1:2] = GridLayout()

    btn_prev = Button(slider_grid[1, 1], label = "<< Prev")
    Label(slider_grid[1, 2], "Depth Layer:", fontsize = 14)
    sl = Slider(slider_grid[1, 3], range = 1:length(z), startvalue = 1, width = 400)
    layer_label = Observable("1 / $(length(z))")
    Label(slider_grid[1, 4], layer_label, fontsize = 14)
    btn_next = Button(slider_grid[1, 5], label = "Next >>")

    button_grid = fig[3, 1:2] = GridLayout()

    btn_first = Button(button_grid[1, 1], label = "|<< First")
    btn_last = Button(button_grid[1, 2], label = "Last >>|")

    btn_label = Observable(withPadding ? "Show Core Model" : "Show Full Model")
    btn_toggle = Button(button_grid[1, 3], label = btn_label)

    btn_reset = Button(button_grid[1, 4], label = "Reset Zoom")

    btn_export = Button(button_grid[1, 5], label = "Export Figure")

    info_grid = fig[4, 1:2] = GridLayout()
    depth_info = Observable("Layer depth range: 0 - $(round(layer_depths[1], digits=1)) m | Cell thickness: $(round(diff(z_edges)[1], digits=1)) m")
    Label(info_grid[1, 1], depth_info, fontsize = 12)

    view_info = Observable(withPadding ? "View: Full Model (with padding)" : "View: Core Model (no padding)")
    Label(info_grid[1, 2], view_info, fontsize = 12, color = :blue)

    current_R_ref = Ref(current_R)

    function reset_zoom!()
        xe = current_x_edges[]
        ye = current_y_edges[]
        limits!(ax, ye[1], ye[end], xe[1], xe[end])
    end

    function update_view!(show_full::Bool)
        if show_full
            current_x_edges[] = x_edges_full
            current_y_edges[] = y_edges_full
            current_R_ref[] = R_full
            btn_label[] = "Show Core Model"
            view_info[] = "View: Full Model (with padding)"
        else
            current_x_edges[] = x_edges_core
            current_y_edges[] = y_edges_core
            current_R_ref[] = R_core
            btn_label[] = "Show Full Model"
            view_info[] = "View: Core Model (no padding)"
        end

        layer_idx = current_layer[]
        slice_data[] = current_R_ref[][:, :, layer_idx]'

        empty!(ax)
        hm = heatmap!(ax, current_y_edges[], current_x_edges[], slice_data,
            colormap = current_colormap,
            colorrange = (cmin, cmax))

        draw_grid!(ax, current_x_edges[], current_y_edges[])

        reset_zoom!()
    end

    function update_slice(layer_idx)
        layer_idx = clamp(layer_idx, 1, length(z))
        current_layer[] = layer_idx
        slice_data[] = current_R_ref[][:, :, layer_idx]'

        layer_label[] = "$layer_idx / $(length(z))"

        if layer_idx == 1
            depth_top = 0.0
        else
            depth_top = layer_depths[layer_idx - 1]
        end
        depth_bottom = layer_depths[layer_idx]
        thickness = diff(z_edges)[layer_idx]

        title_str[] = "Depth Layer $layer_idx / $(length(z)) | Depth: $(round(depth_top, digits=1)) - $(round(depth_bottom, digits=1)) m"
        depth_info[] = "Layer depth range: $(round(depth_top, digits=1)) - $(round(depth_bottom, digits=1)) m | Cell thickness: $(round(thickness, digits=1)) m"
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

    on(btn_first.clicks) do _
        set_close_to!(sl, 1)
    end

    on(btn_last.clicks) do _
        set_close_to!(sl, length(z))
    end

    on(btn_toggle.clicks) do _
        show_full_model[] = !show_full_model[]
        update_view!(show_full_model[])
    end

    on(btn_reset.clicks) do _
        reset_zoom!()
    end

    function export_figure()
        layer_idx = current_layer[]

        if layer_idx == 1
            depth_top = 0.0
        else
            depth_top = layer_depths[layer_idx - 1]
        end
        depth_bottom = layer_depths[layer_idx]

        xe = current_x_edges[]
        ye = current_y_edges[]
        data = slice_data[]
        cmap_val = current_colormap[]

        export_fig = Figure(size = (900, 900), fontsize = 16)

        depth_str = if depth_bottom < 1000
            "Depth: $(round(depth_top, digits=1)) - $(round(depth_bottom, digits=1)) m"
        else
            "Depth: $(round(depth_top/1000, digits=2)) - $(round(depth_bottom/1000, digits=2)) km"
        end
        Label(export_fig[1, 1:2], "Layer $layer_idx | $depth_str",
            fontsize = 20, font = :bold)

        export_ax = Axis(export_fig[2, 1],
            xlabel = axis_xlabel,
            ylabel = axis_ylabel,
            aspect = ax_aspect,
            xlabelsize = 16,
            ylabelsize = 16,
            xticklabelsize = 12,
            yticklabelsize = 12
        )

        export_hm = heatmap!(export_ax, ye, xe, data,
            colormap = cmap_val,
            colorrange = (cmin, cmax))

        if show_grid
            for x_edge in xe
                lines!(export_ax, [ye[1], ye[end]], [x_edge, x_edge],
                    color = (grid_color, grid_alpha),
                    linewidth = grid_linewidth)
            end
            for y_edge in ye
                lines!(export_ax, [y_edge, y_edge], [xe[1], xe[end]],
                    color = (grid_color, grid_alpha),
                    linewidth = grid_linewidth)
            end
        end

        cb_lbl = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
        Colorbar(export_fig[2, 2], export_hm, label = cb_lbl, labelsize = 16, ticklabelsize = 12)

        colsize!(export_fig.layout, 2, Relative(0.05))

        view_type = show_full_model[] ? "full" : "core"
        filename = "$(model_name)_layer$(layer_idx)_$(view_type).png"

        save(filename, export_fig, px_per_unit = 3)

        println("Figure exported: $filename")
        println("  Resolution: 3600 x 3000 pixels")
        println("  Layer: $layer_idx, Depth: $depth_str")

        return filename
    end

    on(btn_export.clicks) do _
        export_figure()
    end

    return fig, (
        ax = ax,
        heatmap = hm,
        slider = sl,
        current_layer = current_layer,
        slice_data = slice_data,
        show_full_model = show_full_model,
        R_full = R_full,
        R_core = R_core,
        x_full = x_full, y_full = y_full,
        x_core = x_core, y_core = y_core,
        z = z,
        x_edges_full = x_edges_full, y_edges_full = y_edges_full,
        x_edges_core = x_edges_core, y_edges_core = y_edges_core,
        z_edges = z_edges,
        layer_depths = layer_depths,
        colorrange = (cmin = cmin, cmax = cmax),
        reset_zoom! = reset_zoom!
    )
end

function get_scope_indices(M; replace_scope::Symbol = :core_only)
    if replace_scope == :full_model
        ix = 1:M.nx
        iy = 1:M.ny
        scope_label = "full model including padding"
    elseif replace_scope == :core_only
        ix = (M.npad[1] + 1):(M.nx - M.npad[1])
        iy = (M.npad[2] + 1):(M.ny - M.npad[2])
        scope_label = "core model only"
    else
        error("replace_scope must be :core_only or :full_model, got $(replace_scope)")
    end

    return ix, iy, scope_label
end

function apply_resistivity_on_scope!(A::Array{<:Real,3}, ix, iy, k::Int, target_resistivity::Real)
    A[ix, iy, k] .= target_resistivity
    return A
end

function blend_resistivity_on_scope!(A::Array{<:Real,3}, ix, iy, k::Int, weight::Real, target_resistivity::Real)
    current_values = A[ix, iy, k]
    blended_values = 10 .^ ((1 - weight) .* log10.(current_values) .+ weight .* log10(target_resistivity))
    A[ix, iy, k] .= blended_values
    return A
end

function modify_deep_resistivity!(A::Array{<:Real,3}, M;
                                   cutoff_layer::Int = 1,
                                   target_resistivity::Real = 1000.0,
                                   blend_previous_percent::Int = 0,
                                   replace_scope::Symbol = :core_only)
    nz = size(A, 3)

    if cutoff_layer < 1 || cutoff_layer > nz
        println("Warning: cutoff_layer $(cutoff_layer) is outside model layers (1:$(nz)). No modification made.")
        return A, cutoff_layer, false, replace_scope
    end

    ix, iy, scope_label = get_scope_indices(M; replace_scope = replace_scope)

    println("Modifying resistivity:")
    println("  Cutoff layer: $(cutoff_layer) / $(nz)")
    println("  Target resistivity: $(target_resistivity) ohm-m")
    println("  Blend previous layer: $(blend_previous_percent)%")
    println("  Blend mode: log10 resistivity")
    println("  Replace scope: $(scope_label)")
    println("  XY cells affected per layer: $(length(ix)) × $(length(iy))")

    if blend_previous_percent > 0 && cutoff_layer > 1
        blend_weight = blend_previous_percent / 100
        previous_layer = cutoff_layer - 1
        blend_resistivity_on_scope!(A, ix, iy, previous_layer, blend_weight, target_resistivity)
        println("  Blended layer $(previous_layer) with weight $(round(blend_weight, digits = 3))")
    end

    for k in cutoff_layer:nz
        apply_resistivity_on_scope!(A, ix, iy, k, target_resistivity)
    end
    println("  Replaced layers $(cutoff_layer) - $(nz) with target value")

    return A, cutoff_layer, true, replace_scope
end

function refresh_viewer_data!(parts, A_modified, M;
                              log10_scale::Bool = true,
                              max_depth = nothing,
                              pad_tol::Real = 0.5)
    M_view = (
        A = A_modified,
        cx = M.cx,
        cy = M.cy,
        cz = M.cz,
        npad = M.npad
    )

    _, _, _, R_full_new, _, _, _ = prepare_model_arrays(M_view;
        log10scale = log10_scale,
        withPadding = true,
        max_depth = max_depth,
        pad_tol = pad_tol
    )

    _, _, _, R_core_new, _, _, _ = prepare_model_arrays(M_view;
        log10scale = log10_scale,
        withPadding = false,
        max_depth = max_depth,
        pad_tol = pad_tol
    )

    parts.R_full .= R_full_new
    parts.R_core .= R_core_new

    layer_idx = parts.current_layer[]
    if parts.show_full_model[]
        parts.slice_data[] = parts.R_full[:, :, layer_idx]'
    else
        parts.slice_data[] = parts.R_core[:, :, layer_idx]'
    end

    return nothing
end

function write_modified_model_with_metadata(outputfile::AbstractString, M, A_modified;
                                            cutoff_layer::Int,
                                            target_resistivity::Real,
                                            blend_previous_percent::Int,
                                            replace_scope::Symbol)
    write_model_modem(outputfile, M.dx, M.dy, M.dz, A_modified, M.origin)

    compact_header = "# 3D MT model written by ModEM in WS format | edit03a layer=$(cutoff_layer) rho=$(target_resistivity) blendprev=$(blend_previous_percent) log10 scope=$(replace_scope)"

    lines = readlines(outputfile)
    if isempty(lines)
        open(outputfile, "w") do io
            println(io, compact_header)
        end
    else
        lines[1] = compact_header
        open(outputfile, "w") do io
            for line in lines
                println(io, line)
            end
        end
    end

    return outputfile
end

function textbox_text(tb)
    if hasproperty(tb, :displayed_string)
        return strip(tb.displayed_string[])
    elseif hasproperty(tb, :stored_string)
        return strip(tb.stored_string[])
    else
        return ""
    end
end

function set_textbox_text!(tb, value::AbstractString)
    if hasproperty(tb, :stored_string)
        tb.stored_string[] = value
    end
    if hasproperty(tb, :displayed_string)
        tb.displayed_string[] = value
    end
    return nothing
end

function EditModelByLayers(model_file::AbstractString;
    target_resistivity::Real = 1000.0,
    blend_previous_percent::Int = 0,
    log10_scale::Bool = true,
    colormap = Reverse(:turbo),
    with_padding::Bool = true,
    max_depth::Union{Nothing, Real} = nothing,
    resistivity_range = (0.0, 4.0),
    show_grid::Bool = true,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.3,
    pad_tol::Real = 0.5,
    viewer_figsize = (1100, 950),
    interactive::Bool = true)

    isempty(model_file) && error("model_file is required")
    isfile(model_file) || error("Model file not found: $model_file")

    println("Loading ModEM model: $model_file")
    M = load_model_modem(model_file)

    println("Model loaded:")
    println("  Grid size: $(M.nx) × $(M.ny) × $(M.nz)")
    println("  Padding cells: $(M.npad)")

    A_original = copy(M.A)
    A_modified = copy(M.A)

    cutoff_layer_used = 1
    modification_applied = false
    replace_scope_used = with_padding ? :full_model : :core_only

    M_modified = (
        A = A_modified,
        x = M.x, y = M.y, z = M.z,
        cx = M.cx, cy = M.cy, cz = M.cz,
        nx = M.nx, ny = M.ny, nz = M.nz,
        npad = M.npad
    )

    model_name = splitext(basename(model_file))[1] * "_modified"

    println("\nCreating viewer for modified model...")

    fig, parts = _edit_depth_slice_viewer(M_modified;
        figsize = viewer_figsize,
        model_name = model_name,
        log10scale = log10_scale,
        cmap = colormap,
        withPadding = with_padding,
        max_depth = max_depth,
        pad_tol = pad_tol,
        resistivity_range = resistivity_range,
        show_grid = show_grid,
        grid_color = grid_color,
        grid_linewidth = grid_linewidth,
        grid_alpha = grid_alpha
    )

    save_grid = fig[4, 1:2] = GridLayout()
    btn_save = Button(save_grid[1, 1], label = "Save Model", fontsize = 14)
    save_label = Observable("")
    Label(save_grid[1, 2], save_label, fontsize = 12, color = :green)

    controls_grid = fig[5, 1:2] = GridLayout()

    Label(controls_grid[1, 1], "Current Layer", fontsize = 12)
    current_layer_label = Observable(string(parts.current_layer[]))
    Label(controls_grid[1, 2], current_layer_label, fontsize = 12)

    Label(controls_grid[1, 3], "Target ρ (Ω·m)", fontsize = 12)
    resistivity_tb = Textbox(controls_grid[1, 4], stored_string = string(Int(round(target_resistivity))), validator = Int, width = 110)

    Label(controls_grid[1, 5], "Blend Previous Layer (%)", fontsize = 12)
    transition_tb = Textbox(controls_grid[1, 6], stored_string = string(blend_previous_percent), validator = Int, width = 80)

    btn_apply = Button(controls_grid[1, 7], label = "Apply Changes")
    btn_reset_edits = Button(controls_grid[1, 8], label = "Reset Edits")

    view_scope_label = Observable(with_padding ? "Current view: full_model -> apply full_model" : "Current view: core_only -> apply core_only")
    apply_status = Observable("No changes are applied until 'Apply Changes' is clicked.")
    Label(controls_grid[2, 1:4], view_scope_label, fontsize = 12, color = :darkgreen)
    Label(controls_grid[2, 5:8], apply_status, fontsize = 12, color = :blue)

    applied_cutoff = Observable(parts.current_layer[])
    applied_target = Observable(target_resistivity)
    applied_transition = Observable(Int(blend_previous_percent))
    applied_scope = Observable(replace_scope_used)

    function update_scope_label!()
        if parts.show_full_model[]
            view_scope_label[] = "Current view: full_model -> apply full_model"
        else
            view_scope_label[] = "Current view: core_only -> apply core_only"
        end
    end

    on(parts.show_full_model) do _
        update_scope_label!()
    end

    function parse_integer_field(value::AbstractString, field_name::AbstractString)
        parsed_value = try
            parse(Int, strip(value))
        catch
            apply_status[] = "$(field_name) must be an integer value."
            return nothing
        end

        return parsed_value
    end

    on(parts.current_layer) do layer_idx
        current_layer_label[] = string(layer_idx)
    end

    function apply_controls!()
        target_value_int = parse_integer_field(textbox_text(resistivity_tb), "Target resistivity")
        isnothing(target_value_int) && return nothing

        transition_value = parse_integer_field(textbox_text(transition_tb), "Blend Previous Layer (%)")
        isnothing(transition_value) && return nothing

        if target_value_int <= 0
            apply_status[] = "Target resistivity must be a positive integer."
            return nothing
        end

        if transition_value < 0 || transition_value > 100
            apply_status[] = "Blend Previous Layer (%) must be an integer between 0 and 100."
            return nothing
        end

        cutoff_value = parts.current_layer[]

        scope_value = parts.show_full_model[] ? :full_model : :core_only
        target_value = Float64(target_value_int)

        A_modified .= A_original
        _, cutoff_used, modification_done, scope_used = modify_deep_resistivity!(
            A_modified,
            M;
            cutoff_layer = cutoff_value,
            target_resistivity = target_value,
            blend_previous_percent = transition_value,
            replace_scope = scope_value
        )

        refresh_viewer_data!(parts, A_modified, M;
            log10_scale = log10_scale,
            max_depth = max_depth,
            pad_tol = pad_tol
        )

        applied_cutoff[] = cutoff_used
        applied_target[] = target_value
        applied_transition[] = transition_value
        applied_scope[] = scope_used
        set_textbox_text!(resistivity_tb, string(target_value_int))
        set_textbox_text!(transition_tb, string(transition_value))

        if modification_done
            apply_status[] = "Applied layer $(cutoff_used), ρ=$(round(target_value, digits = 3)) Ω·m, blend_previous=$(transition_value)% in log10 space, scope=$(scope_used)."
        else
            apply_status[] = "No modification applied. Check cutoff layer."
        end

        return nothing
    end

    on(btn_apply.clicks) do _
        apply_controls!()
    end

    on(btn_reset_edits.clicks) do _
        A_modified .= A_original
        refresh_viewer_data!(parts, A_modified, M;
            log10_scale = log10_scale,
            max_depth = max_depth,
            pad_tol = pad_tol
        )

        applied_cutoff[] = parts.current_layer[]
        applied_target[] = target_resistivity
        applied_transition[] = Int(blend_previous_percent)
        applied_scope[] = parts.show_full_model[] ? :full_model : :core_only
        set_textbox_text!(resistivity_tb, string(Int(round(target_resistivity))))
        set_textbox_text!(transition_tb, string(blend_previous_percent))
        save_label[] = ""
        apply_status[] = "Edits reset. Viewer restored to the original model."
    end

    update_scope_label!()

    on(btn_save.clicks) do _
        cutoff_tag = "layer$(applied_cutoff[])"
        rho_tag = "rho$(Int(round(applied_target[])))"
        trans_tag = "blendprev$(applied_transition[])"
        scope_tag = applied_scope[] == :core_only ? "coreonly" : "withpadding"
        outname = splitext(model_file)[1] * "_modified_$(cutoff_tag)_$(rho_tag)_$(trans_tag)_$(scope_tag).rho"
        write_modified_model_with_metadata(outname, M, A_modified;
            cutoff_layer = applied_cutoff[],
            target_resistivity = applied_target[],
            blend_previous_percent = applied_transition[],
            replace_scope = applied_scope[]
        )
        save_label[] = "Saved: $(basename(outname))"
    end

    println("\nViewer ready!")
    if modification_applied
        println("  Cutoff at layer $(cutoff_layer_used)")
        println("  Replace scope: $(replace_scope_used)")
    else
        println("  No modification applied yet.")
    end
    println("  Use the current view layer and current core/full view, then click 'Apply Changes'.")
    println("  Click 'Save Model' to write the modified model to disk.")

    if interactive
        screen = GLMakie.Screen(fig.scene)
        println("\nClose the figure window to exit...")
        wait(screen)
    end

    return fig, parts, (M = M, A_modified = A_modified)
end

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
                     layers_above=2, layers_below=2, transition_layers=1,
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

function EditModelByDrawing(model_file::AbstractString;
    replacement_resistivity::Real = 10000.0,
    layers_above::Int = 2,
    layers_below::Int = 2,
    transition_layers::Int = 1,
    apply_to_all_depths::Bool = false,
    depth_range = (0.0, 50000.0),
    log10_scale::Bool = true,
    colormap = Reverse(:turbo),
    with_padding::Bool = true,
    max_depth::Union{Nothing, Real} = nothing,
    resistivity_range = (0.0, 4.0),
    show_grid::Bool = true,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.3,
    pad_tol::Real = 0.5,
    viewer_figsize = (1200, 950),
    interactive::Bool = true)

    isempty(model_file) && error("model_file is required")
    isfile(model_file) || error("Model file not found: $model_file")

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
        layers_above=layers_above,
        layers_below=layers_below,
        transition_layers=transition_layers,
        apply_to_all_depths=apply_to_all_depths,
        depth_range=depth_range,
        log10scale=log10_scale,
        cmap=colormap,
        figsize=viewer_figsize,
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

    if interactive
        screen = GLMakie.Screen(fig.scene)
        println("\nClose the figure window to exit...")
        wait(screen)
    end

    return fig, parts
end

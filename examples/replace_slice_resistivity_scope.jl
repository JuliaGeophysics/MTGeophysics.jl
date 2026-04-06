# Example: replace deep resistivity values below a cutoff layer with control over padding edits.
# Author: @pankajkmishra
# This script applies bulk edits to deeper layers and lets you choose whether the edit
# is applied only to the core model or to the full model including lateral padding.

using Pkg

using GLMakie
using Statistics
using Dates

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "CoreUtils3D.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))

model_file = joinpath(@__DIR__, "Cascadia", "cascad_half_inverse.ws")

target_resistivity = 1000.0
blend_previous_percent = 0
replace_scope = :core_only

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

function depth_slice_viewer(
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

function main()
    if !isfile(model_file)
        println("ERROR: Model file not found: $model_file")
        return nothing, nothing
    end

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

    fig, parts = depth_slice_viewer(M_modified;
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

    screen = display(fig)
    println("\nClose the figure window to exit...")
    wait(screen)

    return fig, parts, (M = M, A_modified = A_modified)
end

fig, parts, data = main()
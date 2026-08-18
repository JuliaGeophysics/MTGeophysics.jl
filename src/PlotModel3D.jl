# ---------- Shapefile helper functions ----------

function _is_xy_tuple(c)
    c isa Tuple && length(c) >= 2 && c[1] isa Real && c[2] isa Real
end

function _is_xy_vector(c)
    c isa AbstractVector && length(c) >= 2 && c[1] isa Real && c[2] isa Real
end

# ---------- Shapefile CRS detection and reprojection ----------
# Provided by src/ShapefileOverlay.jl: detect_shapefile_crs,
# shapefile_coord_transform, load_shapefile_geometries, prepare_shapefiles.

function _plot_coords_recursive!(ax, coords;
    coord_transform = (x, y) -> (x, y),
    point_color = :black,
    line_color  = :black,
    point_size  = 7,
    line_width  = 1.5,
    alpha       = 1.0)

    if _is_xy_tuple(coords) || _is_xy_vector(coords)
        x, y = coord_transform(Float64(coords[1]), Float64(coords[2]))
        scatter!(ax, [x], [y], color = (point_color, alpha), markersize = point_size)
        return 1
    end

    if coords isa AbstractVector
        isempty(coords) && return 0
        first_item = first(coords)

        if _is_xy_tuple(first_item) || _is_xy_vector(first_item)
            xraw = Float64[p[1] for p in coords if (_is_xy_tuple(p) || _is_xy_vector(p))]
            yraw = Float64[p[2] for p in coords if (_is_xy_tuple(p) || _is_xy_vector(p))]
            xy = [coord_transform(xraw[i], yraw[i]) for i in eachindex(xraw)]
            xs = Float64[p[1] for p in xy]
            ys = Float64[p[2] for p in xy]
            npts = length(xs)
            if npts == 1
                scatter!(ax, xs, ys, color = (point_color, alpha), markersize = point_size)
            elseif npts > 1
                lines!(ax, xs, ys, color = (line_color, alpha), linewidth = line_width)
            end
            return npts
        else
            count = 0
            for part in coords
                count += _plot_coords_recursive!(ax, part;
                    coord_transform = coord_transform,
                    point_color = point_color,
                    line_color  = line_color,
                    point_size  = point_size,
                    line_width  = line_width,
                    alpha       = alpha)
            end
            return count
        end
    end

    return 0
end

function _plot_geometry!(ax, geom;
    coord_transform = (x, y) -> (x, y),
    point_color = :black,
    line_color  = :black,
    point_size  = 7,
    line_width  = 1.5,
    alpha       = 1.0)

    coords = try
        GeoInterface.coordinates(geom)
    catch
        nothing
    end
    isnothing(coords) && return 0

    return _plot_coords_recursive!(ax, coords;
        coord_transform = coord_transform,
        point_color = point_color,
        line_color  = line_color,
        point_size  = point_size,
        line_width  = line_width,
        alpha       = alpha)
end

function _draw_all_shapefiles!(ax, loaded_shapefiles)
    total = 0
    for shp in loaded_shapefiles
        ct = shp.coord_transform
        for geom in shp.geoms
            total += _plot_geometry!(ax, geom;
                coord_transform = ct,
                point_color = shp.color,
                line_color  = shp.color,
                point_size  = shp.point_size,
                line_width  = shp.line_width,
                alpha       = shp.alpha)
        end
    end
    return total
end

# ---------- Main viewer function ----------

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
    axis_aspect::Union{Nothing, Real} = nothing,
    geographic::Bool = false,
    loaded_shapefiles = [],
    export_dpi::Int = 3,
    export_figsize = (1100, 900),
    custom_extent = nothing,
    export_crs::String = "model",
    gis_output_dir::String = "",
    qml_source_path::String = "",
    lyrx_source_path::String = "",
    site_x::Vector{Float64} = Float64[],
    site_y::Vector{Float64} = Float64[]
)
    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : copy(M.A)

    ix_full = 1:length(x_all)
    iy_full = 1:length(y_all)
    ix_core, iy_core = _core_indices_or_extent(x_all, y_all; pad_tol = pad_tol, extent = custom_extent)

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
        title = "",
        xtickformat = _plain_tickformat,
        ytickformat = _plain_tickformat
    )

    current_layer = Observable(1)
    slice_data = Observable(current_R[:, :, 1]')

    current_colormap = Observable(cmap)

    hm = heatmap!(ax, current_y_edges[], current_x_edges[], slice_data,
                  colormap = current_colormap,
                  colorrange = (cmin, cmax))

    isempty(site_x) || scatter!(ax, site_x, site_y; color = :black, marker = :circle, markersize = 4)

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

    # Draw shapefile overlays on the interactive axis
    _draw_all_shapefiles!(ax, loaded_shapefiles)

    # Constrain initial view to the model extent (shapefiles may extend beyond)
    xe_init = current_x_edges[]
    ye_init = current_y_edges[]
    limits!(ax, ye_init[1], ye_init[end], xe_init[1], xe_init[end])

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
    btn_last  = Button(button_grid[1, 2], label = "Last >>|")

    btn_label  = Observable(withPadding ? "Show Core Model" : "Show Full Model")
    btn_toggle = Button(button_grid[1, 3], label = btn_label)

    btn_reset  = Button(button_grid[1, 4], label = "Reset Zoom")
    btn_export = Button(button_grid[1, 5], label = "Export Figure")
    btn_export_gis = Button(button_grid[1, 6], label = "Export GIS")

    info_grid = fig[4, 1:2] = GridLayout()
    depth_info = Observable("Layer depth range: 0 - $(round(layer_depths[1], digits=1)) m | Cell thickness: $(round(diff(z_edges)[1], digits=1)) m")
    Label(info_grid[1, 1], depth_info, fontsize = 12)

    n_shp = length(loaded_shapefiles)
    full_view_str = n_shp > 0 ? "View: Full Model | Shapefiles loaded: $n_shp" : "View: Full Model (with padding)"
    core_view_str = n_shp > 0 ? "View: Core Model | Shapefiles loaded: $n_shp" : "View: Core Model (no padding)"
    view_info = Observable(withPadding ? full_view_str : core_view_str)
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
            view_info[] = full_view_str
        else
            current_x_edges[] = x_edges_core
            current_y_edges[] = y_edges_core
            current_R_ref[] = R_core
            btn_label[] = "Show Full Model"
            view_info[] = core_view_str
        end

        layer_idx = current_layer[]
        slice_data[] = current_R_ref[][:, :, layer_idx]'

        empty!(ax)
        hm = heatmap!(ax, current_y_edges[], current_x_edges[], slice_data,
                      colormap = current_colormap,
                      colorrange = (cmin, cmax))

        isempty(site_x) || scatter!(ax, site_x, site_y; color = :black, marker = :circle, markersize = 4)

        draw_grid!(ax, current_x_edges[], current_y_edges[])

        # Redraw shapefiles after empty!()
        _draw_all_shapefiles!(ax, loaded_shapefiles)

        # Recompute the distance-faithful aspect for the new extent (degrees axes)
        if geographic
            asp, _ = _distance_based_aspect(collect(current_x_edges[]), collect(current_y_edges[]))
            ax.aspect = AxisAspect(Float64(asp))
        end

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

        export_fig = Figure(size = export_figsize, fontsize = 16)

        depth_str = if depth_bottom < 1000
            "Depth: $(round(depth_top, digits=1)) - $(round(depth_bottom, digits=1)) m"
        else
            "Depth: $(round(depth_top/1000, digits=2)) - $(round(depth_bottom/1000, digits=2)) km"
        end
        Label(export_fig[1, 1:2], "Layer $layer_idx | $depth_str",
              fontsize = 20, font = :bold)

        export_aspect = ax_aspect
        if geographic
            asp, _ = _distance_based_aspect(collect(xe), collect(ye))
            export_aspect = AxisAspect(Float64(asp))
        end

        export_ax = Axis(export_fig[2, 1],
            xlabel = axis_xlabel,
            ylabel = axis_ylabel,
            aspect = export_aspect,
            xlabelsize = 16,
            ylabelsize = 16,
            xticklabelsize = 12,
            yticklabelsize = 12,
            xtickformat = _plain_tickformat,
            ytickformat = _plain_tickformat
        )

        export_hm = heatmap!(export_ax, ye, xe, data,
                            colormap = cmap_val,
                            colorrange = (cmin, cmax))

        isempty(site_x) || scatter!(export_ax, site_x, site_y; color = :black, marker = :circle, markersize = 4)

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

        # Shapefile overlay on export
        _draw_all_shapefiles!(export_ax, loaded_shapefiles)

        # Clip export axes to the current view extent
        limits!(export_ax, ye[1], ye[end], xe[1], xe[end])

        cb_lbl = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
        Colorbar(export_fig[2, 2], export_hm, label = cb_lbl, labelsize = 16, ticklabelsize = 12)

        colsize!(export_fig.layout, 2, Relative(0.05))

        view_type = show_full_model[] ? "full" : "core"
        timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        filename = "$(model_name)_XY_layer$(layer_idx)_$(view_type)_$(timestamp).png"

        save(filename, export_fig, px_per_unit = export_dpi)

        w_px = round(Int, export_figsize[1] * export_dpi)
        h_px = round(Int, export_figsize[2] * export_dpi)

        println("Figure exported: $filename")
        println("  Resolution: $(w_px) × $(h_px) pixels")
        println("  Layer: $layer_idx, $depth_str")

        return filename
    end

    on(btn_export.clicks) do _
        export_figure()
    end

    function export_gis_all()
        # Determine which view to export
        xe = current_x_edges[]
        ye = current_y_edges[]
        R_export = current_R_ref[]

        gis_dir = isempty(gis_output_dir) ? joinpath(pwd(), "$(model_name)-GIS") : gis_output_dir

        _export_all_depth_slices_gis(
            output_dir     = gis_dir,
            model_name     = model_name,
            x_edges        = xe,
            y_edges        = ye,
            R              = R_export,
            z_edges        = z_edges,
            crs            = export_crs,
            log10scale     = log10scale,
            qml_source_path = qml_source_path,
            lyrx_source_path = lyrx_source_path
        )

        view_info[] = "GIS exported: $(basename(gis_dir))/ ($(size(R_export, 3)) layers)"
    end

    on(btn_export_gis.clicks) do _
        try
            export_gis_all()
        catch err
            view_info[] = "GIS export failed: $(err)"
            println("GIS export error: ", err)
        end
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

# ---------- Main ----------

# ---------- Format slider value for display ----------

function _format_slider_val(val::Real, unit::AbstractString)
    if unit == "°"
        return "$(round(val, digits=5))°"
    else
        if abs(val) >= 1_000_000
            return "$(round(val/1000, digits=1)) km"
        elseif abs(val) >= 1000
            return "$(round(val, digits=0)) m"
        else
            return "$(round(val, digits=1)) m"
        end
    end
end

# ---------- Main viewer function ----------

function xz_slice_viewer(
    M;
    model_name::String = "model",
    log10scale::Bool = true,
    cmap = :turbo,
    figsize = (1100, 950),
    withPadding::Bool = false,
    max_depth::Union{Nothing, Real} = nothing,
    pad_tol::Real = 0.2,
    resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}} = nothing,
    show_grid::Bool = false,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.3,
    horiz_label::String = "X (m)",
    slider_label::String = "Y",
    slider_unit::String = "m",
    export_dpi::Int = 3,
    export_figsize = (1100, 900)
)
    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : copy(M.A)

    ix_full = 1:length(x_all)
    iy_full = 1:length(y_all)
    ix_core = core_indices(x_all; tol = pad_tol)
    iy_core = core_indices(y_all; tol = pad_tol)

    if isnothing(max_depth)
        kz = 1:length(z_all)
    else
        kz = z_indices_for_max_depth(z_all, float(max_depth))
    end
    z = z_all[kz]

    x_full = x_all[ix_full]
    y_full = y_all[iy_full]
    x_core = x_all[ix_core]
    y_core = y_all[iy_core]

    R_full = A_all[ix_full, iy_full, kz]
    R_core = A_all[ix_core, iy_core, kz]

    x_edges_full = edges_from_centers(x_full)
    x_edges_core = edges_from_centers(x_core)
    z_edges = edges_from_centers(z)

    # Negate z for plotting so depth increases downward (origin top-left)
    neg_z_edges = -z_edges

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

    current_R = withPadding ? R_full : R_core
    current_y = withPadding ? y_full : y_core

    fig = Figure(size = figsize)

    n_y = length(current_y)
    first_y_str = _format_slider_val(current_y[1], slider_unit)
    title_str = Observable("XZ Slice 1 / $(n_y) | $(slider_label) = $(first_y_str)")
    Label(fig[0, 1:2], title_str, fontsize = 18, font = :bold)

    # Tick formatter: show positive depth values on the negated axis
    depth_formatter(values) = [string(round(abs(v), sigdigits=4)) for v in values]

    ax = Axis(fig[1, 1],
        xlabel = horiz_label,
        ylabel = "Depth (m)",
        title = "",
        ytickformat = depth_formatter
    )

    current_slice_idx = Observable(1)
    slice_data = Observable(current_R[:, 1, :])

    current_colormap = Observable(cmap)

    hm = heatmap!(ax, current_x_edges[], neg_z_edges, slice_data,
                  colormap = current_colormap,
                  colorrange = (cmin, cmax))

    grid_plots = Ref{Vector{Any}}([])

    function draw_grid!(ax, x_edges)
        for p in grid_plots[]
            delete!(ax, p)
        end
        grid_plots[] = []

        if show_grid
            for ze in neg_z_edges
                p = lines!(ax, [x_edges[1], x_edges[end]], [ze, ze],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
                push!(grid_plots[], p)
            end
            for xe in x_edges
                p = lines!(ax, [xe, xe], [neg_z_edges[end], neg_z_edges[1]],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
                push!(grid_plots[], p)
            end
        end
    end

    draw_grid!(ax, current_x_edges[])

    # Constrain initial view to the model extent
    xe_init = current_x_edges[]
    limits!(ax, xe_init[1], xe_init[end], neg_z_edges[end], neg_z_edges[1])

    cb_label = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
    Colorbar(fig[1, 2], hm, label = cb_label)

    slider_grid = fig[2, 1:2] = GridLayout()

    btn_prev = Button(slider_grid[1, 1], label = "<< Prev")
    Label(slider_grid[1, 2], "$(slider_label) Slice:", fontsize = 14)
    sl = Slider(slider_grid[1, 3], range = 1:n_y, startvalue = 1, width = 400)
    layer_label = Observable("1 / $(n_y)")
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
    depth_info = Observable("$(slider_label) position: $(first_y_str)")
    Label(info_grid[1, 1], depth_info, fontsize = 12)

    view_info = Observable(withPadding ? "View: Full Model (with padding)" : "View: Core Model (no padding)")
    Label(info_grid[1, 2], view_info, fontsize = 12, color = :blue)

    current_R_ref = Ref(current_R)
    current_y_ref = Ref(current_y)
    n_y_ref = Ref(n_y)

    function reset_zoom!()
        xe = current_x_edges[]
        limits!(ax, xe[1], xe[end], neg_z_edges[end], neg_z_edges[1])
    end

    function update_view!(show_full::Bool)
        if show_full
            current_x_edges[] = x_edges_full
            current_R_ref[] = R_full
            current_y_ref[] = y_full
            n_y_ref[] = length(y_full)
            btn_label[] = "Show Core Model"
            view_info[] = "View: Full Model (with padding)"
        else
            current_x_edges[] = x_edges_core
            current_R_ref[] = R_core
            current_y_ref[] = y_core
            n_y_ref[] = length(y_core)
            btn_label[] = "Show Full Model"
            view_info[] = "View: Core Model (no padding)"
        end

        set_close_to!(sl, clamp(current_slice_idx[], 1, n_y_ref[]))

        slice_idx = clamp(current_slice_idx[], 1, n_y_ref[])
        current_slice_idx[] = slice_idx
        slice_data[] = current_R_ref[][:, slice_idx, :]

        empty!(ax)
        hm = heatmap!(ax, current_x_edges[], neg_z_edges, slice_data,
                      colormap = current_colormap,
                      colorrange = (cmin, cmax))

        draw_grid!(ax, current_x_edges[])

        reset_zoom!()
    end

    function update_slice(slice_idx)
        n = n_y_ref[]
        slice_idx = clamp(slice_idx, 1, n)
        current_slice_idx[] = slice_idx
        slice_data[] = current_R_ref[][:, slice_idx, :]

        layer_label[] = "$slice_idx / $(n)"

        y_pos = current_y_ref[][slice_idx]
        y_str = _format_slider_val(y_pos, slider_unit)

        title_str[] = "XZ Slice $slice_idx / $(n) | $(slider_label) = $(y_str)"
        depth_info[] = "$(slider_label) position: $(y_str)"
    end

    on(sl.value) do val
        update_slice(val)
    end

    on(btn_prev.clicks) do _
        set_close_to!(sl, max(1, sl.value[] - 1))
    end

    on(btn_next.clicks) do _
        set_close_to!(sl, min(n_y_ref[], sl.value[] + 1))
    end

    on(btn_first.clicks) do _
        set_close_to!(sl, 1)
    end

    on(btn_last.clicks) do _
        set_close_to!(sl, n_y_ref[])
    end

    on(btn_toggle.clicks) do _
        show_full_model[] = !show_full_model[]
        update_view!(show_full_model[])
    end

    on(btn_reset.clicks) do _
        reset_zoom!()
    end

    function export_figure()
        slice_idx = current_slice_idx[]
        y_pos = current_y_ref[][slice_idx]
        y_str = _format_slider_val(y_pos, slider_unit)

        xe = current_x_edges[]
        data = slice_data[]
        cmap_val = current_colormap[]

        export_fig = Figure(size = export_figsize, fontsize = 16)

        Label(export_fig[1, 1:2], "XZ Slice $slice_idx | $(slider_label) = $(y_str)",
              fontsize = 20, font = :bold)

        export_ax = Axis(export_fig[2, 1],
            xlabel = horiz_label,
            ylabel = "Depth (m)",
            xlabelsize = 16,
            ylabelsize = 16,
            xticklabelsize = 12,
            yticklabelsize = 12,
            ytickformat = depth_formatter
        )

        export_hm = heatmap!(export_ax, xe, neg_z_edges, data,
                            colormap = cmap_val,
                            colorrange = (cmin, cmax))

        if show_grid
            for ze in neg_z_edges
                lines!(export_ax, [xe[1], xe[end]], [ze, ze],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
            end
            for x_edge in xe
                lines!(export_ax, [x_edge, x_edge], [neg_z_edges[end], neg_z_edges[1]],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
            end
        end

        # Clip export to current model extent
        limits!(export_ax, xe[1], xe[end], neg_z_edges[end], neg_z_edges[1])

        cb_lbl = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
        Colorbar(export_fig[2, 2], export_hm, label = cb_lbl, labelsize = 16, ticklabelsize = 12)

        colsize!(export_fig.layout, 2, Relative(0.05))

        view_type = show_full_model[] ? "full" : "core"
        timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        filename = "$(model_name)_XZ_slice$(slice_idx)_$(view_type)_$(timestamp).png"

        save(filename, export_fig, px_per_unit = export_dpi)

        w_px = round(Int, export_figsize[1] * export_dpi)
        h_px = round(Int, export_figsize[2] * export_dpi)

        println("Figure exported: $filename")
        println("  Resolution: $(w_px) × $(h_px) pixels")
        println("  XZ Slice: $slice_idx, $(slider_label) = $(y_str)")

        return filename
    end

    on(btn_export.clicks) do _
        export_figure()
    end

    return fig, (
        ax = ax,
        heatmap = hm,
        slider = sl,
        current_slice_idx = current_slice_idx,
        slice_data = slice_data,
        show_full_model = show_full_model,
        R_full = R_full,
        R_core = R_core,
        x_full = x_full, y_full = y_full,
        x_core = x_core, y_core = y_core,
        z = z,
        x_edges_full = x_edges_full,
        x_edges_core = x_edges_core,
        z_edges = z_edges,
        colorrange = (cmin = cmin, cmax = cmax),
        reset_zoom! = reset_zoom!
    )
end

# ---------- Main viewer function ----------

function yz_slice_viewer(
    M;
    model_name::String = "model",
    log10scale::Bool = true,
    cmap = :turbo,
    figsize = (1100, 950),
    withPadding::Bool = false,
    max_depth::Union{Nothing, Real} = nothing,
    pad_tol::Real = 0.2,
    resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}} = nothing,
    show_grid::Bool = false,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.3,
    horiz_label::String = "Y (m)",
    slider_label::String = "X",
    slider_unit::String = "m",
    export_dpi::Int = 3,
    export_figsize = (1100, 900)
)
    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : copy(M.A)

    ix_full = 1:length(x_all)
    iy_full = 1:length(y_all)
    ix_core = core_indices(x_all; tol = pad_tol)
    iy_core = core_indices(y_all; tol = pad_tol)

    if isnothing(max_depth)
        kz = 1:length(z_all)
    else
        kz = z_indices_for_max_depth(z_all, float(max_depth))
    end
    z = z_all[kz]

    x_full = x_all[ix_full]
    y_full = y_all[iy_full]
    x_core = x_all[ix_core]
    y_core = y_all[iy_core]

    R_full = A_all[ix_full, iy_full, kz]
    R_core = A_all[ix_core, iy_core, kz]

    y_edges_full = edges_from_centers(y_full)
    y_edges_core = edges_from_centers(y_core)
    z_edges = edges_from_centers(z)

    # Negate z for plotting so depth increases downward (origin top-left)
    neg_z_edges = -z_edges

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

    current_y_edges = Observable(withPadding ? y_edges_full : y_edges_core)

    current_R = withPadding ? R_full : R_core
    current_x = withPadding ? x_full : x_core

    fig = Figure(size = figsize)

    n_x = length(current_x)
    first_x_str = _format_slider_val(current_x[1], slider_unit)
    title_str = Observable("YZ Slice 1 / $(n_x) | $(slider_label) = $(first_x_str)")
    Label(fig[0, 1:2], title_str, fontsize = 18, font = :bold)

    # Tick formatter: show positive depth values on the negated axis
    depth_formatter(values) = [string(round(abs(v), sigdigits=4)) for v in values]

    ax = Axis(fig[1, 1],
        xlabel = horiz_label,
        ylabel = "Depth (m)",
        title = "",
        ytickformat = depth_formatter
    )

    current_slice_idx = Observable(1)
    slice_data = Observable(current_R[1, :, :])

    current_colormap = Observable(cmap)

    hm = heatmap!(ax, current_y_edges[], neg_z_edges, slice_data,
                  colormap = current_colormap,
                  colorrange = (cmin, cmax))

    grid_plots = Ref{Vector{Any}}([])

    function draw_grid!(ax, y_edges)
        for p in grid_plots[]
            delete!(ax, p)
        end
        grid_plots[] = []

        if show_grid
            for ze in neg_z_edges
                p = lines!(ax, [y_edges[1], y_edges[end]], [ze, ze],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
                push!(grid_plots[], p)
            end
            for ye in y_edges
                p = lines!(ax, [ye, ye], [neg_z_edges[end], neg_z_edges[1]],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
                push!(grid_plots[], p)
            end
        end
    end

    draw_grid!(ax, current_y_edges[])

    # Constrain initial view to the model extent
    ye_init = current_y_edges[]
    limits!(ax, ye_init[1], ye_init[end], neg_z_edges[end], neg_z_edges[1])

    cb_label = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
    Colorbar(fig[1, 2], hm, label = cb_label)

    slider_grid = fig[2, 1:2] = GridLayout()

    btn_prev = Button(slider_grid[1, 1], label = "<< Prev")
    Label(slider_grid[1, 2], "$(slider_label) Slice:", fontsize = 14)
    sl = Slider(slider_grid[1, 3], range = 1:n_x, startvalue = 1, width = 400)
    layer_label = Observable("1 / $(n_x)")
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
    depth_info = Observable("$(slider_label) position: $(first_x_str)")
    Label(info_grid[1, 1], depth_info, fontsize = 12)

    view_info = Observable(withPadding ? "View: Full Model (with padding)" : "View: Core Model (no padding)")
    Label(info_grid[1, 2], view_info, fontsize = 12, color = :blue)

    current_R_ref = Ref(current_R)
    current_x_ref = Ref(current_x)
    n_x_ref = Ref(n_x)

    function reset_zoom!()
        ye = current_y_edges[]
        limits!(ax, ye[1], ye[end], neg_z_edges[end], neg_z_edges[1])
    end

    function update_view!(show_full::Bool)
        if show_full
            current_y_edges[] = y_edges_full
            current_R_ref[] = R_full
            current_x_ref[] = x_full
            n_x_ref[] = length(x_full)
            btn_label[] = "Show Core Model"
            view_info[] = "View: Full Model (with padding)"
        else
            current_y_edges[] = y_edges_core
            current_R_ref[] = R_core
            current_x_ref[] = x_core
            n_x_ref[] = length(x_core)
            btn_label[] = "Show Full Model"
            view_info[] = "View: Core Model (no padding)"
        end

        set_close_to!(sl, clamp(current_slice_idx[], 1, n_x_ref[]))

        slice_idx = clamp(current_slice_idx[], 1, n_x_ref[])
        current_slice_idx[] = slice_idx
        slice_data[] = current_R_ref[][slice_idx, :, :]

        empty!(ax)
        hm = heatmap!(ax, current_y_edges[], neg_z_edges, slice_data,
                      colormap = current_colormap,
                      colorrange = (cmin, cmax))

        draw_grid!(ax, current_y_edges[])

        reset_zoom!()
    end

    function update_slice(slice_idx)
        n = n_x_ref[]
        slice_idx = clamp(slice_idx, 1, n)
        current_slice_idx[] = slice_idx
        slice_data[] = current_R_ref[][slice_idx, :, :]

        layer_label[] = "$slice_idx / $(n)"

        x_pos = current_x_ref[][slice_idx]
        x_str = _format_slider_val(x_pos, slider_unit)

        title_str[] = "YZ Slice $slice_idx / $(n) | $(slider_label) = $(x_str)"
        depth_info[] = "$(slider_label) position: $(x_str)"
    end

    on(sl.value) do val
        update_slice(val)
    end

    on(btn_prev.clicks) do _
        set_close_to!(sl, max(1, sl.value[] - 1))
    end

    on(btn_next.clicks) do _
        set_close_to!(sl, min(n_x_ref[], sl.value[] + 1))
    end

    on(btn_first.clicks) do _
        set_close_to!(sl, 1)
    end

    on(btn_last.clicks) do _
        set_close_to!(sl, n_x_ref[])
    end

    on(btn_toggle.clicks) do _
        show_full_model[] = !show_full_model[]
        update_view!(show_full_model[])
    end

    on(btn_reset.clicks) do _
        reset_zoom!()
    end

    function export_figure()
        slice_idx = current_slice_idx[]
        x_pos = current_x_ref[][slice_idx]
        x_str = _format_slider_val(x_pos, slider_unit)

        ye = current_y_edges[]
        data = slice_data[]
        cmap_val = current_colormap[]

        export_fig = Figure(size = export_figsize, fontsize = 16)

        Label(export_fig[1, 1:2], "YZ Slice $slice_idx | $(slider_label) = $(x_str)",
              fontsize = 20, font = :bold)

        export_ax = Axis(export_fig[2, 1],
            xlabel = horiz_label,
            ylabel = "Depth (m)",
            xlabelsize = 16,
            ylabelsize = 16,
            xticklabelsize = 12,
            yticklabelsize = 12,
            ytickformat = depth_formatter
        )

        export_hm = heatmap!(export_ax, ye, neg_z_edges, data,
                            colormap = cmap_val,
                            colorrange = (cmin, cmax))

        if show_grid
            for ze in neg_z_edges
                lines!(export_ax, [ye[1], ye[end]], [ze, ze],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
            end
            for y_edge in ye
                lines!(export_ax, [y_edge, y_edge], [neg_z_edges[end], neg_z_edges[1]],
                       color = (grid_color, grid_alpha),
                       linewidth = grid_linewidth)
            end
        end

        # Clip export to current model extent
        limits!(export_ax, ye[1], ye[end], neg_z_edges[end], neg_z_edges[1])

        cb_lbl = log10scale ? "log₁₀ ρ (Ω·m)" : "ρ (Ω·m)"
        Colorbar(export_fig[2, 2], export_hm, label = cb_lbl, labelsize = 16, ticklabelsize = 12)

        colsize!(export_fig.layout, 2, Relative(0.05))

        view_type = show_full_model[] ? "full" : "core"
        timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
        filename = "$(model_name)_YZ_slice$(slice_idx)_$(view_type)_$(timestamp).png"

        save(filename, export_fig, px_per_unit = export_dpi)

        w_px = round(Int, export_figsize[1] * export_dpi)
        h_px = round(Int, export_figsize[2] * export_dpi)

        println("Figure exported: $filename")
        println("  Resolution: $(w_px) × $(h_px) pixels")
        println("  YZ Slice: $slice_idx, $(slider_label) = $(x_str)")

        return filename
    end

    on(btn_export.clicks) do _
        export_figure()
    end

    return fig, (
        ax = ax,
        heatmap = hm,
        slider = sl,
        current_slice_idx = current_slice_idx,
        slice_data = slice_data,
        show_full_model = show_full_model,
        R_full = R_full,
        R_core = R_core,
        x_full = x_full, y_full = y_full,
        x_core = x_core, y_core = y_core,
        z = z,
        y_edges_full = y_edges_full,
        y_edges_core = y_edges_core,
        z_edges = z_edges,
        colorrange = (cmin = cmin, cmax = cmax),
        reset_zoom! = reset_zoom!
    )
end

function draw_north_and_scale!(ax;
    xv,
    yv,
    z_fixed::Real,
    target_crs::AbstractString,
    north_axis::Symbol = :y,
    show_north::Bool = true,
    show_scale::Bool = true,
    color = :black,
    line_width::Real = 2.0)

    isempty(xv) && return
    isempty(yv) && return

    xmin, xmax = extrema(xv)
    ymin, ymax = extrema(yv)
    dx = xmax - xmin
    dy = ymax - ymin
    (dx <= 0 || dy <= 0) && return

    z = Float64(z_fixed)
    x_north = xmin + 0.86 * dx
    y_north = ymin + 0.10 * dy
    x_scale = xmin + 0.68 * dx
    y_scale = ymin + 0.10 * dy

    if show_north
        arrow_len = 0.06 * min(dx, dy)
        if north_axis == :x
            x1 = x_north + arrow_len
            y1 = y_north
            ah = 0.025 * min(dx, dy)
            lines!(ax, [x_north, x1], [y_north, y1], [z, z], color = color, linewidth = line_width)
            lines!(ax, [x1 - ah, x1], [y1 - ah, y1], [z, z], color = color, linewidth = line_width)
            lines!(ax, [x1 - ah, x1], [y1 + ah, y1], [z, z], color = color, linewidth = line_width)
            text!(ax, [x1 + 0.02 * dx], [y1], [z], text = ["N"], color = color, fontsize = 16)
        else
            x1 = x_north
            y1 = y_north + arrow_len
            ah = 0.025 * min(dx, dy)
            lines!(ax, [x_north, x1], [y_north, y1], [z, z], color = color, linewidth = line_width)
            lines!(ax, [x1 - ah, x1], [y1 - ah, y1], [z, z], color = color, linewidth = line_width)
            lines!(ax, [x1 + ah, x1], [y1 - ah, y1], [z, z], color = color, linewidth = line_width)
            text!(ax, [x1], [y1 + 0.02 * dy], [z], text = ["N"], color = color, fontsize = 16)
        end
    end

    if show_scale
        target_len = 0.10 * dx
        scale_len = _nice_scale_length(target_len)
        xb0 = x_scale
        xb1 = xb0 + scale_len
        yb = y_scale

        lines!(ax, [xb0, xb1], [yb, yb], [z, z], color = color, linewidth = line_width)

        tick = 0.008 * dy
        lines!(ax, [xb0, xb0], [yb - tick, yb + tick], [z, z], color = color, linewidth = line_width)
        lines!(ax, [xb1, xb1], [yb - tick, yb + tick], [z, z], color = color, linewidth = line_width)

        unit_label = uppercase(strip(target_crs)) == "EPSG:4326" ? "deg" : "m"
        label_val = unit_label == "m" && scale_len >= 1000 ? "$(round(scale_len / 1000; digits = 2)) km" : "$(round(scale_len; sigdigits = 3)) $unit_label"
        text!(ax, [0.5 * (xb0 + xb1)], [yb + 0.02 * dy], [z], text = [label_val], color = color, fontsize = 12)
    end
end

function plot_shapefile_on_3d!(ax, shapefile_path;
    z_fixed::Real = 0.0,
    point_color = :black,
    line_color = :black,
    point_size = 7,
    line_width = 1.5,
    auto_reproject_to_wgs84 = true,
    post_transform = (x, y) -> (x, y),
    xlim::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    ylim::Union{Nothing, Tuple{<:Real, <:Real}} = nothing)
    if isnothing(shapefile_path) || !isfile(shapefile_path)
        isnothing(shapefile_path) || @warn "Shapefile not found: $shapefile_path"
        return 0
    end
    table = Shapefile.Table(shapefile_path)
    rows = collect(table)
    crs_type, _ = detect_shapefile_crs(shapefile_path)
    coord_transform, _, _ = _make_coord_transform(shapefile_path, crs_type; auto_reproject_to_wgs84 = auto_reproject_to_wgs84)
    function _is_xy(c)
        return (c isa Tuple || c isa AbstractVector) && length(c) >= 2 && c[1] isa Real && c[2] isa Real
    end

    xlo, xhi = isnothing(xlim) ? (-Inf, Inf) : (min(Float64(xlim[1]), Float64(xlim[2])), max(Float64(xlim[1]), Float64(xlim[2])))
    ylo, yhi = isnothing(ylim) ? (-Inf, Inf) : (min(Float64(ylim[1]), Float64(ylim[2])), max(Float64(ylim[1]), Float64(ylim[2])))
    inside = (x, y) -> (x >= xlo && x <= xhi && y >= ylo && y <= yhi)

    function plot_coords_recursive!(coords)
        if _is_xy(coords)
            return 0
        elseif coords isa AbstractVector
            isempty(coords) && return 0
            first_item = first(coords)
            if _is_xy(first_item)
                xy = Tuple{Float64, Float64}[]
                for p in coords
                    if _is_xy(p)
                        x0, y0 = coord_transform(Float64(p[1]), Float64(p[2]))
                        x, y = post_transform(x0, y0)
                        push!(xy, (x, y))
                    end
                end

                segments_plotted = 0
                run_x = Float64[]
                run_y = Float64[]
                for (x, y) in xy
                    if inside(x, y)
                        push!(run_x, x)
                        push!(run_y, y)
                    else
                        if length(run_x) >= 2
                            lines!(ax, run_x, run_y, fill(z_fixed, length(run_x)), color = line_color, linewidth = line_width)
                            segments_plotted += 1
                        end
                        empty!(run_x)
                        empty!(run_y)
                    end
                end

                if length(run_x) >= 2
                    lines!(ax, run_x, run_y, fill(z_fixed, length(run_x)), color = line_color, linewidth = line_width)
                    segments_plotted += 1
                end

                return segments_plotted
            else
                count = 0
                for part in coords
                    count += plot_coords_recursive!(part)
                end
                return count
            end
        end
        return 0
    end
    total_points = 0
    for row in rows
        geom = GeoInterface.geometry(row)
        coords = try
            GeoInterface.coordinates(geom)
        catch
            nothing
        end
        if !isnothing(coords)
            total_points += plot_coords_recursive!(coords)
        end
    end
    return total_points
end

function modem_3d_viewer(
    M;
    log10scale::Bool = true,
    cmap = :Spectral,
    figsize = (1800, 920),
    withPadding::Bool = true,
    max_depth::Union{Nothing, Real} = nothing,
    pad_tol::Real = 0.2,
    resistivity_range::Union{Nothing, Tuple{<:Real,<:Real}} = nothing,
    target_crs::AbstractString,
    overlay_transform = (x, y) -> (x, y),
    north_axis::Symbol = :y,
    site_e::Vector{Float64} = Float64[],
    site_n::Vector{Float64} = Float64[],
    shapefile_path = nothing,
    overlay_z_fixed::Real = 0.0,
    overlay_auto_reproject_to_wgs84::Bool = true,
    overlay_point_color = :black,
    overlay_line_color = :black,
    overlay_point_size = 7,
    overlay_line_width = 1.5,
    show_north_arrow::Bool = true,
    show_scale_bar::Bool = true,
    annotation_color = :black,
    annotation_line_width::Real = 2.0,
    default_view_direction = Vec3f(-1.05f0, -0.80f0, 0.72f0),
    default_view_scale = 1.12f0
)

    x_all = M.cx
    y_all = M.cy
    z_all = M.cz
    A_all = log10scale ? log10.(M.A) : M.A
    axis_meta = _axis_metadata_for_crs(target_crs)

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
    !isnothing(shapefile_path) && plot_shapefile_on_3d!(ax.scene, shapefile_path;
        z_fixed = overlay_z_fixed,
        point_color = overlay_point_color,
        line_color = overlay_line_color,
        point_size = overlay_point_size,
        line_width = overlay_line_width,
        auto_reproject_to_wgs84 = overlay_auto_reproject_to_wgs84,
        post_transform = overlay_transform,
        xlim = extrema(current_x[]),
        ylim = extrema(current_y[]))
    draw_north_and_scale!(ax.scene;
        xv = current_x[],
        yv = current_y[],
        z_fixed = overlay_z_fixed,
        target_crs = target_crs,
        north_axis = north_axis,
        show_north = show_north_arrow,
        show_scale = show_scale_bar,
        color = annotation_color,
        line_width = annotation_line_width)
    isempty(site_e) || scatter!(ax.scene, site_e, site_n, fill(Float64(overlay_z_fixed), length(site_e)); color = :black, marker = :circle, markersize = 7)
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
    lbl_yz = Label(controls[1, 6], "$(round(x[sl_yz.value[]], digits=3)) $(axis_meta.unit)", fontsize = 10, color = :gray35)

    Label(controls[1, 7], "Y:", halign = :right, fontsize = 12)
    btn_prev_xz = Button(controls[1, 8], label = "Prev", fontsize = 10)
    sl_xz = Slider(controls[1, 9], range = 1:length(y), startvalue = round(Int, length(y) / 2), width = 180)
    btn_next_xz = Button(controls[1, 10], label = "Next", fontsize = 10)
    tog_xz = Toggle(controls[1, 11], active = true)
    lbl_xz = Label(controls[1, 12], "$(round(y[sl_xz.value[]], digits=3)) $(axis_meta.unit)", fontsize = 10, color = :gray35)

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
    Label(axis_bar[1, 2], "CRS: $(target_crs)", fontsize = 11, color = :gray35, halign = :right)

    rowgap!(controls, 2)
    colgap!(controls, 6)

    function update_axis_info!()
        xv = current_x[]
        yv = current_y[]
        zv = current_z[]
        sx = xv[sl_yz.value[]]
        sy = yv[sl_xz.value[]]
        sz = -zv[sl_xy.value[]]
        axis_info[] = "$(axis_meta.x_name): [$(round(minimum(xv), digits=3)), $(round(maximum(xv), digits=3))] $(axis_meta.unit)   $(axis_meta.y_name): [$(round(minimum(yv), digits=3)), $(round(maximum(yv), digits=3))] $(axis_meta.unit)   Depth: [0, $(round(maximum(-zv), digits=0))] m   |   Slice: $(axis_meta.x_name)=$(round(sx, digits=3)) $(axis_meta.unit), $(axis_meta.y_name)=$(round(sy, digits=3)) $(axis_meta.unit), Depth=$(round(sz, digits=0)) m"
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
        !isnothing(shapefile_path) && plot_shapefile_on_3d!(ax.scene, shapefile_path;
            z_fixed = overlay_z_fixed,
            point_color = overlay_point_color,
            line_color = overlay_line_color,
            point_size = overlay_point_size,
            line_width = overlay_line_width,
            auto_reproject_to_wgs84 = overlay_auto_reproject_to_wgs84,
            post_transform = overlay_transform,
            xlim = extrema(new_x),
            ylim = extrema(new_y))
        draw_north_and_scale!(ax.scene;
            xv = new_x,
            yv = new_y,
            z_fixed = overlay_z_fixed,
            target_crs = target_crs,
            north_axis = north_axis,
            show_north = show_north_arrow,
            show_scale = show_scale_bar,
            color = annotation_color,
            line_width = annotation_line_width)
        isempty(site_e) || scatter!(ax.scene, site_e, site_n, fill(Float64(overlay_z_fixed), length(site_e)); color = :black, marker = :circle, markersize = 7)

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
        lbl_yz.text[] = "$(round(new_x[mid_x], digits=3)) $(axis_meta.unit)"
        lbl_xz.text[] = "$(round(new_y[mid_y], digits=3)) $(axis_meta.unit)"
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

        isempty(site_e) || scatter!(export_ax, site_e, site_n, fill(Float64(overlay_z_fixed), length(site_e)); color = :black, marker = :circle, markersize = 7)

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

function PlotModelXY(model_file::AbstractString, data_file::AbstractString = "";
    crs::AbstractString = "model",
    shapefiles = [],
    log10_scale::Bool = true,
    colormap = :Spectral,
    with_padding::Bool = false,
    max_depth::Union{Nothing, Real} = nothing,
    resistivity_range = (0.0, 4.0),
    show_grid::Bool = true,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.1,
    pad_tol::Real = 0.5,
    custom_extent = nothing,
    viewer_figsize = (1100, 950),
    export_dpi::Int = 3,
    export_figsize = (1100, 900),
    gis_output_dir::AbstractString = "",
    qml_style_file::AbstractString = joinpath(@__DIR__, "MTSlices.qml"),
    lyrx_style_file::AbstractString = joinpath(@__DIR__, "MTSlices.lyrx"),
    interactive::Bool = true)

    isempty(model_file) && error("model_file is required")
    isfile(model_file) || error("Model file not found: $model_file")

    println("Loading ModEM model: $model_file")
    M = load_model_modem(model_file)

    println("Model loaded successfully!")
    println("  Grid size: $(M.nx) × $(M.ny) × $(M.nz)")
    println("  X range: $(minimum(M.x)) to $(maximum(M.x)) m")
    println("  Y range: $(minimum(M.y)) to $(maximum(M.y)) m")
    println("  Z range: $(minimum(M.z)) to $(maximum(M.z)) m")
    println("  Padding cells: $(M.npad)")

    crs_up = uppercase(strip(crs))
    need_data = crs_up != "MODEL"

    d = nothing
    if !isempty(data_file) && isfile(data_file)
        println("\nLoading ModEM data: $data_file")
        d = load_data_modem(data_file)
    elseif need_data
        isempty(data_file) && error("A ModEM data file is required for crs = \"$crs\"; pass data_file, or use crs = \"model\"")
        error("Data file not found: $data_file")
    end

    println("\nCoordinate conversion:")
    crs_extent = _convert_extent_to_crs(custom_extent, crs)
    if !isnothing(custom_extent)
        if isnothing(crs_extent)
            println("  Custom extent: ignored (not supported for this CRS)")
        else
            println("  Custom extent (lat/lon → $(crs)):")
            println("    x-axis: $(crs_extent.min_x) .. $(crs_extent.max_x)")
            println("    y-axis: $(crs_extent.min_y) .. $(crs_extent.max_y)")
        end
    end
    M_crs, xlabel, ylabel, aspect = _build_model_in_crs(M, d, crs, crs_extent;
        pad_tol = pad_tol, with_padding = with_padding)

    loaded_shapefiles = []
    if !isempty(shapefiles)
        println("\nShapefile overlay:")
        loaded_shapefiles = prepare_shapefiles(shapefiles, crs)
    end

    site_x, site_y = _stations_in_crs(d, crs)

    println("\nCreating interactive depth-slice viewer...")

    model_name = splitext(basename(model_file))[1]
    if crs_up != "MODEL"
        model_name *= "_$(replace(lowercase(strip(crs)), ":" => ""))"
    end

    fig, parts = depth_slice_viewer(M_crs;
        model_name      = model_name,
        log10scale      = log10_scale,
        cmap            = colormap,
        figsize         = viewer_figsize,
        withPadding     = with_padding,
        max_depth       = max_depth,
        pad_tol         = pad_tol,
        resistivity_range = resistivity_range,
        show_grid       = show_grid,
        grid_color      = grid_color,
        grid_linewidth  = grid_linewidth,
        grid_alpha      = grid_alpha,
        axis_xlabel     = xlabel,
        axis_ylabel     = ylabel,
        axis_aspect     = aspect,
        geographic      = (crs_up == "EPSG:4326"),
        loaded_shapefiles = loaded_shapefiles,
        export_dpi      = export_dpi,
        export_figsize  = export_figsize,
        custom_extent   = crs_extent,
        export_crs      = crs,
        gis_output_dir  = gis_output_dir,
        qml_source_path = qml_style_file,
        lyrx_source_path = lyrx_style_file,
        site_x          = site_x,
        site_y          = site_y)

    println("\nViewer ready!")
    println("  - Coordinate system: $crs")
    println("  - Use the slider to navigate through depth layers")
    println("  - Use Previous/Next buttons for step-by-step navigation")
    println("  - Click 'Export Figure' for high-resolution PNG export")
    isempty(loaded_shapefiles) || println("  - $(length(loaded_shapefiles)) shapefile(s) overlaid on the depth slice")

    if interactive
        Makie.update_state_before_display!(fig)
        screen = GLMakie.Screen(fig.scene)
        if screen isa GLMakie.Screen && GLMakie.requires_update(screen)
            GLMakie.render_frame(screen)
        end
        println("\nClose the figure window to exit...")
        wait(screen)
    end

    return fig, parts
end

function _plot_model_section(model_file::AbstractString, data_file::AbstractString, crs::AbstractString)
    isempty(model_file) && error("model_file is required")
    isfile(model_file) || error("Model file not found: $model_file")

    println("Loading ModEM model: $model_file")
    M = load_model_modem(model_file)

    println("Model loaded successfully!")
    println("  Grid size: $(M.nx) × $(M.ny) × $(M.nz)")
    println("  X range: $(minimum(M.x)) to $(maximum(M.x)) m")
    println("  Y range: $(minimum(M.y)) to $(maximum(M.y)) m")
    println("  Z range: $(minimum(M.z)) to $(maximum(M.z)) m")
    println("  Padding cells: $(M.npad)")

    crs_up = uppercase(strip(crs))

    d = nothing
    if crs_up != "MODEL"
        isempty(data_file) && error("A ModEM data file is required for crs = \"$crs\"; pass data_file, or use crs = \"model\"")
        isfile(data_file) || error("Data file not found: $data_file")
        println("\nLoading ModEM data for georeferencing: $data_file")
        d = load_data_modem(data_file)
    end

    model_name = splitext(basename(model_file))[1]
    if crs_up != "MODEL"
        model_name *= "_$(replace(lowercase(strip(crs)), ":" => ""))"
    end

    return M, d, model_name
end

function PlotModelXZ(model_file::AbstractString, data_file::AbstractString = "";
    crs::AbstractString = "model",
    log10_scale::Bool = true,
    colormap = :Spectral,
    with_padding::Bool = false,
    max_depth::Union{Nothing, Real} = 200000,
    resistivity_range = (0.0, 4.0),
    show_grid::Bool = false,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.1,
    pad_tol::Real = 0.5,
    viewer_figsize = (1100, 950),
    export_dpi::Int = 3,
    export_figsize = (1100, 900),
    interactive::Bool = true)

    M, d, model_name = _plot_model_section(model_file, data_file, crs)

    println("\nCoordinate conversion:")
    M_crs, horiz_label, slider_label, slider_unit = _build_xz_model_in_crs(M, d, crs)

    println("\nCreating interactive XZ cross-section viewer...")

    fig, parts = xz_slice_viewer(M_crs;
        model_name      = model_name,
        log10scale      = log10_scale,
        cmap            = colormap,
        figsize         = viewer_figsize,
        withPadding     = with_padding,
        max_depth       = max_depth,
        pad_tol         = pad_tol,
        resistivity_range = resistivity_range,
        show_grid       = show_grid,
        grid_color      = grid_color,
        grid_linewidth  = grid_linewidth,
        grid_alpha      = grid_alpha,
        horiz_label     = horiz_label,
        slider_label    = slider_label,
        slider_unit     = slider_unit,
        export_dpi      = export_dpi,
        export_figsize  = export_figsize)

    println("\nViewer ready!")
    println("  - Coordinate system: $crs")
    println("  - Horizontal axis: $horiz_label")
    println("  - Slider axis: $slider_label ($slider_unit)")
    println("  - Use the slider to navigate through $(slider_label) slices (XZ cross-sections)")
    println("  - Use Previous/Next buttons for step-by-step navigation")
    println("  - Click 'Export Figure' for high-resolution PNG export")

    if interactive
        screen = GLMakie.Screen(fig.scene)
        println("\nClose the figure window to exit...")
        wait(screen)
    end

    return fig, parts
end

function PlotModelYZ(model_file::AbstractString, data_file::AbstractString = "";
    crs::AbstractString = "model",
    log10_scale::Bool = true,
    colormap = :Spectral,
    with_padding::Bool = false,
    max_depth::Union{Nothing, Real} = 50000,
    resistivity_range = (0.0, 4.0),
    show_grid::Bool = false,
    grid_color = :black,
    grid_linewidth::Real = 0.5,
    grid_alpha::Real = 0.1,
    pad_tol::Real = 0.5,
    viewer_figsize = (1100, 950),
    export_dpi::Int = 3,
    export_figsize = (1100, 900),
    interactive::Bool = true)

    M, d, model_name = _plot_model_section(model_file, data_file, crs)

    println("\nCoordinate conversion:")
    M_crs, horiz_label, slider_label, slider_unit = _build_yz_model_in_crs(M, d, crs)

    println("\nCreating interactive YZ cross-section viewer...")

    fig, parts = yz_slice_viewer(M_crs;
        model_name      = model_name,
        log10scale      = log10_scale,
        cmap            = colormap,
        figsize         = viewer_figsize,
        withPadding     = with_padding,
        max_depth       = max_depth,
        pad_tol         = pad_tol,
        resistivity_range = resistivity_range,
        show_grid       = show_grid,
        grid_color      = grid_color,
        grid_linewidth  = grid_linewidth,
        grid_alpha      = grid_alpha,
        horiz_label     = horiz_label,
        slider_label    = slider_label,
        slider_unit     = slider_unit,
        export_dpi      = export_dpi,
        export_figsize  = export_figsize)

    println("\nViewer ready!")
    println("  - Coordinate system: $crs")
    println("  - Horizontal axis: $horiz_label")
    println("  - Slider axis: $slider_label ($slider_unit)")
    println("  - Use the slider to navigate through $(slider_label) slices (YZ cross-sections)")
    println("  - Use Previous/Next buttons for step-by-step navigation")
    println("  - Click 'Export Figure' for high-resolution PNG export")

    if interactive
        screen = GLMakie.Screen(fig.scene)
        println("\nClose the figure window to exit...")
        wait(screen)
    end

    return fig, parts
end

function PlotModelXYZ(model_file::AbstractString, data_file::AbstractString;
    crs::AbstractString = "PROJECT",
    shapefile_path = nothing,
    log10_scale::Bool = true,
    colormap = :Spectral,
    resistivity_range = (0.0, 4.0),
    max_depth::Union{Nothing, Real} = 5000.0,
    show_padding::Bool = false,
    pad_tolerance::Real = 0.2,
    viewer_figsize = (1800, 920),
    default_view_direction = Vec3f(-1.05f0, -0.80f0, 0.72f0),
    default_view_scale = 1.12f0,
    overlay_z_fixed::Real = 0.0,
    overlay_auto_reproject_to_wgs84::Bool = true,
    overlay_point_color = :black,
    overlay_line_color = :black,
    overlay_point_size = 7,
    overlay_line_width = 1.5,
    show_north_arrow::Bool = true,
    show_scale_bar::Bool = true,
    annotation_color = :black,
    annotation_line_width::Real = 2.0,
    interactive::Bool = true)

    isempty(model_file) && error("model_file is required")
    isfile(model_file) || error("Model file not found: $model_file")
    isempty(data_file) && error("data_file is required for georeferencing")
    isfile(data_file) || error("Data file not found: $data_file")

    println("Loading ModEM model from: $model_file")
    M = load_model_modem(model_file)
    println("Loading ModEM data for georeference from: $data_file")
    d = load_data_modem(data_file)

    resolved_target_crs = _normalize_target_crs(crs)
    _warn_if_target_crs_looks_unexpected(d, resolved_target_crs)

    if _is_project_local_crs(resolved_target_crs)
        println("Using project-local coordinates from data origin:")
        println("  Coordinate system: local transverse Mercator / survey metres")
    end

    x_target, y_target, lat0, lon0, shiftlat, shiftlon, lat_ref, lon_ref, mismatch_dim_consistent, station_tx, station_ty =
        model_xy_to_target_crs_centers(M, d, resolved_target_crs)
    wgs84_to_target = _resolve_wgs84_to_target_xy_transform(resolved_target_crs, lat0, lon0)
    overlay_transform = (lon, lat) -> begin
        lon_aligned = Float64(lon) + Float64(shiftlon)
        lat_aligned = Float64(lat) + Float64(shiftlat)
        e, n = wgs84_to_target(lon_aligned, lat_aligned)
        return e, n
    end
    north_axis = :y

    M_target = (
        A = permutedims(M.A, (2, 1, 3)),
        cx = x_target,
        cy = y_target,
        cz = M.cz
    )

    println("Model dimensions: $(size(M.A))")
    println("  X cells: $(length(M.cx))")
    println("  Y cells: $(length(M.cy))")
    println("  Z cells: $(length(M.cz))")
    println("Target CRS plotting:")
    println("  Target CRS: $resolved_target_crs")
    axis_meta = _axis_metadata_for_crs(resolved_target_crs)
    println("  Axis convention: X=$(axis_meta.x_name), Y=$(axis_meta.y_name)")
    println("  Axis consistency mismatch: $(round(mismatch_dim_consistent, digits=4))")
    println("  X range: [$(round(minimum(x_target), digits=3)), $(round(maximum(x_target), digits=3))]")
    println("  Y range: [$(round(minimum(y_target), digits=3)), $(round(maximum(y_target), digits=3))]")
    println("  Reference lat/lon: ($(round(lat_ref, digits = 6)), $(round(lon_ref, digits = 6)))")
    println("Model georeference:")
    println("  Origin (lat, lon): ($(round(lat0, digits = 6)), $(round(lon0, digits = 6)))")
    println("  Data alignment shift: Δlat=$(shiftlat), Δlon=$(shiftlon)")

    println("\nLaunching 3D viewer...")
    println("  Log10 scale: $log10_scale")
    println("  Colormap: $colormap")
    println("  Resistivity range: $resistivity_range")
    println("  Show padding: $show_padding")
    if !isnothing(max_depth)
        println("  Max depth: $max_depth m")
    end

    fig, parts = modem_3d_viewer(M_target;
        log10scale = log10_scale,
        cmap = colormap,
        figsize = viewer_figsize,
        withPadding = show_padding,
        max_depth = max_depth,
        pad_tol = pad_tolerance,
        resistivity_range = resistivity_range,
        target_crs = resolved_target_crs,
        overlay_transform = overlay_transform,
        north_axis = north_axis,
        site_e = collect(Float64.(station_tx)),
        site_n = collect(Float64.(station_ty)),
        shapefile_path = shapefile_path,
        overlay_z_fixed = overlay_z_fixed,
        overlay_auto_reproject_to_wgs84 = overlay_auto_reproject_to_wgs84,
        overlay_point_color = overlay_point_color,
        overlay_line_color = overlay_line_color,
        overlay_point_size = overlay_point_size,
        overlay_line_width = overlay_line_width,
        show_north_arrow = show_north_arrow,
        show_scale_bar = show_scale_bar,
        annotation_color = annotation_color,
        annotation_line_width = annotation_line_width,
        default_view_direction = default_view_direction,
        default_view_scale = default_view_scale)

    println("\n3D Viewer Controls:")
    println("  - Drag: Rotate view")
    println("  - Scroll: Zoom in/out")
    println("  - Right-drag: Pan")
    println("  - Prev/Next + sliders: Step or move slice planes")
    println("  - Toggles: Show/hide planes")
    println("  - 'Show Core/Full': Toggle padding cells")
    println("  - 'Reset View': Reset camera angle")
    println("  - 'Export Figure': Save high-quality PNG")

    if interactive
        screen = GLMakie.Screen(fig.scene)
        println("\nViewer is open. Close the window to exit.")
        wait(screen)
    end

    return fig, parts
end

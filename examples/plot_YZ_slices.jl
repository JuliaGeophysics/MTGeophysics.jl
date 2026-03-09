# Example: interactive YZ cross-section viewer for a ModEM model.
# Author: @pankajkmishra
# This script loads a model, builds 2D YZ cross-sections (slicing along the X axis),
# and provides controls for browsing/export.
# Use it for fast inspection of model structure along the X direction.
#
# The user chooses the visualisation coordinate system at the top of the script:
#   coordinate_system = "EPSG:3067"   — ETRS-TM35FIN (default, metres)
#   coordinate_system = "EPSG:4326"   — WGS 84 lat/lon (degrees)
#   coordinate_system = "model"       — raw model XY (metres, no reprojection)
#
# When using a projected or geographic CRS the script needs a ModEM data file
# (data_file) for georeferencing (origin lat/lon + station coordinates).

using Pkg

using GLMakie
using Statistics
using Dates
using Proj

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "Data.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))

# ---------- Model & data files ----------
model_file = joinpath(@__DIR__, ".", "I_NLCG_140.rho")
data_file  = joinpath(@__DIR__, ".", "I_NLCG_140.dat")   # needed for EPSG:3067 / EPSG:4326

# ---------- Coordinate system ----------
#   "EPSG:3067"  — ETRS-TM35FIN (Easting / Northing in metres)  ← default
#   "EPSG:4326"  — WGS 84 geographic (Longitude / Latitude in degrees)
#   "model"      — raw model XY (no conversion; data_file is not needed)
coordinate_system = "EPSG:3067"

# ---------- Visualisation controls ----------
log10_scale       = true
colormap          = :Spectral
with_padding      = false          # false → start with core view (default)
max_depth         = 50000          # Maximum depth (m) to visualise (default 50 km)
resistivity_range = (0.0, 4.0)
show_grid         = false
grid_color        = :black
grid_linewidth    = 0.5
grid_alpha        = 0.1
pad_tol           = 0.5

# ---------- Export settings ----------
export_dpi     = 3                # px_per_unit  (3 → ≈ 3300×2700 px at default size)
export_figsize = (1100, 900)      # export figure size (points)

# ---------- Coordinate transformation helpers ----------

function _local_tm_to_wgs84_transform(lat0::Real, lon0::Real)
    src = "+proj=tmerc +lat_0=$(float(lat0)) +lon_0=$(float(lon0)) +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs"
    return Proj.Transformation(src, "EPSG:4326"; always_xy = true)
end

function _convert_station_xy_to_latlon(d, trans)
    ns = size(d.loc, 1)
    lon_pred = zeros(ns)
    lat_pred = zeros(ns)
    for i in 1:ns
        ll = trans((Float64(d.y[i]) + 500000.0, Float64(d.x[i])))
        lon_pred[i] = Float64(ll[1])
        lat_pred[i] = Float64(ll[2])
    end
    return lat_pred, lon_pred
end

function model_xy_to_latlon_centers(M, d)
    lat_vals = d.loc[:, 1]
    lon_vals = d.loc[:, 2]

    lat0 = (length(d.origin) >= 1 && isfinite(d.origin[1])) ? Float64(d.origin[1]) : mean(lat_vals)
    lon0 = (length(d.origin) >= 2 && isfinite(d.origin[2])) ? Float64(d.origin[2]) : mean(lon_vals)

    trans = _local_tm_to_wgs84_transform(lat0, lon0)

    lon_centers = Float64[]
    for y in M.cy
        ll = trans((Float64(y) + 500000.0, 0.0))
        push!(lon_centers, Float64(ll[1]))
    end

    lat_centers = Float64[]
    for x in M.cx
        ll = trans((500000.0, Float64(x)))
        push!(lat_centers, Float64(ll[2]))
    end

    lat_pred, lon_pred = _convert_station_xy_to_latlon(d, trans)
    shiftlat = mean(lat_pred .- lat_vals)
    shiftlon = mean(lon_pred .- lon_vals)

    lat_centers .-= shiftlat
    lon_centers .-= shiftlon

    return lat_centers, lon_centers, lat0, lon0, shiftlat, shiftlon
end

function _latlon_to_epsg3067(lat_centers, lon_centers)
    trans = Proj.Transformation("EPSG:4326", "EPSG:3067"; always_xy = true)
    ref_lon = mean(lon_centers)
    ref_lat = mean(lat_centers)

    northing = Float64[]
    for lat in lat_centers
        p = trans((ref_lon, lat))
        push!(northing, Float64(p[2]))
    end

    easting = Float64[]
    for lon in lon_centers
        p = trans((lon, ref_lat))
        push!(easting, Float64(p[1]))
    end

    return northing, easting
end

# ---------- Build model in chosen CRS ----------
# For YZ cross-sections:
#   horizontal axis = Y (cols) → Easting  / Longitude / model-Y
#   slider axis     = X (rows) → Northing / Latitude  / model-X

function _build_yz_model_in_crs(M, d, crs::AbstractString)
    crs_up = uppercase(strip(crs))

    if crs_up == "MODEL"
        M_out = M
        horiz_label = "Y (m)"
        slider_label = "X"
        slider_unit = "m"
        println("  Coordinate system: raw model XY (metres)")
        return M_out, horiz_label, slider_label, slider_unit
    end

    lat_centers, lon_centers, lat0, lon0, shiftlat, shiftlon = model_xy_to_latlon_centers(M, d)

    println("  Georeferencing from data file:")
    println("    Origin (lat, lon): ($(round(lat0, digits=6)), $(round(lon0, digits=6)))")
    println("    Data alignment shift: Δlat=$(round(shiftlat, digits=8)), Δlon=$(round(shiftlon, digits=8))")

    if crs_up == "EPSG:4326"
        M_out = (A = M.A, cx = lat_centers, cy = lon_centers, cz = M.cz)
        horiz_label = "Longitude (°)"
        slider_label = "Latitude"
        slider_unit = "°"
        println("  Coordinate system: WGS 84 (EPSG:4326)")
        return M_out, horiz_label, slider_label, slider_unit

    elseif crs_up == "EPSG:3067"
        northing, easting = _latlon_to_epsg3067(lat_centers, lon_centers)
        M_out = (A = M.A, cx = northing, cy = easting, cz = M.cz)
        horiz_label = "Easting (m)"
        slider_label = "Northing"
        slider_unit = "m"
        println("  Coordinate system: ETRS-TM35FIN (EPSG:3067)")
        println("    Northing range: $(minimum(northing)) to $(maximum(northing))")
        println("    Easting range:  $(minimum(easting)) to $(maximum(easting))")
        return M_out, horiz_label, slider_label, slider_unit

    else
        error("Unsupported coordinate_system: \"$crs\".  Use \"EPSG:3067\", \"EPSG:4326\", or \"model\".")
    end
end

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

# ---------- Main ----------

function main()
    if !isfile(model_file)
        println("="^60)
        println("ERROR: Model file not found!")
        println("Please edit this script and set 'model_file' to your ModEM model path.")
        println("Current path: $model_file")
        println("="^60)
        return nothing, nothing
    end

    println("Loading ModEM model: $model_file")
    M = load_model_modem(model_file)

    println("Model loaded successfully!")
    println("  Grid size: $(M.nx) × $(M.ny) × $(M.nz)")
    println("  X range: $(minimum(M.x)) to $(maximum(M.x)) m")
    println("  Y range: $(minimum(M.y)) to $(maximum(M.y)) m")
    println("  Z range: $(minimum(M.z)) to $(maximum(M.z)) m")
    println("  Padding cells: $(M.npad)")

    # Build model in the chosen coordinate system
    crs_up = uppercase(strip(coordinate_system))
    need_data = crs_up != "MODEL"

    d = nothing
    if need_data
        if !isfile(data_file)
            println("="^60)
            println("ERROR: Data file not found!")
            println("  A ModEM data file is required for coordinate_system = \"$coordinate_system\".")
            println("  Set data_file at the top of this script, or use coordinate_system = \"model\".")
            println("  Current path: $data_file")
            println("="^60)
            return nothing, nothing
        end
        println("\nLoading ModEM data for georeferencing: $data_file")
        d = load_data_modem(data_file)
    end

    println("\nCoordinate conversion:")
    M_crs, horiz_label, slider_label, slider_unit = _build_yz_model_in_crs(M, d, coordinate_system)

    println("\nCreating interactive YZ cross-section viewer...")

    model_name = splitext(basename(model_file))[1]
    if crs_up != "MODEL"
        model_name *= "_$(replace(lowercase(strip(coordinate_system)), ":" => ""))"
    end

    fig, parts = yz_slice_viewer(M_crs;
        model_name      = model_name,
        log10scale      = log10_scale,
        cmap            = colormap,
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
        export_figsize  = export_figsize
    )

    println("\nViewer ready!")
    println("  - Coordinate system: $coordinate_system")
    println("  - Horizontal axis: $horiz_label")
    println("  - Slider axis: $slider_label ($slider_unit)")
    println("  - Use the slider to navigate through $(slider_label) slices (YZ cross-sections)")
    println("  - Use Previous/Next buttons for step-by-step navigation")
    println("  - Click 'Export Figure' for high-resolution PNG export")

    screen = display(fig)

    println("\nClose the figure window to exit...")
    wait(screen)

    return fig, parts
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

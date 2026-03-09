# Example: interactive 3D volume-slice model viewer.
# Author: @pankajkmishra
# This script visualizes XY/XZ/YZ slices with controls for depth, padding, and export
# Use it to explore overall 3D model geometry and resistivity trends

using Pkg
#Pkg.activate(dirname(@__DIR__))
using GLMakie
using Statistics
using Dates
using Shapefile
using GeoInterface
using Proj

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "Data.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))

# =========================
# User controls (edit here)
# =========================
model_file = joinpath(@__DIR__, "I_NLCG_070.rho")
data_file = joinpath(@__DIR__, "I_NLCG_140.dat")
shapefile_path = joinpath(@__DIR__, "gis", "Tnew", "Tnew.shp")

log10_scale = true
colormap = :Spectral
resistivity_range = (1.0, 4.0)
max_depth = 50000.0
show_padding = false
pad_tolerance = 0.2

default_view_direction = Vec3f(-1.05f0, -0.80f0, 0.72f0)
default_view_scale = 1.12f0

overlay_z_fixed = 0.0
overlay_auto_reproject_to_wgs84 = true
overlay_point_color = :black
overlay_line_color = :black
overlay_point_size = 7
overlay_line_width = 1.5
target_crs = "EPSG:3067"
show_north_arrow = true
show_scale_bar = true
annotation_color = :black
annotation_line_width = 2.0

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

function _resolve_wgs84_to_target_xy_transform(target_crs::AbstractString)
    crs = uppercase(strip(target_crs))
    if crs == "EPSG:4326"
        return (lon, lat) -> (Float64(lon), Float64(lat))
    end

    trans = Proj.Transformation("EPSG:4326", strip(target_crs); always_xy = true)
    return (lon, lat) -> begin
        p = trans((Float64(lon), Float64(lat)))
        return Float64(p[1]), Float64(p[2])
    end
end

function model_xy_to_target_crs_centers(M, d, target_crs::AbstractString)
    lat_centers, lon_centers, lat0, lon0, shiftlat, shiftlon = model_xy_to_latlon_centers(M, d)
    lat_ref = mean(lat_centers)
    lon_ref = mean(lon_centers)

    src_local_tm = "+proj=tmerc +lat_0=$(float(lat0)) +lon_0=$(float(lon0)) +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs"
    local_tm_to_target = Proj.Transformation(src_local_tm, strip(target_crs); always_xy = true)

    # Standard GIS orientation for plotting:
    # X axis = Easting (from model local y)
    # Y axis = Northing (from model local x)
    x_target = Float64[]
    for y_local in M.cy
        p = local_tm_to_target((Float64(y_local) + 500000.0, 0.0))
        push!(x_target, Float64(p[1]))
    end

    y_target = Float64[]
    for x_local in M.cx
        p = local_tm_to_target((500000.0, Float64(x_local)))
        push!(y_target, Float64(p[2]))
    end

    station_tx = Float64[]
    station_ty = Float64[]
    for i in eachindex(d.x)
        p = local_tm_to_target((Float64(d.y[i]) + 500000.0, Float64(d.x[i])))
        push!(station_tx, Float64(p[1]))
        push!(station_ty, Float64(p[2]))
    end

    epsv = eps(Float64)
    span_x_model = max(maximum(x_target) - minimum(x_target), epsv)   # Easting axis
    span_y_model = max(maximum(y_target) - minimum(y_target), epsv)   # Northing axis
    span_e_sta = max(maximum(station_tx) - minimum(station_tx), epsv) # Easting
    span_n_sta = max(maximum(station_ty) - minimum(station_ty), epsv) # Northing

    mismatch_dim_consistent = abs(log(span_x_model / span_e_sta)) + abs(log(span_y_model / span_n_sta))

    return x_target, y_target, lat0, lon0, shiftlat, shiftlon, lat_ref, lon_ref, mismatch_dim_consistent
end

function _nice_scale_length(target::Float64)
    target <= 0 && return 1.0
    expo = floor(log10(target))
    base = 10.0^expo
    frac = target / base
    nice_frac = if frac < 1.5
        1.0
    elseif frac < 3.5
        2.0
    elseif frac < 7.5
        5.0
    else
        10.0
    end
    return nice_frac * base
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

"""
Overlay shapefile geometries on a 3D plot at fixed z height.
shapefile_path: path to .shp file
z_fixed: height at which to plot shapes (default 0)
point_color, line_color, point_size, line_width: style controls
auto_reproject_to_wgs84: reproject projected CRS to WGS84
"""

# --- CRS detection and coordinate transform logic from 07 ---
function _prj_path_from_shp(shp::AbstractString)
    root, _ = splitext(shp)
    return root * ".prj"
end

function _detect_crs_type(shp::AbstractString)
    prj_path = _prj_path_from_shp(shp)
    if !isfile(prj_path)
        return :unknown, ""
    end
    wkt = read(prj_path, String)
    wktu = uppercase(wkt)
    if occursin("PROJCS", wktu) || occursin("PROJCRS", wktu)
        return :projected, wkt
    elseif occursin("GEOGCS", wktu) || occursin("GEOGRAPHICCRS", wktu)
        return :geographic, wkt
    else
        return :unknown, wkt
    end
end

function _make_coord_transform(shp_path::AbstractString, crs_type::Symbol; auto_reproject_to_wgs84::Bool = true)
    if !(auto_reproject_to_wgs84 && crs_type == :projected)
        return (x, y) -> (x, y), false, "none"
    end
    prj_path = _prj_path_from_shp(shp_path)
    if !isfile(prj_path)
        return (x, y) -> (x, y), false, "missing .prj"
    end
    wkt = read(prj_path, String)
    try
        trans = Proj.Transformation(wkt, "EPSG:4326"; always_xy = true)
        f = (x, y) -> begin
            ll = trans((x, y))
            return Float64(ll[1]), Float64(ll[2])
        end
        return f, true, "WKT -> EPSG:4326"
    catch
        return (x, y) -> (x, y), false, "transformation failed"
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
    if !isfile(shapefile_path)
        @warn "Shapefile not found: $shapefile_path"
        return 0
    end
    table = Shapefile.Table(shapefile_path)
    rows = collect(table)
    crs_type, _ = _detect_crs_type(shapefile_path)
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
    overlay_transform = (x, y) -> (x, y),
    north_axis::Symbol = :y
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
    plot_shapefile_on_3d!(ax.scene, shapefile_path;
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
        plot_shapefile_on_3d!(ax.scene, shapefile_path;
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
    println("Loading ModEM data for georeference from: $data_file")
    d = load_data_modem(data_file)

    x_target, y_target, lat0, lon0, shiftlat, shiftlon, lat_ref, lon_ref, mismatch_dim_consistent =
        model_xy_to_target_crs_centers(M, d, target_crs)
    wgs84_to_target = _resolve_wgs84_to_target_xy_transform(target_crs)
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
    println("  Target CRS: $target_crs")
    println("  Axis convention: X=Easting, Y=Northing (GIS-standard)")
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
        withPadding = show_padding,
        max_depth = max_depth,
        pad_tol = pad_tolerance,
        resistivity_range = resistivity_range,
        overlay_transform = overlay_transform,
        north_axis = north_axis
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

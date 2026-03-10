# Example: interactive depth-slice (XY) viewer for a ModEM model with shapefile overlay.
# Author: @pankajkmishra
# This script loads a model, builds 2D map slices by depth, and provides controls
# for browsing/export.  Optionally overlays one or more shapefiles on the map view.
#
# The user chooses the visualisation coordinate system at the top of the script:
#   coordinate_system = "EPSG:3067"   — ETRS-TM35FIN (default, metres)
#   coordinate_system = "EPSG:4326"   — WGS 84 lat/lon (degrees)
#   coordinate_system = "model"       — raw model XY (metres, no reprojection)
#
# When using a projected or geographic CRS the script needs a ModEM data file
# (data_file) for georeferencing (origin lat/lon + station coordinates).
# Shapefiles are automatically reprojected to the chosen CRS using their .prj
# sidecar files (via Proj.jl).

using Pkg

using GLMakie
using Statistics
using Dates
using Shapefile
using GeoInterface
using Proj

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "Data.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))

# ---------- Model & data files ----------
model_file = joinpath(@__DIR__, "Cascadia", "cascad_half_inverse.ws")
data_file  = joinpath(@__DIR__, "Cascadia", "cascad_errfl5.dat")   # needed for EPSG:3067 / EPSG:4326

# ---------- Coordinate system ----------
# Choose how the map axes are labelled and how model centres are converted:
#   "EPSG:3067"  — ETRS-TM35FIN (Easting / Northing in metres)  ← default
#   "EPSG:4326"  — WGS 84 geographic (Longitude / Latitude in degrees)
#   "model"      — raw model XY (no conversion; data_file is not needed)
coordinate_system = "EPSG:4326"

# ---------- Visualisation controls ----------
log10_scale       = true
colormap          = :Spectral
with_padding      = false          # false → start with core/extent view (default)
max_depth         = nothing
resistivity_range = (0.0, 4.0)
show_grid         = true
grid_color        = :black
grid_linewidth    = 0.5
grid_alpha        = 0.1
pad_tol           = 0.5

# ---------- Custom core extent (optional) ----------
# By default the "core model" is detected automatically from the grid spacing
# (cells whose width is within `pad_tol` of the median width).
#
# To define the core extent manually, set `custom_extent` to a NamedTuple with
# the bounding box in **geographic coordinates (lat/lon, WGS 84)**.
# The script converts it automatically to the chosen `coordinate_system`.
# Set to `nothing` to use the automatic detection.



custom_extent = nothing

# custom_extent = (min_lat = 63.0, max_lat = 65.0, min_lon = 25.0, max_lon = 28.0)

# ---------- Shapefile overlay configuration ----------
# Define shapefiles to overlay on the depth-slice map view.
# Each entry is a NamedTuple with the following fields:
#
#   path        – absolute or relative path to the .shp file
#   enabled     – true / false — set to false to skip without deleting the entry
#   color       – any Makie colour  (e.g. :red, :black, "#00FF00", RGBf(0,0,1))
#   alpha       – transparency  0.0 (invisible) … 1.0 (opaque)
#   point_size  – marker size for point geometries
#   line_width  – line width for line / polygon geometries
#
# Shapefiles are automatically reprojected from their native CRS (detected
# from the .prj sidecar file) into the chosen `coordinate_system`.
#
# Leave the list empty if you don't want any shapefile overlay:
#   shapefiles = []

shapefiles = [
   # (path = joinpath(@__DIR__, "gis", "2023MeasPts", "2023MeasPts.shp"), enabled = true,  color = :black, alpha = 0.9, point_size = 7,  line_width = 1.5),
    (path = joinpath(@__DIR__, "gis", "Tnew", "Tnew.shp"),              enabled = true,  color = :black,   alpha = 0.8, point_size = 5,  line_width = 2.0),
]

# ---------- Export settings ----------
export_dpi     = 3                # px_per_unit  (3 → ≈ 3300×2700 px at default size)
export_figsize = (1100, 900)      # export figure size (points)

# ---------- Shapefile helper functions ----------

function _is_xy_tuple(c)
    c isa Tuple && length(c) >= 2 && c[1] isa Real && c[2] isa Real
end

function _is_xy_vector(c)
    c isa AbstractVector && length(c) >= 2 && c[1] isa Real && c[2] isa Real
end

# ---------- Shapefile CRS detection and reprojection (ported from 02c) ----------

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

function _make_coord_transform_to_crs(shp_path::AbstractString, target_crs::AbstractString)
    crs_up = uppercase(strip(target_crs))
    # For raw model coordinates we cannot reproject shapefiles
    if crs_up == "MODEL"
        return (x, y) -> (x, y), false, "model coordinates (no reprojection)"
    end

    prj_path = _prj_path_from_shp(shp_path)
    if !isfile(prj_path)
        @warn "No .prj file for $shp_path — coordinates used as-is"
        return (x, y) -> (x, y), false, "missing .prj"
    end

    wkt = read(prj_path, String)
    try
        trans = Proj.Transformation(wkt, crs_up; always_xy = true)
        f = (x, y) -> begin
            p = trans((x, y))
            return Float64(p[1]), Float64(p[2])
        end
        return f, true, "WKT -> $crs_up"
    catch e
        @warn "Shapefile reprojection failed for $shp_path: $e"
        return (x, y) -> (x, y), false, "transformation failed"
    end
end

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

function _load_shapefile_geometries(shp_path::AbstractString)
    if !isfile(shp_path)
        @warn "Shapefile not found, skipping: $shp_path"
        return Any[]
    end
    table = Shapefile.Table(shp_path)
    return Any[GeoInterface.geometry(r) for r in collect(table)]
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

function _prepare_shapefiles(defs, target_crs::AbstractString)
    result = []
    for s in defs
        if !s.enabled
            println("  Shapefile DISABLED: $(s.path)")
            continue
        end
        geoms = _load_shapefile_geometries(s.path)
        if isempty(geoms)
            println("  Shapefile empty or not found: $(s.path)")
            continue
        end
        # Build coordinate transform from the shapefile's native CRS to target_crs
        crs_type, _ = _detect_crs_type(s.path)
        coord_transform, transformed, transform_info = _make_coord_transform_to_crs(s.path, target_crs)
        status = transformed ? "reprojected ($transform_info)" : "as-is ($transform_info)"
        println("  Shapefile loaded: $(s.path)  ($(length(geoms)) features, CRS: $crs_type, $status)")
        push!(result, (
            geoms           = geoms,
            color           = s.color,
            alpha           = s.alpha,
            point_size      = s.point_size,
            line_width      = s.line_width,
            path            = s.path,
            coord_transform = coord_transform
        ))
    end
    return result
end

# ---------- Custom / automatic core index selection ----------

function _indices_in_range(centers::AbstractVector, lo::Real, hi::Real)
    # Return the index range of cell centres that fall within [lo, hi].
    first_idx = findfirst(c -> c >= lo, centers)
    last_idx  = findlast(c -> c <= hi, centers)
    if isnothing(first_idx) || isnothing(last_idx) || first_idx > last_idx
        @warn "Custom extent ($lo .. $hi) does not overlap the model centres; falling back to full range."
        return 1:length(centers)
    end
    return first_idx:last_idx
end

function _core_indices_or_extent(x_all, y_all; pad_tol::Real = 0.2, extent = nothing)
    # x_all = row centres (plotted on y-axis)
    # y_all = col centres (plotted on x-axis)
    # extent (if not nothing) must already be in plot-CRS coords:
    #   (min_x, max_x) → column axis (y_all),  (min_y, max_y) → row axis (x_all)
    if isnothing(extent)
        ix = core_indices(x_all; tol = pad_tol)
        iy = core_indices(y_all; tol = pad_tol)
    else
        iy = _indices_in_range(y_all, extent.min_x, extent.max_x)
        ix = _indices_in_range(x_all, extent.min_y, extent.max_y)
    end
    return ix, iy
end

function _convert_extent_to_crs(extent, crs::AbstractString)
    # Convert a lat/lon extent NamedTuple to the chosen coordinate system.
    # Input:  (min_lat, max_lat, min_lon, max_lon)  in WGS 84
    # Output: (min_x, max_x, min_y, max_y) in the plot CRS, where
    #         x = column axis (Easting / Longitude) and y = row axis (Northing / Latitude)
    isnothing(extent) && return nothing

    crs_up = uppercase(strip(crs))

    if crs_up == "MODEL"
        @warn "custom_extent in lat/lon cannot be used with coordinate_system=\"model\". Ignoring."
        return nothing

    elseif crs_up == "EPSG:4326"
        # Already lat/lon — just remap field names
        return (min_x = extent.min_lon, max_x = extent.max_lon,
                min_y = extent.min_lat, max_y = extent.max_lat)

    else
        # Project the four corners of the bounding box
        trans = Proj.Transformation("EPSG:4326", crs_up; always_xy = true)
        corners = [
            trans((extent.min_lon, extent.min_lat)),
            trans((extent.max_lon, extent.min_lat)),
            trans((extent.min_lon, extent.max_lat)),
            trans((extent.max_lon, extent.max_lat)),
        ]
        xs = Float64[c[1] for c in corners]
        ys = Float64[c[2] for c in corners]
        return (min_x = minimum(xs), max_x = maximum(xs),
                min_y = minimum(ys), max_y = maximum(ys))
    end
end

# ---------- Coordinate transformation helpers (from 02a) ----------

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

function _wrap180(lon::Real)
    return mod(float(lon) + 180.0, 360.0) - 180.0
end

function _minimal_lon_span_deg(lons::AbstractVector{<:Real})
    vals = sort(mod.(Float64.(lons), 360.0))
    n = length(vals)
    n <= 1 && return 0.0
    max_gap = -Inf
    for i in 1:(n - 1)
        max_gap = max(max_gap, vals[i + 1] - vals[i])
    end
    max_gap = max(max_gap, (vals[1] + 360.0) - vals[end])
    return clamp(360.0 - max_gap, 0.0, 360.0)
end

function _haversine_km(lat1::Real, lon1::Real, lat2::Real, lon2::Real)
    R = 6371.0088
    φ1 = deg2rad(float(lat1))
    φ2 = deg2rad(float(lat2))
    Δφ = deg2rad(float(lat2 - lat1))
    Δλ = deg2rad(_wrap180(float(lon2 - lon1)))
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    c = 2 * atan(sqrt(max(a, 0.0)), sqrt(max(1.0 - a, 0.0)))
    return R * c
end

function _distance_based_aspect(lat_vals::AbstractVector{<:Real}, lon_vals::AbstractVector{<:Real})
    lat_min = minimum(lat_vals)
    lat_max = maximum(lat_vals)
    lat_ref = clamp(mean(lat_vals), -89.999, 89.999)
    lon_span = _minimal_lon_span_deg(lon_vals)
    width_km  = _haversine_km(lat_ref, 0.0, lat_ref, lon_span)
    height_km = _haversine_km(lat_min, 0.0, lat_max, 0.0)
    return width_km / max(height_km, eps(Float64)), lat_ref
end

function _latlon_to_epsg3067(lat_centers, lon_centers)
    # Convert independent lat/lon axis arrays to ETRS-TM35FIN (EPSG:3067).
    # Each axis is converted independently using a reference value for the other.
    trans = Proj.Transformation("EPSG:4326", "EPSG:3067"; always_xy = true)
    ref_lon = mean(lon_centers)
    ref_lat = mean(lat_centers)

    # Convert each lat center at reference longitude → extract northing
    northing = Float64[]
    for lat in lat_centers
        p = trans((ref_lon, lat))
        push!(northing, Float64(p[2]))
    end

    # Convert each lon center at reference latitude → extract easting
    easting = Float64[]
    for lon in lon_centers
        p = trans((lon, ref_lat))
        push!(easting, Float64(p[1]))
    end

    return northing, easting   # northing ↔ X (row), easting ↔ Y (col)
end

# ---------- Build model in chosen CRS ----------

function _build_model_in_crs(M, d, crs::AbstractString, crs_extent = nothing)
    crs_up = uppercase(strip(crs))

    if crs_up == "MODEL"
        # Raw model coordinates — no transformation needed
        M_out = M
        xlabel = "Y (m)"
        ylabel = "X (m)"
        aspect = nothing
        println("  Coordinate system: raw model XY (metres)")
        return M_out, xlabel, ylabel, aspect
    end

    # Geographic transform via data file
    lat_centers, lon_centers, lat0, lon0, shiftlat, shiftlon = model_xy_to_latlon_centers(M, d)

    println("  Georeferencing from data file:")
    println("    Origin (lat, lon): ($(round(lat0, digits=6)), $(round(lon0, digits=6)))")
    println("    Data alignment shift: Δlat=$(round(shiftlat, digits=8)), Δlon=$(round(shiftlon, digits=8))")

    if crs_up == "EPSG:4326"
        # WGS 84 lat/lon
        M_out = (A = M.A, cx = lat_centers, cy = lon_centers, cz = M.cz)
        xlabel = "Longitude (°)"
        ylabel = "Latitude (°)"

        ix_core, iy_core = _core_indices_or_extent(lat_centers, lon_centers; pad_tol = pad_tol, extent = crs_extent)
        lat_for_asp = with_padding ? lat_centers : lat_centers[ix_core]
        lon_for_asp = with_padding ? lon_centers : lon_centers[iy_core]
        map_asp, lat_ref = _distance_based_aspect(lat_for_asp, lon_for_asp)
        aspect = map_asp

        println("  Coordinate system: WGS 84 (EPSG:4326)")
        println("    Latitude range:  $(minimum(lat_centers)) to $(maximum(lat_centers))")
        println("    Longitude range: $(minimum(lon_centers)) to $(maximum(lon_centers))")
        println("    Distance-based aspect (W/H): $(round(map_asp, digits=4)) at lat̄=$(round(lat_ref, digits=3))°")
        return M_out, xlabel, ylabel, aspect

    elseif crs_up == "EPSG:3067"
        # ETRS-TM35FIN
        northing, easting = _latlon_to_epsg3067(lat_centers, lon_centers)
        M_out = (A = M.A, cx = northing, cy = easting, cz = M.cz)
        xlabel = "Easting (m)"
        ylabel = "Northing (m)"
        aspect = nothing   # DataAspect — metres in both axes

        println("  Coordinate system: ETRS-TM35FIN (EPSG:3067)")
        println("    Easting range:  $(minimum(easting)) to $(maximum(easting))")
        println("    Northing range: $(minimum(northing)) to $(maximum(northing))")
        return M_out, xlabel, ylabel, aspect

    else
        error("Unsupported coordinate_system: \"$crs\".  Use \"EPSG:3067\", \"EPSG:4326\", or \"model\".")
    end
end

# ---------- WKT strings for .prj files ----------

function _default_wgs84_wkt()
    return "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]]"
end

function _default_tm35fin_wkt()
    return "PROJCS[\"ETRS89 / TM35FIN(E,N)\",GEOGCS[\"ETRS89\",DATUM[\"European_Terrestrial_Reference_System_1989\",SPHEROID[\"GRS 1980\",6378137,298.257222101,AUTHORITY[\"EPSG\",\"7019\"]],AUTHORITY[\"EPSG\",\"6258\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4258\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",27],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"3067\"]]"
end

function _resolve_prj_wkt_for_crs(crs::AbstractString)
    c = uppercase(strip(crs))
    c == "EPSG:4326" && return _default_wgs84_wkt()
    c == "EPSG:3067" && return _default_tm35fin_wkt()
    return ""
end

# ---------- GIS shapefile export for all depth slices ----------

function _export_all_depth_slices_gis(;
    output_dir::AbstractString,
    model_name::AbstractString,
    x_edges::AbstractVector,   # row edges (Northing / Lat / X)  → shapefile Y
    y_edges::AbstractVector,   # col edges (Easting  / Lon / Y)  → shapefile X
    R::AbstractArray,          # nx × ny × nz  (already log10 if log10scale)
    z_edges::AbstractVector,   # depth edges
    crs::AbstractString,
    log10scale::Bool,
    qml_source_path::AbstractString,
    lyrx_source_path::AbstractString = "",
    vis_range::Tuple{<:Real,<:Real} = (0.0, 4.0))

    mkpath(output_dir)

    wkt = _resolve_prj_wkt_for_crs(crs)
    nx, ny, nz = size(R)

    vis_lo = Float64(min(vis_range[1], vis_range[2]))
    vis_hi = Float64(max(vis_range[1], vis_range[2]))

    layer_depths = cumsum(diff(z_edges))

    println("\nExporting $nz depth slices as shapefiles to:")
    println("  $output_dir")

    for k in 1:nz
        depth_top = k == 1 ? 0.0 : layer_depths[k - 1]
        depth_bot = layer_depths[k]

        # Build polygon features for this slice
        polygons = Shapefile.Polygon[]
        row_vals     = Int32[]
        col_vals     = Int32[]
        rho_vals     = Float64[]
        rho_ohmm_vals = Float64[]
        rho_vis_vals = Float64[]
        layer_vals   = Int32[]
        ztop_vals    = Float64[]
        zbot_vals    = Float64[]

        for i in 1:nx
            for j in 1:ny
                val = Float64(R[i, j, k])
                isfinite(val) || continue

                # Shapefile X = col axis (Easting/Lon/Y), Y = row axis (Northing/Lat/X)
                x1 = Float64(y_edges[j])
                x2 = Float64(y_edges[j + 1])
                y1 = Float64(x_edges[i])
                y2 = Float64(x_edges[i + 1])

                pts = Shapefile.Point[
                    Shapefile.Point(x1, y1),
                    Shapefile.Point(x2, y1),
                    Shapefile.Point(x2, y2),
                    Shapefile.Point(x1, y2),
                    Shapefile.Point(x1, y1)   # close polygon
                ]

                rect = Shapefile.Rect(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
                poly = Shapefile.Polygon(rect, Int32[0], pts)
                push!(polygons, poly)

                v_log10 = log10scale ? val : log10(max(val, eps(Float64)))
                v_ohmm  = log10scale ? 10.0 ^ val : val
                v_vis   = clamp(v_log10, vis_lo, vis_hi)

                push!(row_vals, Int32(i))
                push!(col_vals, Int32(j))
                push!(rho_vals, round(v_log10, digits = 6))
                push!(rho_ohmm_vals, round(v_ohmm, digits = 6))
                push!(rho_vis_vals, round(v_vis, digits = 6))
                push!(layer_vals, Int32(k))
                push!(ztop_vals, round(depth_top, digits = 2))
                push!(zbot_vals, round(depth_bot, digits = 2))
            end
        end

        feats = (
            row      = row_vals,
            col      = col_vals,
            rho      = rho_vals,
            rho_ohmm = rho_ohmm_vals,
            rho_vis  = rho_vis_vals,
            layer    = layer_vals,
            ztop_m   = ztop_vals,
            zbot_m   = zbot_vals
        )

        depth_str = "$(round(Int, depth_top))m_$(round(Int, depth_bot))m"
        filename = "$(model_name)_depth_$(depth_str).shp"
        outpath = joinpath(output_dir, filename)

        writer = Shapefile.Writer(polygons, feats)
        Shapefile.write(outpath, writer; force = true)

        # Write .prj and .qpj sidecar files
        if !isempty(wkt)
            base = splitext(outpath)[1]
            open(base * ".prj", "w") do io; write(io, wkt); end
            open(base * ".qpj", "w") do io; write(io, wkt); end
        end

        println("  Layer $k/$nz: $filename  ($(length(polygons)) cells, depth $(round(depth_top, digits=1))–$(round(depth_bot, digits=1)) m)")
    end

    # Copy QML style file (QGIS)
    if isfile(qml_source_path)
        dst = joinpath(output_dir, basename(qml_source_path))
        cp(qml_source_path, dst; force = true)
        println("  Style file copied: $(basename(qml_source_path))  (QGIS)")
    else
        @warn "QML style file not found: $qml_source_path"
    end

    # Copy LYRX style file (ArcGIS Pro)
    if !isempty(lyrx_source_path) && isfile(lyrx_source_path)
        dst = joinpath(output_dir, basename(lyrx_source_path))
        cp(lyrx_source_path, dst; force = true)
        println("  Style file copied: $(basename(lyrx_source_path))  (ArcGIS Pro)")
    elseif !isempty(lyrx_source_path)
        @warn "LYRX style file not found: $lyrx_source_path"
    end

    println("GIS export complete: $nz shapefiles written.")
    return nz
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
    loaded_shapefiles = [],
    export_dpi::Int = 3,
    export_figsize = (1100, 900),
    custom_extent = nothing,
    export_crs::String = "model",
    qml_source_path::String = "",
    lyrx_source_path::String = ""
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
    shp_info_str = n_shp > 0 ? "Shapefiles loaded: $n_shp" : "No shapefiles"
    view_info = Observable(withPadding ? "View: Full Model | $shp_info_str" : "View: Core Model | $shp_info_str")
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
            view_info[] = "View: Full Model | $shp_info_str"
        else
            current_x_edges[] = x_edges_core
            current_y_edges[] = y_edges_core
            current_R_ref[] = R_core
            btn_label[] = "Show Full Model"
            view_info[] = "View: Core Model | $shp_info_str"
        end

        layer_idx = current_layer[]
        slice_data[] = current_R_ref[][:, :, layer_idx]'

        empty!(ax)
        hm = heatmap!(ax, current_y_edges[], current_x_edges[], slice_data,
                      colormap = current_colormap,
                      colorrange = (cmin, cmax))

        draw_grid!(ax, current_x_edges[], current_y_edges[])

        # Redraw shapefiles after empty!()
        _draw_all_shapefiles!(ax, loaded_shapefiles)

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

        gis_dir = joinpath(pwd(), "$(model_name)-GIS")

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
    # Convert user-defined lat/lon extent to the plot coordinate system
    crs_extent = _convert_extent_to_crs(custom_extent, coordinate_system)
    if !isnothing(custom_extent)
        if isnothing(crs_extent)
            println("  Custom extent: ignored (not supported for this CRS)")
        else
            println("  Custom extent (lat/lon → $(coordinate_system)):")
            println("    x-axis: $(crs_extent.min_x) .. $(crs_extent.max_x)")
            println("    y-axis: $(crs_extent.min_y) .. $(crs_extent.max_y)")
        end
    end
    M_crs, xlabel, ylabel, aspect = _build_model_in_crs(M, d, coordinate_system, crs_extent)

    # Load enabled shapefiles
    println("\nShapefile overlay:")
    if isempty(shapefiles)
        println("  No shapefiles defined.")
    end
    loaded_shapefiles = _prepare_shapefiles(shapefiles, coordinate_system)

    println("\nCreating interactive depth-slice viewer...")

    model_name = splitext(basename(model_file))[1]
    if crs_up != "MODEL"
        model_name *= "_$(replace(lowercase(strip(coordinate_system)), ":" => ""))"
    end

    # Resolve paths to style files
    qml_path  = joinpath(dirname(@__DIR__), "src", "MTSlices.qml")
    lyrx_path = joinpath(dirname(@__DIR__), "src", "MTSlices.lyrx")

    fig, parts = depth_slice_viewer(M_crs;
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
        axis_xlabel     = xlabel,
        axis_ylabel     = ylabel,
        axis_aspect     = aspect,
        loaded_shapefiles = loaded_shapefiles,
        export_dpi      = export_dpi,
        export_figsize  = export_figsize,
        custom_extent   = crs_extent,
        export_crs      = coordinate_system,
        qml_source_path = qml_path,
        lyrx_source_path = lyrx_path
    )

    println("\nViewer ready!")
    println("  - Coordinate system: $coordinate_system")
    println("  - Use the slider to navigate through depth layers")
    println("  - Use Previous/Next buttons for step-by-step navigation")
    println("  - Click 'Export Figure' for high-resolution PNG export")
    if !isempty(loaded_shapefiles)
        println("  - $(length(loaded_shapefiles)) shapefile(s) overlaid on the depth slice")
    end

    screen = display(fig)

    println("\nClose the figure window to exit...")
    wait(screen)

    return fig, parts
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

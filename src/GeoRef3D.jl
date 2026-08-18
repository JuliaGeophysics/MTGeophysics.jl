function _looks_like_coordinate_system_arg(arg::AbstractString)
    arg_up = uppercase(strip(arg))
    return arg_up == "MODEL" || startswith(arg_up, "EPSG:")
end

function _local_tm_proj_string(lat0::Real, lon0::Real)
    return "+proj=tmerc +lat_0=$(float(lat0)) +lon_0=$(float(lon0)) +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs"
end

_is_project_local_crs(crs::AbstractString) = uppercase(strip(crs)) == "PROJECT"

function _local_tm_to_wgs84_transform(lat0::Real, lon0::Real)
    src = _local_tm_proj_string(lat0, lon0)
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

function _latlon_to_projected(lat_centers, lon_centers, crs::AbstractString)
    # Convert independent lat/lon axis arrays to any projected CRS.
    # Each axis is converted independently using a reference value for the other.
    trans = Proj.Transformation("EPSG:4326", uppercase(strip(crs)); always_xy = true)
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

function _build_model_in_crs(M, d, crs::AbstractString, crs_extent = nothing;
    pad_tol::Real = 0.5, with_padding::Bool = false)
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

    else
        # Any projected CRS (EPSG:3067, EPSG:32610, EPSG:26910, EPSG:3005, ...)
        northing, easting = _latlon_to_projected(lat_centers, lon_centers, crs_up)
        M_out = (A = M.A, cx = northing, cy = easting, cz = M.cz)
        xlabel = "Easting (m)"
        ylabel = "Northing (m)"
        aspect = nothing   # DataAspect — metres in both axes

        println("  Coordinate system: $crs_up (projected)")
        println("    Easting range:  $(minimum(easting)) to $(maximum(easting))")
        println("    Northing range: $(minimum(northing)) to $(maximum(northing))")
        return M_out, xlabel, ylabel, aspect
    end
end

# ---------- Build model in chosen CRS ----------
# For XZ cross-sections:
#   horizontal axis = X (rows) → Northing / Latitude / model-X
#   slider axis     = Y (cols) → Easting  / Longitude / model-Y

function _build_xz_model_in_crs(M, d, crs::AbstractString)
    crs_up = uppercase(strip(crs))

    if crs_up == "MODEL"
        M_out = M
        horiz_label = "X (m)"
        slider_label = "Y"
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
        horiz_label = "Latitude (°)"
        slider_label = "Longitude"
        slider_unit = "°"
        println("  Coordinate system: WGS 84 (EPSG:4326)")
        return M_out, horiz_label, slider_label, slider_unit

    else
        northing, easting = _latlon_to_projected(lat_centers, lon_centers, crs_up)
        M_out = (A = M.A, cx = northing, cy = easting, cz = M.cz)
        horiz_label = "Northing (m)"
        slider_label = "Easting"
        slider_unit = "m"
        println("  Coordinate system: $crs_up (projected)")
        println("    Northing range: $(minimum(northing)) to $(maximum(northing))")
        println("    Easting range:  $(minimum(easting)) to $(maximum(easting))")
        return M_out, horiz_label, slider_label, slider_unit
    end
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

    else
        northing, easting = _latlon_to_projected(lat_centers, lon_centers, crs_up)
        M_out = (A = M.A, cx = northing, cy = easting, cz = M.cz)
        horiz_label = "Easting (m)"
        slider_label = "Northing"
        slider_unit = "m"
        println("  Coordinate system: $crs_up (projected)")
        println("    Northing range: $(minimum(northing)) to $(maximum(northing))")
        println("    Easting range:  $(minimum(easting)) to $(maximum(easting))")
        return M_out, horiz_label, slider_label, slider_unit
    end
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

# ---------- station locations in the chosen plot CRS ----------

function _stations_in_crs(d, crs::AbstractString)
    isnothing(d) && return (Float64[], Float64[])
    crs_up = uppercase(strip(crs))
    if crs_up == "MODEL"
        return (collect(Float64.(d.y)), collect(Float64.(d.x)))
    elseif crs_up == "EPSG:4326"
        return (collect(Float64.(d.loc[:, 2])), collect(Float64.(d.loc[:, 1])))
    else
        northing, easting = _latlon_to_projected(collect(Float64.(d.loc[:, 1])), collect(Float64.(d.loc[:, 2])), crs_up)
        return (collect(Float64.(easting)), collect(Float64.(northing)))
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
    c == "MODEL" && return ""
    c == "EPSG:4326" && return _default_wgs84_wkt()
    c == "EPSG:3067" && return _default_tm35fin_wkt()
    try
        return Proj.proj_as_wkt(Proj.CRS(c), Proj.PJ_WKT1_ESRI)
    catch e
        @warn "Could not derive .prj WKT for $crs ($e); sidecar files will be omitted"
        return ""
    end
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



function _plain_tick_label(v::Real)
    av = abs(v)
    if av >= 1000 || isinteger(v)
        return string(round(Int, v))
    elseif av >= 1
        return string(round(v, digits = 3))
    else
        return string(round(v, digits = 5))
    end
end

_plain_tickformat(values) = [_plain_tick_label(v) for v in values]

function _resolve_wgs84_to_target_xy_transform(target_crs::AbstractString, lat0::Real, lon0::Real)
    crs = uppercase(strip(target_crs))
    if crs == "EPSG:4326"
        return (lon, lat) -> (Float64(lon), Float64(lat))
    elseif _is_project_local_crs(crs)
        trans = Proj.Transformation("EPSG:4326", _local_tm_proj_string(lat0, lon0); always_xy = true)
        return (lon, lat) -> begin
            p = trans((Float64(lon), Float64(lat)))
            return Float64(p[1]), Float64(p[2])
        end
    end

    trans = Proj.Transformation("EPSG:4326", strip(target_crs); always_xy = true)
    return (lon, lat) -> begin
        p = trans((Float64(lon), Float64(lat)))
        return Float64(p[1]), Float64(p[2])
    end
end

function _suggest_utm_crs(lat::Real, lon::Real)
    zone = clamp(floor(Int, (Float64(lon) + 180.0) / 6.0) + 1, 1, 60)
    prefix = lat >= 0 ? 326 : 327
    return "EPSG:$(prefix)$(lpad(string(zone), 2, '0'))"
end

_is_finland_extent(lat::Real, lon::Real) = 58.0 <= lat <= 71.5 && 18.0 <= lon <= 33.5

function _normalize_target_crs(requested_crs::AbstractString)
    crs = strip(requested_crs)
    crs_up = uppercase(crs)
    if _is_project_local_crs(crs_up)
        return "PROJECT"
    end
    startswith(crs_up, "EPSG:") && return crs_up
    error("Unsupported target_crs: \"$requested_crs\". Use \"PROJECT\" for the safe local survey coordinates, or pass an explicit EPSG code such as \"EPSG:3067\" or \"EPSG:32610\".")
end

function _warn_if_target_crs_looks_unexpected(d, target_crs::AbstractString)
    size(d.loc, 1) == 0 && return

    crs = uppercase(strip(target_crs))
    _is_project_local_crs(crs) && return
    lat_mean = mean(d.loc[:, 1])
    lon_mean = mean(d.loc[:, 2])
    in_finland = _is_finland_extent(lat_mean, lon_mean)

    if crs == "EPSG:3067" && !in_finland
        utm_hint = _suggest_utm_crs(lat_mean, lon_mean)
        @warn "target_crs=EPSG:3067 is Finland-specific and does not match the data-file coordinates. Use a projected CRS for the survey region instead." mean_lat=lat_mean mean_lon=lon_mean suggested_projected_crs=utm_hint
    elseif crs == "EPSG:4326"
        projected_hint = in_finland ? "EPSG:3067" : _suggest_utm_crs(lat_mean, lon_mean)
        @warn "3D viewing in EPSG:4326 mixes degree-based horizontal axes with metre-based depth, so the model can look distorted. Prefer a projected CRS." mean_lat=lat_mean mean_lon=lon_mean suggested_projected_crs=projected_hint
    end
end

function _axis_metadata_for_crs(target_crs::AbstractString)
    crs = uppercase(strip(target_crs))
    if crs == "EPSG:4326"
        return (x_name = "Longitude", y_name = "Latitude", unit = "deg")
    end
    return (x_name = "Easting", y_name = "Northing", unit = "m")
end

function model_xy_to_target_crs_centers(M, d, target_crs::AbstractString)
    lat_centers, lon_centers, lat0, lon0, shiftlat, shiftlon = model_xy_to_latlon_centers(M, d)
    lat_ref = mean(lat_centers)
    lon_ref = mean(lon_centers)

    x_target = Float64[]
    y_target = Float64[]
    station_tx = Float64[]
    station_ty = Float64[]

    if _is_project_local_crs(target_crs)
        # Project coordinates are the local TM metres implied by the ModEM
        # origin and station XY offsets: Easting from local y, Northing from local x.
        x_target = collect(Float64.(M.cy .+ 500000.0))
        y_target = collect(Float64.(M.cx))
        station_tx = collect(Float64.(d.y .+ 500000.0))
        station_ty = collect(Float64.(d.x))
    else
        src_local_tm = _local_tm_proj_string(lat0, lon0)
        local_tm_to_target = Proj.Transformation(src_local_tm, strip(target_crs); always_xy = true)

        # Standard GIS orientation for plotting:
        # X axis = Easting (from model local y)
        # Y axis = Northing (from model local x)
        for y_local in M.cy
            p = local_tm_to_target((Float64(y_local) + 500000.0, 0.0))
            push!(x_target, Float64(p[1]))
        end

        for x_local in M.cx
            p = local_tm_to_target((500000.0, Float64(x_local)))
            push!(y_target, Float64(p[2]))
        end

        for i in eachindex(d.x)
            p = local_tm_to_target((Float64(d.y[i]) + 500000.0, Float64(d.x[i])))
            push!(station_tx, Float64(p[1]))
            push!(station_ty, Float64(p[2]))
        end
    end

    epsv = eps(Float64)
    span_x_model = max(maximum(x_target) - minimum(x_target), epsv)   # Easting axis
    span_y_model = max(maximum(y_target) - minimum(y_target), epsv)   # Northing axis
    span_e_sta = max(maximum(station_tx) - minimum(station_tx), epsv) # Easting
    span_n_sta = max(maximum(station_ty) - minimum(station_ty), epsv) # Northing

    mismatch_dim_consistent = abs(log(span_x_model / span_e_sta)) + abs(log(span_y_model / span_n_sta))

    return x_target, y_target, lat0, lon0, shiftlat, shiftlon, lat_ref, lon_ref, mismatch_dim_consistent, station_tx, station_ty
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

function _make_coord_transform(shp_path::AbstractString, crs_type::Symbol; auto_reproject_to_wgs84::Bool = true)
    if !(auto_reproject_to_wgs84 && crs_type == :projected)
        return (x, y) -> (x, y), false, "none"
    end
    prj_path = prj_path_from_shp(shp_path)
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

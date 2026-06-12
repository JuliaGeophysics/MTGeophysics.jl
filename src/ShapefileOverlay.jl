"""
    prj_path_from_shp(shp)

Path of the `.prj` sidecar file next to a `.shp` file.
"""
prj_path_from_shp(shp::AbstractString) = first(splitext(shp)) * ".prj"

"""
    detect_shapefile_crs(shp)

Return `(crs_type, wkt)` where `crs_type` is `:projected`, `:geographic`, or
`:unknown`, based on the shapefile's `.prj` sidecar contents.
"""
function detect_shapefile_crs(shp::AbstractString)
    prj_path = prj_path_from_shp(shp)
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

"""
    shapefile_coord_transform(shp_path, target_crs)

Return `(f, transformed, info)` where `f(x, y)` maps shapefile coordinates into
`target_crs` using the `.prj` sidecar. Falls back to the identity transform
(with `transformed = false`) when `target_crs` is `"model"`, the `.prj` file is
missing, or the projection cannot be constructed.
"""
function shapefile_coord_transform(shp_path::AbstractString, target_crs::AbstractString)
    crs_up = uppercase(strip(target_crs))
    if crs_up == "MODEL"
        return (x, y) -> (x, y), false, "model coordinates (no reprojection)"
    end

    prj_path = prj_path_from_shp(shp_path)
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

"""
    load_shapefile_geometries(shp_path)

All geometries of a shapefile as a `Vector{Any}`; empty (with a warning) when
the file does not exist.
"""
function load_shapefile_geometries(shp_path::AbstractString)
    if !isfile(shp_path)
        @warn "Shapefile not found, skipping: $shp_path"
        return Any[]
    end
    table = Shapefile.Table(shp_path)
    return Any[GeoInterface.geometry(r) for r in collect(table)]
end

"""
    prepare_shapefiles(defs, target_crs)

Load and reproject a list of shapefile overlay definitions (NamedTuples with
fields `path`, `enabled`, `color`, `alpha`, `point_size`, `line_width`) into
ready-to-draw entries carrying the geometries and a coordinate transform into
`target_crs`. Disabled, empty, and missing entries are skipped with a message.
"""
function prepare_shapefiles(defs, target_crs::AbstractString)
    result = []
    for s in defs
        if !s.enabled
            println("  Shapefile DISABLED: $(s.path)")
            continue
        end
        geoms = load_shapefile_geometries(s.path)
        if isempty(geoms)
            println("  Shapefile empty or not found: $(s.path)")
            continue
        end
        crs_type, _ = detect_shapefile_crs(s.path)
        ct, transformed, transform_info = shapefile_coord_transform(s.path, target_crs)
        status = transformed ? "reprojected ($transform_info)" : "as-is ($transform_info)"
        println("  Shapefile loaded: $(s.path)  ($(length(geoms)) features, CRS: $crs_type, $status)")
        push!(result, (
            geoms           = geoms,
            color           = s.color,
            alpha           = s.alpha,
            point_size      = s.point_size,
            line_width      = s.line_width,
            path            = s.path,
            coord_transform = ct
        ))
    end
    return result
end

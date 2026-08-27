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

# ---------- interactive shapefile chooser ----------

const _SHAPEFILE_PICKER_PS = """
Add-Type -AssemblyName System.Windows.Forms
\$f = New-Object System.Windows.Forms.OpenFileDialog
\$f.Filter = 'Shapefile (*.shp)|*.shp|All files (*.*)|*.*'
if (\$f.ShowDialog() -eq 'OK') { Write-Output \$f.FileName }
"""

"""
    pick_shapefile()

Open the platform file chooser for a shapefile, for viewers that import overlays
while running.
In:  nothing.
Out: the chosen path, or "" if cancelled or no chooser is available.
"""
function pick_shapefile()
    try
        if Sys.isapple()
            return readchomp(`osascript -e $("POSIX path of (choose file with prompt \"Select a shapefile\" of type {\"shp\"})")`)
        elseif Sys.iswindows()
            return readchomp(`powershell -NoProfile -Command $_SHAPEFILE_PICKER_PS`)
        elseif !isnothing(Sys.which("zenity"))
            return readchomp(`zenity --file-selection --title=$("Select a shapefile") --file-filter=$("Shapefiles | *.shp")`)
        elseif !isnothing(Sys.which("kdialog"))
            return readchomp(`kdialog --getopenfilename $(pwd()) $("*.shp")`)
        elseif !isnothing(Sys.which("yad"))
            return readchomp(`yad --file --file-filter=$("*.shp")`)
        end
    catch
    end
    return ""
end

# ---------- clipping overlays to a plot extent ----------

"""
    _clip_segment_to_box(x1, y1, x2, y2, xmin, xmax, ymin, ymax)

Liang-Barsky clip of one segment against an axis-aligned box.
In:  segment endpoints and box bounds.
Out: clipped (x1, y1, x2, y2), or `nothing` if the segment misses the box.
"""
function _clip_segment_to_box(x1, y1, x2, y2, xmin, xmax, ymin, ymax)
    dx, dy = x2 - x1, y2 - y1
    t0, t1 = 0.0, 1.0
    for (pk, qk) in ((-dx, x1 - xmin), (dx, xmax - x1), (-dy, y1 - ymin), (dy, ymax - y1))
        if pk == 0
            qk < 0 && return nothing
        else
            r = qk / pk
            if pk < 0
                r > t1 && return nothing
                r > t0 && (t0 = r)
            else
                r < t0 && return nothing
                r < t1 && (t1 = r)
            end
        end
    end
    return (x1 + t0*dx, y1 + t0*dy, x1 + t1*dx, y1 + t1*dy)
end

"""
    _clip_polyline_to_box(xs, ys, lim)

Clip a polyline to a box, splitting it where it leaves and re-enters.
In:  polyline vertices and lim = (xmin, xmax, ymin, ymax).
Out: vector of (xs, ys) pieces inside the box; empty if none are.
"""
function _clip_polyline_to_box(xs, ys, lim)
    xmin, xmax, ymin, ymax = lim
    pieces = Tuple{Vector{Float64},Vector{Float64}}[]
    cx = Float64[]; cy = Float64[]
    flush!() = (length(cx) > 1 && push!(pieces, (copy(cx), copy(cy))); empty!(cx); empty!(cy))
    for k in 1:length(xs)-1
        seg = _clip_segment_to_box(xs[k], ys[k], xs[k+1], ys[k+1], xmin, xmax, ymin, ymax)
        if isnothing(seg)
            flush!()
            continue
        end
        ax1, ay1, ax2, ay2 = seg
        if isempty(cx) || cx[end] != ax1 || cy[end] != ay1
            flush!()
            push!(cx, ax1); push!(cy, ay1)
        end
        push!(cx, ax2); push!(cy, ay2)
    end
    flush!()
    return pieces
end

"""
    _inside_box(x, y, lim)

Whether a point lies inside lim = (xmin, xmax, ymin, ymax).
"""
_inside_box(x, y, lim) = lim[1] <= x <= lim[2] && lim[3] <= y <= lim[4]

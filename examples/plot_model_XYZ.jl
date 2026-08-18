using MTGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : ""
data_file  = length(ARGS) >= 2 ? ARGS[2] : ""
target_crs = length(ARGS) >= 3 ? ARGS[3] : "PROJECT"

shapefile_path = nothing

log10_scale       = true
colormap          = :Spectral
resistivity_range = (0.0, 4.0)
max_depth         = 5000.0
show_padding      = false
pad_tolerance     = 0.2
viewer_figsize    = (1800, 920)

default_view_direction = (-1.05, -0.80, 0.72)
default_view_scale     = 1.12

overlay_z_fixed                 = 0.0
overlay_auto_reproject_to_wgs84 = true
overlay_point_color             = :black
overlay_line_color              = :black
overlay_point_size              = 7
overlay_line_width              = 1.5
show_north_arrow                = true
show_scale_bar                  = true
annotation_color                = :black
annotation_line_width           = 2.0

PlotModelXYZ(model_file, data_file;
    crs                             = target_crs,
    shapefile_path                  = shapefile_path,
    log10_scale                     = log10_scale,
    colormap                        = colormap,
    resistivity_range               = resistivity_range,
    max_depth                       = max_depth,
    show_padding                    = show_padding,
    pad_tolerance                   = pad_tolerance,
    viewer_figsize                  = viewer_figsize,
    default_view_direction          = default_view_direction,
    default_view_scale              = default_view_scale,
    overlay_z_fixed                 = overlay_z_fixed,
    overlay_auto_reproject_to_wgs84 = overlay_auto_reproject_to_wgs84,
    overlay_point_color             = overlay_point_color,
    overlay_line_color              = overlay_line_color,
    overlay_point_size              = overlay_point_size,
    overlay_line_width              = overlay_line_width,
    show_north_arrow                = show_north_arrow,
    show_scale_bar                  = show_scale_bar,
    annotation_color                = annotation_color,
    annotation_line_width           = annotation_line_width,
)

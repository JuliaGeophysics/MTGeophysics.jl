using MTGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : ""
data_file  = length(ARGS) >= 2 ? ARGS[2] : ""
coordinate_system = length(ARGS) >= 3 ? ARGS[3] : "EPSG:32610"

shapefile_dir = normpath(@__DIR__, "cascadia", "shapefiles")

shapefiles = [
    (path = joinpath(@__DIR__, "cascadia", "shapefiles", "cascadia_coastline.shp"), enabled = true, color = :black, alpha = 0.95, point_size = 7, line_width = 1.2),
    (path = joinpath(@__DIR__, "cascadia", "shapefiles", "cascadia_state_boundaries.shp"), enabled = true, color = :black, alpha = 0.85, point_size = 7, line_width = 1.0),
    (path = joinpath(@__DIR__, "cascadia", "shapefiles", "cascadia_country_boundaries.shp"), enabled = true, color = :black, alpha = 0.85, point_size = 7, line_width = 1.0),
]

log10_scale       = true
colormap          = :Spectral
with_padding      = false
max_depth         = nothing
resistivity_range = (0.0, 5.0)
show_grid         = true
grid_color        = :black
grid_linewidth    = 0.5
grid_alpha        = 0.1
pad_tol           = 0.5

custom_extent = nothing

viewer_figsize = (1100, 950)
export_dpi     = 3
export_figsize = (1100, 900)

gis_output_dir = ""

PlotModelXY(model_file, data_file;
    crs               = coordinate_system,
    shapefiles        = shapefiles,
    log10_scale       = log10_scale,
    colormap          = colormap,
    with_padding      = with_padding,
    max_depth         = max_depth,
    resistivity_range = resistivity_range,
    show_grid         = show_grid,
    grid_color        = grid_color,
    grid_linewidth    = grid_linewidth,
    grid_alpha        = grid_alpha,
    pad_tol           = pad_tol,
    custom_extent     = custom_extent,
    viewer_figsize    = viewer_figsize,
    export_dpi        = export_dpi,
    export_figsize    = export_figsize,
    gis_output_dir    = gis_output_dir,
)

using MTGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : ""
data_file  = length(ARGS) >= 2 ? ARGS[2] : ""
coordinate_system = length(ARGS) >= 3 ? ARGS[3] : "model"

log10_scale       = true
colormap          = :Spectral
with_padding      = false
max_depth         = 50000
resistivity_range = (0.0, 4.0)
show_grid         = false
grid_color        = :black
grid_linewidth    = 0.5
grid_alpha        = 0.1
pad_tol           = 0.5

viewer_figsize = (1100, 950)
export_dpi     = 3
export_figsize = (1100, 900)

PlotModelYZ(model_file, data_file;
    crs               = coordinate_system,
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
    viewer_figsize    = viewer_figsize,
    export_dpi        = export_dpi,
    export_figsize    = export_figsize,
)

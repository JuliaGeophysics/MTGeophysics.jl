using GLMakie
using MTGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : ""

target_resistivity     = 1000.0
blend_previous_percent = 0

log10_scale       = true
colormap          = Reverse(:turbo)
with_padding      = true
max_depth         = nothing
resistivity_range = (0.0, 4.0)
show_grid         = true
grid_color        = :black
grid_linewidth    = 0.5
grid_alpha        = 0.3
pad_tol           = 0.5

viewer_figsize = (1100, 950)

EditModelByLayers(model_file;
    target_resistivity     = target_resistivity,
    blend_previous_percent = blend_previous_percent,
    log10_scale            = log10_scale,
    colormap               = colormap,
    with_padding           = with_padding,
    max_depth              = max_depth,
    resistivity_range      = resistivity_range,
    show_grid              = show_grid,
    grid_color             = grid_color,
    grid_linewidth         = grid_linewidth,
    grid_alpha             = grid_alpha,
    pad_tol                = pad_tol,
    viewer_figsize         = viewer_figsize,
)

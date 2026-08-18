using GLMakie
using MTGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : ""

replacement_resistivity = 10000.0
layers_above            = 2
layers_below            = 2
transition_layers       = 1
apply_to_all_depths     = false
depth_range             = (0.0, 50000.0)

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

viewer_figsize = (1200, 950)

EditModelByDrawing(model_file;
    replacement_resistivity = replacement_resistivity,
    layers_above            = layers_above,
    layers_below            = layers_below,
    transition_layers       = transition_layers,
    apply_to_all_depths     = apply_to_all_depths,
    depth_range             = depth_range,
    log10_scale             = log10_scale,
    colormap                = colormap,
    with_padding            = with_padding,
    max_depth               = max_depth,
    resistivity_range       = resistivity_range,
    show_grid               = show_grid,
    grid_color              = grid_color,
    grid_linewidth          = grid_linewidth,
    grid_alpha              = grid_alpha,
    pad_tol                 = pad_tol,
    viewer_figsize          = viewer_figsize,
)

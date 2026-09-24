# 2D inversion inputs from a data file: mesh, topography, water, covariance, mask and controls
# Author: @pankajkmishra
# Opens the GLMakie mesh window by default; mode "batch" (or MTGEO_MESH_MODE=batch) writes the inputs straight away
# Usage: julia --project=. examples/make_mesh2D.jl [data.dat] [topo.dat] [out_dir] [gui|batch]

using MTGeophysics

data_file = get(ARGS, 1, joinpath(@__DIR__, "data", "2D-IV", "data.dat"))
topo_file = get(ARGS, 2, joinpath(dirname(data_file), "topo.dat"))
out_dir   = get(ARGS, 3, joinpath(dirname(data_file), "mesh"))
mode      = Symbol(lowercase(get(ARGS, 4, get(ENV, "MTGEO_MESH_MODE", "gui"))))

# lakes and sea in the padding: profile range (m, the data's local y) and water level (m a.s.l.)
water = [(y_range = (-30_000.0, -18_000.0), level = 95.0)]

MakeMesh2D(data_file;
    out_dir,
    topo_path         = isfile(topo_file) ? topo_file : "",
    water             = isfile(topo_file) ? water : NamedTuple[],
    mode,
    inv_ctrl          = joinpath(@__DIR__, "ctrl", "2D", "InvCtrl.GN"),
    cell_width_frac   = 0.5,      # core cell width / median station spacing
    core_margin_cells = 4,        # core cells beyond the outer stations
    n_pad             = 12,       # padding cells each side
    pad_factor        = 1.3,
    first_layer_div   = 5.0,      # first layer = δ(f_max) / this
    vertical_factor   = 1.1,
    depth_mult        = 4.0,      # model depth = this × δ(f_min)
    air_layers        = 10,
    air_thickness     = 50_000.0,
    air_growth        = 2.0,
    dipole_length     = 100.0,
    cov_smoothing     = 0.3,
    water_resistivity = 100.0,
)

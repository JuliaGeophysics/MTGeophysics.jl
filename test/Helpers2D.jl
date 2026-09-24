# Shared 2D test fixtures
# Author: @pankajkmishra
# The hill mesh, mesh edits, Jyväskylä station coordinates and the topography case that several 2D test files use

using Test

# a small hill: no topographic air in the middle, one cell on the flanks, two outside;
# the receivers sit on three surface rows, which share node rows
function _hill_mesh()
    base = BuildMesh2D(frequencies = [0.3, 3.0], y_core_range = (-750.0, 750.0), y_core_cell = 250.0,
                       y_padding = 600.0, air_cells = 2, air_top = -1000.0, ground_layers = [50.0, 50.0, 100.0, 300.0, 900.0],
                       receiver_positions = [-625.0, -375.0, -125.0, 125.0, 400.0, 600.0])
    yc = (base.y_nodes[1:end-1] .+ base.y_nodes[2:end]) ./ 2
    topo = [abs(y) < 300 ? 0 : abs(y) < 600 ? 1 : 2 for y in yc]
    _remesh(base; topo_air = topo)
end

_remesh(m; kw...) = MT2DMesh(; (k => getfield(m, k) for k in fieldnames(MT2DMesh) if !haskey(kw, k))..., kw...)

# E-W stations near Jyväskylä at profile positions y
function _wgs(y)
    trans = MTGeophysics.Proj.Transformation(MTGeophysics._local_tm_proj_string(62.25, 25.75), "EPSG:4326"; always_xy = true)
    p = [trans(500_000.0 + yi, 0.0) for yi in y]
    (latitudes = [q[2] for q in p], longitudes = [q[1] for q in p])
end

# five stations over a hill, a narrow valley and a lake, on a 20 × 15 earth model
function _topo_case()
    stations = collect(-400.0:200.0:400.0)
    elevation(y) = 100 + 60exp(-(y / 300)^2) - 40exp(-((y + 800) / 100)^2)
    ty = collect(-1500.0:50.0:1500.0)
    topo = Topo2D(; _wgs(ty)..., elevations = elevation.(ty))
    data = (; receivers = stations, frequencies = [1.0, 10.0], z_positions = zeros(5),
              site_names = ["TK" * lpad(i, 2, '0') for i in 1:5], _wgs(stations)..., origin = [62.25, 25.75])
    model = ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = fill(100.0, 20),
                        z_cell_sizes = vcat(fill(20.0, 10), [50.0, 100.0, 200.0, 400.0, 800.0]),
                        resistivity = fill(100.0, 15, 20), n_air_cells = 0, origin = [0.0, -1000.0, 0.0],
                        rotation = 0.0, format = "LOGE")
    lake = (y_range = (-950.0, -700.0), level = 105.0)
    (; stations, elevation, ty, topo, data, model, lake)
end

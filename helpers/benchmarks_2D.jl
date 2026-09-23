# 2D COMEMI benchmark generator
# Author: @pankajkmishra
# Writes each case as a ModEM-style input set in examples/data/<case>: noisy synthetic data from a
# fine mesh, and start/prior models and cov.ctrl on a coarser inversion mesh
# The control files are shipped in examples/ctrl/2D; FwdCtrl there also sets the air of both meshes
# Usage: julia --project=. helpers/benchmarks_2D.jl [2D-I] [2D-II] [2D-III]      (default below)

using MTGeophysics
using Printf

const Proj = MTGeophysics.Proj

#---------- settings ----------

const BENCHMARK_CASES = Dict(
    "2D-I"   => "comemi2d_case1_dyke",
    "2D-II"  => "comemi2d_case2_resistive_blocks",
    "2D-III" => "comemi2d_case3_mixed",
)
# 2D-IV (2D-III with topography) takes over as the default once topography is supported
const DEFAULT_CASES = ["2D-III"]

# E-W profile near Jyväskylä, central Finland; local x north, y east, sites on the flat surface
const SURVEY = (
    origin         = (62.25, 25.75),                  # WGS84 lat, lon of local (0, 0)
    crs            = "EPSG:3067",                     # ETRS-TM35FIN, metric grid for the site positions
    site_prefix    = "JYV",
    receivers      = collect(-8000.0:1600.0:8000.0),
    frequencies    = 10 .^ range(-1, 3, length = 17), # 0.1-1000 Hz, 4 per decade
    error_fraction = 0.05,                            # impedance error, fraction of |Z|
    rng_seed       = 20260308,
)

# the synthetic data come from a finer mesh than the one the inversion uses; layers start at
# skin depth / first_layer_div for the highest frequency and grow by vertical_factor
const DATA_MESH = (
    y_core_range = (-9000.0, 9000.0), y_core_cell = 250.0, y_padding = 40_000.0, pad_factor = 1.3,
    background_resistivity = 100.0, first_layer_div = 10.0, vertical_factor = 1.05, depth_mult = 4.0,
)
const INVERSION_MESH = (
    y_core_range = (-9000.0, 9000.0), y_core_cell = 500.0, y_padding = 40_000.0, pad_factor = 1.3,
    background_resistivity = 100.0, first_layer_div = 5.0, vertical_factor = 1.1, depth_mult = 4.0,
)

const CTRL_DIR = joinpath(dirname(@__DIR__), "examples", "ctrl", "2D")
const FWD_PATH = joinpath(CTRL_DIR, "FwdCtrl")
const FWD = ReadFwdCtrl2D(FWD_PATH)

#---------- survey ----------

# WGS84 lat/lon of profile positions y (m east of the origin)
function survey_coordinates(y::AbstractVector{<:Real})
    to_grid = Proj.Transformation("EPSG:4326", SURVEY.crs; always_xy = true)
    to_wgs = Proj.Transformation(SURVEY.crs, "EPSG:4326"; always_xy = true)
    e0, n0 = to_grid(SURVEY.origin[2], SURVEY.origin[1])
    lonlat = [to_wgs(e0 + yi, n0) for yi in y]
    (latitudes = [p[2] for p in lonlat], longitudes = [p[1] for p in lonlat])
end

# zero impedances with fractional errors, the template ForwardSolve2D fills in
function survey_template()
    f = collect(Float64, SURVEY.frequencies)
    nf, ns = length(f), length(SURVEY.receivers)
    coords = survey_coordinates(SURVEY.receivers)
    nan = fill(NaN, nf, ns)
    DataFile2D(title = "2D survey template", periods = 1 ./ f, frequencies = f,
               site_names = [@sprintf("%s%03d", SURVEY.site_prefix, i) for i in 1:ns],
               receivers = SURVEY.receivers, x_positions = zeros(ns), z_positions = zeros(ns),
               z_xy = zeros(ComplexF64, nf, ns), z_xy_error = fill(SURVEY.error_fraction, nf, ns),
               z_yx = zeros(ComplexF64, nf, ns), z_yx_error = fill(SURVEY.error_fraction, nf, ns),
               z_xx = complex.(nan), z_xx_error = copy(nan), z_yy = complex.(nan), z_yy_error = copy(nan),
               rho_xy = copy(nan), phase_xy = copy(nan), rho_yx = copy(nan), phase_yx = copy(nan),
               latitudes = coords.latitudes, longitudes = coords.longitudes, origin = collect(SURVEY.origin))
end

function survey_mesh(m)
    layers = mt2d_geometric_layers(SURVEY.frequencies; m.background_resistivity, m.first_layer_div,
                                   m.vertical_factor, m.depth_mult)
    BuildMesh2D(; m.y_core_range, m.y_core_cell, m.y_padding, m.pad_factor, frequencies = SURVEY.frequencies,
                receiver_positions = SURVEY.receivers, ground_layers = layers, air_top = -FWD.air_thickness,
                air_cells = FWD.air_layers)
end

#---------- generator ----------

"""
    SaveBenchmarks2D(; output_root=examples/data, cases=DEFAULT_CASES) -> records

Write each COMEMI case (`"2D-I"`, `"2D-II"`, `"2D-III"`) into `output_root/<case>`:
- `model.true`: the true model on the fine data mesh (16 m surface layer), for plots only
- `data.dat`: noisy synthetic impedances from the fine mesh, errors 5% of |Z|
- `model.start`, `model.prior`: 100 ohm m halfspace on the inversion mesh (32 m surface layer)
- `cov.ctrl`: model covariance for the inversion mesh, every cell free

The forward and inversion controls come from `examples/ctrl/2D`.
"""
function SaveBenchmarks2D(;
    output_root::AbstractString = joinpath(dirname(@__DIR__), "examples", "data"),
    cases::AbstractVector{<:AbstractString} = DEFAULT_CASES,
)
    unknown = setdiff(cases, keys(BENCHMARK_CASES))
    isempty(unknown) || error("unknown cases $unknown, choose from $(sort(collect(keys(BENCHMARK_CASES))))")
    data_mesh, inv_mesh = survey_mesh(DATA_MESH), survey_mesh(INVERSION_MESH)
    truths = Dict(m.name => m for m in MTGeophysics.build_mt2d_comemi_models(data_mesh))
    halfspace = build_mt2d_halfspace_model(inv_mesh; background_resistivity = INVERSION_MESH.background_resistivity)
    nz, ny = length(inv_mesh.z_cell_sizes) - inv_mesh.n_air_cells, length(inv_mesh.y_cell_sizes)
    records = NamedTuple[]
    for case in cases
        dir = joinpath(output_root, case)
        mkpath(dir)
        paths = (
            case_dir = dir,
            true_model_path = WriteModel2D(joinpath(dir, "model.true"), data_mesh, truths[BENCHMARK_CASES[case]].resistivity),
            start_model_path = WriteModel2D(joinpath(dir, "model.start"), inv_mesh, halfspace),
            prior_model_path = WriteModel2D(joinpath(dir, "model.prior"), inv_mesh, halfspace),
            cov_path = WriteCov2D(joinpath(dir, "cov.ctrl"), Cov2D(nz, ny)),
        )
        data_path = mktempdir() do tmp
            template = write_data2d(joinpath(tmp, "data.template"), survey_template())
            ForwardSolve2D(paths.true_model_path, template, FWD_PATH; output_path = joinpath(dir, "data.dat"),
                           add_noise = true, rng_seed = SURVEY.rng_seed)
        end
        push!(records, merge(paths, (; data_path)))
    end
    records
end

#---------- script entry ----------

function main(args::AbstractVector{<:AbstractString} = ARGS)
    records = SaveBenchmarks2D(cases = isempty(args) ? DEFAULT_CASES : collect(args))
    for r in records
        println(r.case_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

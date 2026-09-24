# 2D COMEMI benchmark generator
# Author: @pankajkmishra
# Writes each case as a ModEM-style input set in examples/data/<case>: noisy synthetic data from a
# fine mesh, and start/prior models, cov.ctrl (GN, NLCG) and mask.ctrl (VFSA) on a coarser inversion mesh
# The control files are shipped in examples/ctrl/2D; FwdCtrl there also sets the air of both meshes
# 2D-IV is 2D-III under synthetic Finnish relief with a lake and its bathymetry among the stations, written with its
# topo.dat; its data mesh splits every inversion cell, so both meshes share one topography staircase and differ only
# below the ground
# Usage: julia --project=. helpers/benchmarks_2D.jl [2D-I] [2D-II] [2D-III] [2D-IV]      (default below)

using MTGeophysics
using Printf

const Proj = MTGeophysics.Proj

#---------- settings ----------

const BENCHMARK_CASES = Dict(
    "2D-I"   => "comemi2d_case1_dyke",
    "2D-II"  => "comemi2d_case2_resistive_blocks",
    "2D-III" => "comemi2d_case3_mixed",
    "2D-IV"  => "comemi2d_case3_mixed",
)
const TOPOGRAPHY_CASES = ["2D-IV"]
const DEFAULT_CASES = ["2D-IV"]

# E-W profile near Jyväskylä, central Finland; local x north, y east; sites every km at inversion cell centres,
# and the topographic cases leave out those on the lake
const SURVEY = (
    origin         = (62.25, 25.75),                  # WGS84 lat, lon of local (0, 0)
    crs            = "EPSG:3067",                     # ETRS-TM35FIN, metric grid for the site positions
    site_prefix    = "Fin",
    receivers      = collect(-8250.0:1000.0:8250.0),
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

# topography cases: the data mesh is the inversion mesh with every cell split y × z, air and water
# inherited from the parent cell; independent staircases leave a TM modelling floor near rms 100
const SHARED_SURFACE_SPLIT = (y = 2, z = 2)

# smooth relief of central Finland, 100-225 m a.s.l.: rolling till with a ridge, and a 3 km lake at 100 m a.s.l.
# over the buried conductor, with a steep western shore and two basins (58 and 36 m deep) split by a shoal;
# fresh humic lake water, 5 mS/m; elevations m a.s.l., the lake bottom from its bathymetry
const TOPOGRAPHY = (
    y = collect(-60_000.0:100.0:60_000.0),
    lake = (y_range = (-5500.0, -2500.0), level = 100.0),
    water_resistivity = 200.0,
)
land(y) = 140 + 40 * sin(2π * y / 11_000) + 45 * exp(-((y - 3000) / 1600)^2)
function bathymetry(y)
    a, b = TOPOGRAPHY.lake.y_range
    a < y < b || return 0.0
    sin(π * (y - a) / (b - a))^(1 / 3) * (8 + 50 * exp(-((y + 4700) / 600)^2) + 30 * exp(-((y + 3300) / 450)^2))
end
relief(y) = bathymetry(y) > 0 ? TOPOGRAPHY.lake.level - bathymetry(y) : land(y)

# stations of a case: none over the lake in the topographic ones
stations(topographic::Bool) = topographic ?
    filter(y -> !(TOPOGRAPHY.lake.y_range[1] <= y <= TOPOGRAPHY.lake.y_range[2]), SURVEY.receivers) : SURVEY.receivers

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

# zero impedances with fractional errors at stations y, the template ForwardSolve2D fills in
function survey_template(y::AbstractVector{<:Real} = SURVEY.receivers)
    f = collect(Float64, SURVEY.frequencies)
    nf, ns = length(f), length(y)
    coords = survey_coordinates(y)
    nan = fill(NaN, nf, ns)
    DataFile2D(title = "2D survey template", periods = 1 ./ f, frequencies = f,
               site_names = [@sprintf("%s%03d", SURVEY.site_prefix, i) for i in 1:ns],
               receivers = collect(Float64, y), x_positions = zeros(ns), z_positions = zeros(ns),
               z_xy = zeros(ComplexF64, nf, ns), z_xy_error = fill(SURVEY.error_fraction, nf, ns),
               z_yx = zeros(ComplexF64, nf, ns), z_yx_error = fill(SURVEY.error_fraction, nf, ns),
               z_xx = complex.(nan), z_xx_error = copy(nan), z_yy = complex.(nan), z_yy_error = copy(nan),
               rho_xy = copy(nan), phase_xy = copy(nan), rho_yx = copy(nan), phase_yx = copy(nan),
               latitudes = coords.latitudes, longitudes = coords.longitudes, origin = collect(SURVEY.origin))
end

function survey_mesh(m, y = SURVEY.receivers)
    layers = mt2d_geometric_layers(SURVEY.frequencies; m.background_resistivity, m.first_layer_div,
                                   m.vertical_factor, m.depth_mult)
    BuildMesh2D(; m.y_core_range, m.y_core_cell, m.y_padding, m.pad_factor, frequencies = SURVEY.frequencies,
                receiver_positions = collect(Float64, y), ground_layers = layers, air_top = -FWD.air_thickness,
                air_cells = FWD.air_layers)
end

function subdivide_mesh(m::MT2DMesh, ky::Int, kz::Int)
    na = m.n_air_cells
    dy = repeat(m.y_cell_sizes ./ ky, inner = ky)
    dz = vcat(m.z_cell_sizes[1:na], repeat(m.z_cell_sizes[na+1:end] ./ kz, inner = kz))
    MT2DMesh(y_nodes = m.y_nodes[1] .+ vcat(0.0, cumsum(dy)), z_nodes = m.z_nodes[1] .+ vcat(0.0, cumsum(dz)),
             y_cell_sizes = dy, z_cell_sizes = dz, receiver_positions = m.receiver_positions,
             frequencies = m.frequencies, n_air_cells = na, air_resistivity = m.air_resistivity)
end

# fine earth model with the coarse mask's air (1e17) and water cut in, cell by parent cell
function inherit_surface(ρ::AbstractMatrix, mask::AbstractMatrix{<:Integer}, ky::Int, kz::Int)
    out = copy(ρ)
    for iy in axes(out, 2), iz in axes(out, 1)
        parent = mask[(iz - 1) ÷ kz + 1, (iy - 1) ÷ ky + 1]
        parent == 0 && (out[iz, iy] = MTGeophysics.MT2D_AIR_TAG)
        parent == MTGeophysics.MT2D_MASK_WATER && (out[iz, iy] = TOPOGRAPHY.water_resistivity)
    end
    out
end

function survey_topography()
    coords = survey_coordinates(TOPOGRAPHY.y)
    Topo2D(latitudes = coords.latitudes, longitudes = coords.longitudes, elevations = relief.(TOPOGRAPHY.y))
end

# the model file with the topography cut in, and the data template with its station depths
function cut_topography(model_path, template, topo)
    t = Topography2D(ReadModel2D(model_path), template, topo; water = [TOPOGRAPHY.lake],
                     water_resistivity = TOPOGRAPHY.water_resistivity)
    WriteModel2D(model_path, t.model.y_cell_sizes, t.model.z_cell_sizes, t.model.resistivity)
    t
end

#---------- generator ----------

"""
    SaveBenchmarks2D(; output_root=examples/data, cases=DEFAULT_CASES) -> records

Write each COMEMI case (`"2D-I"`, `"2D-II"`, `"2D-III"`, `"2D-IV"`) into `output_root/<case>`:
- `model.true`: the true model on the fine data mesh, for plots only
- `data.dat`: noisy synthetic impedances from the fine mesh, errors 5% of |Z|
- `model.start`, `model.prior`: 100 ohm m halfspace on the inversion mesh (32 m surface layer)
- `cov.ctrl` (GN, NLCG) and `mask.ctrl` (VFSA): every inversion cell free

With topography (2D-IV) the models hold the topographic air (1e17 ohm m) and the lake
water, `cov.ctrl` and `mask.ctrl` give them mask 0 and 9, data Z is the depth below the
model top (the highest station), and `topo.dat` holds the relief with the lake bottom in WGS84 lat, lon,
elevation. Its data mesh splits every inversion cell `SHARED_SURFACE_SPLIT` and takes the
air and water of the parent cell, so the two meshes share one topography staircase.
The forward and inversion controls come from `examples/ctrl/2D`.
"""
function SaveBenchmarks2D(;
    output_root::AbstractString = joinpath(dirname(@__DIR__), "examples", "data"),
    cases::AbstractVector{<:AbstractString} = DEFAULT_CASES,
)
    unknown = setdiff(cases, keys(BENCHMARK_CASES))
    isempty(unknown) || error("unknown cases $unknown, choose from $(sort(collect(keys(BENCHMARK_CASES))))")
    inv_mesh = survey_mesh(INVERSION_MESH)
    ky, kz = SHARED_SURFACE_SPLIT
    halfspace = build_mt2d_halfspace_model(inv_mesh; background_resistivity = INVERSION_MESH.background_resistivity)
    nz, ny = length(inv_mesh.z_cell_sizes) - inv_mesh.n_air_cells, length(inv_mesh.y_cell_sizes)
    records = NamedTuple[]
    for case in cases
        topographic = case in TOPOGRAPHY_CASES
        data_mesh = topographic ? subdivide_mesh(inv_mesh, ky, kz) : survey_mesh(DATA_MESH)
        truth = only(m for m in MTGeophysics.build_mt2d_comemi_models(data_mesh) if m.name == BENCHMARK_CASES[case])
        dir = joinpath(output_root, case)
        mkpath(dir)
        paths = (
            case_dir = dir,
            true_model_path = WriteModel2D(joinpath(dir, "model.true"), data_mesh, truth.resistivity),
            start_model_path = WriteModel2D(joinpath(dir, "model.start"), inv_mesh, halfspace),
            prior_model_path = WriteModel2D(joinpath(dir, "model.prior"), inv_mesh, halfspace),
            cov_path = WriteCov2D(joinpath(dir, "cov.ctrl"), Cov2D(nz, ny)),
            mask_path = WriteMask2D(joinpath(dir, "mask.ctrl"), ones(Int, nz, ny)),
        )
        template = survey_template(stations(topographic))
        if topographic
            topo = survey_topography()
            WriteTopo2D(joinpath(dir, "topo.dat"), topo)
            t = cut_topography(paths.start_model_path, template, topo)
            cut_topography(paths.prior_model_path, template, topo)
            template = t.data
            c = Cov2D(nz, ny)
            WriteCov2D(paths.cov_path, Cov2D(sy = c.sy, sz = c.sz, n_smooth = c.n_smooth, mask = t.mask))
            WriteMask2D(paths.mask_path, t.mask)
            fine = ReadModel2D(paths.true_model_path)
            WriteModel2D(paths.true_model_path, fine.y_cell_sizes, fine.z_cell_sizes,
                         inherit_surface(fine.resistivity, t.mask, ky, kz))
            e = relief.(stations(true))
            @printf("%s: datum %.1f m a.s.l., station relief %.1f m, %d air and %d water cells in the inversion model\n",
                    case, t.datum, maximum(e) - minimum(e), count(==(0), t.mask), count(==(9), t.mask))
        end
        data_path = mktempdir() do tmp
            ForwardSolve2D(paths.true_model_path, write_data2d(joinpath(tmp, "data.template"), template), FWD_PATH;
                           output_path = joinpath(dir, "data.dat"), add_noise = true, rng_seed = SURVEY.rng_seed)
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

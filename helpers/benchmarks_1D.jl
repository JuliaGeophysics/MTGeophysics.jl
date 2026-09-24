# 1D benchmark generator
# Author: @pankajkmishra
# Writes each case into examples/data/<case>: noisy synthetic data from the exact layered truth, and the truth for
# plots; the inversion needs only the data and a control file from examples/ctrl/1D (MakeMesh1D lays out the layers)
# Usage: julia --project=. helpers/benchmarks_1D.jl [1D-I]

using MTGeophysics
using Printf

#---------- settings ----------

# 1D-I: five layers under one site near Jyväskylä
const CASES_1D = Dict(
    "1D-I" => (thicknesses = [120.0, 280.0, 650.0, 1400.0], resistivities = [100.0, 20.0, 350.0, 40.0, 800.0]),
)
const SURVEY_1D = (
    site = "UA01", latitude = 62.25, longitude = 25.75,
    frequencies = 10 .^ range(-2.5, 3, length = 23),   # 0.003-1000 Hz, 4 per decade
    error_fraction = 0.05, rng_seed = 20260308,
)

#---------- generator ----------

# zero impedances with fractional errors at the one site, the template ForwardSolve1D fills in
function survey_template_1d()
    f = collect(Float64, SURVEY_1D.frequencies)
    nf = length(f)
    nan = fill(NaN, nf, 1)
    DataFile2D(title = "1D survey template", periods = 1 ./ f, frequencies = f, site_names = [SURVEY_1D.site],
               receivers = [0.0], x_positions = [0.0], z_positions = [0.0],
               z_xy = zeros(ComplexF64, nf, 1), z_xy_error = fill(SURVEY_1D.error_fraction, nf, 1),
               z_yx = zeros(ComplexF64, nf, 1), z_yx_error = fill(SURVEY_1D.error_fraction, nf, 1),
               z_xx = complex.(nan), z_xx_error = copy(nan), z_yy = complex.(nan), z_yy_error = copy(nan),
               rho_xy = copy(nan), phase_xy = copy(nan), rho_yx = copy(nan), phase_yx = copy(nan),
               latitudes = [SURVEY_1D.latitude], longitudes = [SURVEY_1D.longitude],
               origin = [SURVEY_1D.latitude, SURVEY_1D.longitude])
end

"""
    SaveBenchmarks1D(; output_root=examples/data, cases=["1D-I"]) -> records

Write each 1D case into `output_root/<case>`:
- `model.true`: the layered truth (last layer the halfspace), for plots only
- `data.dat`: noisy synthetic impedances, errors 5% of |Z|, ZXY = Z and ZYX = -Z
"""
function SaveBenchmarks1D(; output_root::AbstractString = joinpath(dirname(@__DIR__), "examples", "data"),
                          cases::AbstractVector{<:AbstractString} = ["1D-I"])
    unknown = setdiff(cases, keys(CASES_1D))
    isempty(unknown) || error("unknown cases $unknown, choose from $(sort(collect(keys(CASES_1D))))")
    records = NamedTuple[]
    for case in cases
        c = CASES_1D[case]
        dir = mkpath(joinpath(output_root, case))
        # the last layer is the halfspace; its file thickness only sets how deep the plots draw it
        truth = vcat(c.thicknesses, 10 * c.thicknesses[end])
        paths = (
            case_dir = dir,
            true_model_path = WriteModel2D(joinpath(dir, "model.true"), [1.0], truth, reshape(c.resistivities, :, 1)),
        )
        data_path = mktempdir() do tmp
            ForwardSolve1D(paths.true_model_path, write_data2d(joinpath(tmp, "data.template"), survey_template_1d());
                           output_path = joinpath(dir, "data.dat"), add_noise = true, rng_seed = SURVEY_1D.rng_seed)
        end
        push!(records, merge(paths, (; data_path)))
    end
    records
end

function main(args::AbstractVector{<:AbstractString} = ARGS)
    for r in SaveBenchmarks1D(cases = isempty(args) ? ["1D-I"] : collect(args))
        println(r.case_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# 2D COMEMI benchmark generator
# Author: @pankajkmishra
# Writes true model, halfspace start model, error template, and noisy observed data per case
# Usage: julia --project=. helpers/benchmarks_2D.jl      (settings below)

using MTGeophysics

#---------- script settings ----------

# cases: 1 = dyke, 2 = resistive blocks, 3 = mixed
const BENCHMARK_CASES = [3]

# survey and skin-depth mesh used when run as a script; the inversion examples read
# the files written here, so this is the one place the synthetic survey is defined
const BENCHMARK_MESH = (
    frequencies            = 10 .^ range(-1, 3, length = 17),   # 0.1-1000 Hz, 4 per decade
    receiver_positions     = collect(-8000.0:1600.0:8000.0),
    background_resistivity = 100.0,
    y_core_range           = (-9000.0, 9000.0),
    y_core_cell            = 500.0,
    y_padding              = 40_000.0,
    pad_factor             = 1.3,
    air_top                = -40_000.0,
    air_cells              = 8,
    max_core_layers        = 50,        # uniform dz down to 1 skin depth of f_min
)

#---------- generator ----------

"""
    SaveBenchmarks2D(; output_root=joinpath(dirname(@__DIR__), "examples"), cases=[1],
                     mesh=build_default_mt2d_mesh(), start_resistivity=100.0,
                     error_fraction=0.05, rng_seed=20260308)

Inputs:
- `output_root`: Parent directory under which per-case benchmark directories are created.
- `cases`: COMEMI2D cases to generate (1 = dyke, 2 = resistive blocks, 3 = mixed).
- `mesh`: Mesh and survey (frequencies, receivers) the models and data live on.
- `start_resistivity`: Halfspace resistivity of the start model.
- `error_fraction`: Impedance error as a fraction of |Z|.
- `rng_seed`: Noise seed.

Output:
- Vector of named tuples with the written true-model, start-model, reference-data, and observed-data paths.

Description:
- Writes each case into its own directory, e.g. `examples/0COMEMI2D-III/Comemi2D3.{true,ini,ref,obs}`.
- The COMEMI bodies are defined in metres, so any mesh that covers them works.
"""
function SaveBenchmarks2D(;
    output_root::AbstractString = joinpath(dirname(@__DIR__), "examples"),
    cases::AbstractVector{Int} = [1],
    mesh::MT2DMesh = build_default_mt2d_mesh(),
    start_resistivity::Real = 100.0,
    error_fraction::Real = 0.05,
    rng_seed::Integer = 20260308,
)
    all_cases = [
        ("comemi2d_case1_dyke", "Comemi2D1", "0COMEMI2D-I"),
        ("comemi2d_case2_resistive_blocks", "Comemi2D2", "0COMEMI2D-II"),
        ("comemi2d_case3_mixed", "Comemi2D3", "0COMEMI2D-III"),
    ]
    all(c -> 1 <= c <= length(all_cases), cases) && !isempty(cases) ||
        error("cases must be a non-empty subset of 1:$(length(all_cases)), got $cases")

    models = Dict(m.name => m for m in MTGeophysics.build_mt2d_comemi_models(mesh))
    saved = NamedTuple[]
    for (case_key, case_label, dir_name) in all_cases[cases]
        case_dir = joinpath(output_root, dir_name)
        mkpath(case_dir)
        model = models[case_key]

        model_path = write_model2d(joinpath(case_dir, "$(case_label).true"), mesh, model.resistivity;
                                   title = model.label)
        start_model_path = write_model2d(joinpath(case_dir, "$(case_label).ini"), mesh,
            build_mt2d_halfspace_model(mesh; background_resistivity = start_resistivity);
            title = "Halfspace starting model for $(case_label)")
        reference_path = write_mt2d_data_template(joinpath(case_dir, "$(case_label).ref"), mesh;
                                                  impedance_error_fraction = error_fraction)
        observed_path = ForwardSolve2D(model_path, reference_path; add_noise = true,
                                       output_path = joinpath(case_dir, "$(case_label).obs"), rng_seed = rng_seed)
        push!(saved, (; case_dir, model_path, start_model_path, reference_path, observed_path))
    end
    saved
end

"""
    save_benchmarks2d(; kwargs...)

Lowercase alias for `SaveBenchmarks2D`.
"""
save_benchmarks2d(; kwargs...) = SaveBenchmarks2D(; kwargs...)

#---------- script entry ----------

function main(args::AbstractVector{<:AbstractString} = ARGS)
    isempty(args) || error("usage: julia --project=. helpers/benchmarks_2D.jl")
    mesh = BuildMesh2D(; BENCHMARK_MESH...)
    saved = SaveBenchmarks2D(cases = BENCHMARK_CASES, mesh = mesh,
                             start_resistivity = BENCHMARK_MESH.background_resistivity)
    println("SavedCases = ", length(saved))
    for s in saved
        println("  ", s.case_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

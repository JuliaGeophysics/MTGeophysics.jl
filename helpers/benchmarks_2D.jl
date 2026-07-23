# This helper script writes the standard 2D benchmark true-model, start-model, reference-data, and observed-data files.

using MTGeophysics

# cases to generate when run as a script: 1 = dyke, 2 = resistive blocks, 3 = mixed
const BENCHMARK_CASES = [1]

"""
    SaveBenchmarks2D(; output_root=joinpath(dirname(@__DIR__), "examples"), cases=[1])

Inputs:
- `output_root`: Parent directory under which per-case benchmark directories are created.
- `cases`: Which COMEMI2D cases to generate (1 = dyke, 2 = resistive blocks, 3 = mixed).

Output:
- Vector of named tuples containing the written true-model, start-model, reference-data, and observed-data paths.

Description:
- Writes each selected 2D benchmark into its own sub-directory (e.g. `examples/0COMEMI2D-I/`).
- Writes a homogeneous halfspace starting model (`.ini`) on the same grid as each true model.
"""
function SaveBenchmarks2D(;
    output_root::AbstractString = joinpath(dirname(@__DIR__), "examples"),
    cases::AbstractVector{Int} = [1],
)
    all_cases = [
        ("comemi2d_case1_dyke", "Comemi2D1", "0COMEMI2D-I"),
        ("comemi2d_case2_resistive_blocks", "Comemi2D2", "0COMEMI2D-II"),
        ("comemi2d_case3_mixed", "Comemi2D3", "0COMEMI2D-III"),
    ]
    all(c -> 1 <= c <= length(all_cases), cases) && !isempty(cases) ||
        error("cases must be a non-empty subset of 1:$(length(all_cases)), got $cases")
    case_names = all_cases[cases]

    # Build mesh once; write true models into first case dir as temp staging.
    first_case_dir = joinpath(output_root, case_names[1][3])
    mkpath(first_case_dir)
    mesh_result = MakeMesh2D(output_dir = first_case_dir)

    saved = NamedTuple[]
    for (case_key, case_label, dir_name) in case_names
        case_dir = joinpath(output_root, dir_name)
        mkpath(case_dir)

        source_model_path = mesh_result.model_paths[case_key]
        target_model_path = joinpath(case_dir, "$(case_label).true")
        abspath(source_model_path) == abspath(target_model_path) || cp(source_model_path, target_model_path; force = true)

        start_model_path = joinpath(case_dir, "$(case_label).ini")
        halfspace_resistivity = build_mt2d_halfspace_model(mesh_result.mesh)
        write_model2d(start_model_path, mesh_result.mesh, halfspace_resistivity; title = "Halfspace starting model for $(case_label)")

        reference_path = write_mt2d_data_template(
            joinpath(case_dir, "$(case_label).ref"),
            mesh_result.mesh;
            impedance_error_fraction = 0.05,
        )
        observed_path = ForwardSolve2D(
            target_model_path,
            reference_path;
            add_noise = true,
            output_path = joinpath(case_dir, "$(case_label).obs"),
            rng_seed = 20260308,
        )
        push!(saved, (
            case_dir = case_dir,
            model_path = target_model_path,
            start_model_path = start_model_path,
            reference_path = reference_path,
            observed_path = observed_path,
        ))
    end

    # Clean up stale model files that MakeMesh2D staged into the first case dir,
    # including .true files of cases that were not selected.
    for (case_key, case_label, dir_name) in all_cases
        stale = joinpath(first_case_dir, "$(case_label).true")
        target = joinpath(output_root, dir_name, "$(case_label).true")
        keep = abspath(stale) == abspath(target) && (case_key, case_label, dir_name) in case_names
        if !keep && isfile(stale)
            rm(stale)
        end
    end

    saved
end

"""
    save_benchmarks2d(; kwargs...)

Inputs:
- Keyword arguments accepted by `SaveBenchmarks2D`.

Output:
- Benchmark file records.

Description:
- Lowercase alias for `SaveBenchmarks2D`.
"""
save_benchmarks2d(; kwargs...) = SaveBenchmarks2D(; kwargs...)

function main(args::AbstractVector{<:AbstractString} = ARGS)
    isempty(args) || error("usage: julia --project=. helpers/benchmarks_2D.jl")
    saved = SaveBenchmarks2D(cases = BENCHMARK_CASES)
    println("SavedCases = ", length(saved))
    for s in saved
        println("  ", s.case_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

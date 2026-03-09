# This helper script writes the standard 2D benchmark true-model, start-model, reference-data, and observed-data files.

using MTGeophysics

"""
    SaveBenchmarks2D(; output_root=joinpath(dirname(@__DIR__), "Examples"))

Inputs:
- `output_root`: Parent directory under which per-case benchmark directories are created.

Output:
- Vector of named tuples containing the written true-model, start-model, reference-data, and observed-data paths.

Description:
- Writes each 2D benchmark into its own sub-directory (e.g. `Examples/0COMEMI2D-I/`).
- Writes a homogeneous halfspace starting model (`.ini`) on the same grid as each true model.
"""
function SaveBenchmarks2D(;
    output_root::AbstractString = joinpath(dirname(@__DIR__), "Examples"),
)
    case_names = [
        ("comemi2d_case1_dyke", "Comemi2D1", "0COMEMI2D-I"),
        ("comemi2d_case2_resistive_blocks", "Comemi2D2", "0COMEMI2D-II"),
        ("comemi2d_case3_mixed", "Comemi2D3", "0COMEMI2D-III"),
    ]

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

    # Clean up any stale model files that MakeMesh2D wrote into the first case dir
    # for the other cases (case2, case3 .rho files).
    for (case_key, case_label, dir_name) in case_names
        stale = joinpath(first_case_dir, "$(case_label).true")
        target = joinpath(output_root, dir_name, "$(case_label).true")
        if abspath(stale) != abspath(target) && isfile(stale)
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
    isempty(args) || error("usage: julia --project=. Helpers/benchmarks_2d.jl")
    saved = SaveBenchmarks2D()
    println("SavedCases = ", length(saved))
    for s in saved
        println("  ", s.case_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

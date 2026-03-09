# This helper script writes the standard 1D benchmark true-model, reference-data, and observed-data files.

using MTGeophysics

"""
    SaveBenchmarks1D(; output_root=joinpath(dirname(@__DIR__), "Examples"))

Inputs:
- `output_root`: Parent directory under which per-case benchmark directories are created.

Output:
- Vector of named tuples containing the written true-model, reference-data, and observed-data paths.

Description:
- Writes each 1D benchmark into its own sub-directory (e.g. `Examples/0Layered1D/`).
"""
function SaveBenchmarks1D(;
    output_root::AbstractString = joinpath(dirname(@__DIR__), "Examples"),
)
    cases = [
        (
            label = "Layered1D",
            dir_name = "0Layered1D",
            thicknesses = [120.0, 280.0, 650.0, 1400.0],
            resistivities = [100.0, 20.0, 350.0, 40.0, 800.0],
            frequencies = collect(10 .^ range(-2.5, 2.5, length = 24)),
        ),
    ]

    saved = NamedTuple[]
    for case in cases
        case_dir = joinpath(output_root, case.dir_name)
        mkpath(case_dir)

        model_path = MakeMesh1D(
            output_dir = case_dir,
            model_name = "$(case.label).true",
            thicknesses = case.thicknesses,
            resistivities = case.resistivities,
        ).model_path
        reference_path = write_mt1d_data_template(
            joinpath(case_dir, "$(case.label).ref"),
            case.frequencies;
            rho_error_fraction = 0.05,
            phase_error_deg = 1.5,
        )
        specification = load_mt1d_data_spec(reference_path)
        response = solve_mt1d_analytical(case.frequencies, case.resistivities, case.thicknesses)
        observed = MTGeophysics._apply_mt1d_noise(response, specification; rng_seed = 20260308)
        observed_path = write_mt1d_observed_data(
            joinpath(case_dir, "$(case.label).obs"),
            specification,
            observed,
        )
        push!(saved, (
            case_dir = case_dir,
            model_path = model_path,
            reference_path = reference_path,
            observed_path = observed_path,
        ))
    end

    saved
end

"""
    save_benchmarks1d(; kwargs...)

Inputs:
- Keyword arguments accepted by `SaveBenchmarks1D`.

Output:
- Benchmark file records.

Description:
- Lowercase alias for `SaveBenchmarks1D`.
"""
save_benchmarks1d(; kwargs...) = SaveBenchmarks1D(; kwargs...)

function main(args::AbstractVector{<:AbstractString} = ARGS)
    isempty(args) || error("usage: julia --project=. Helpers/benchmarks_1d.jl")
    saved = SaveBenchmarks1D()
    println("SavedCases = ", length(saved))
    for s in saved
        println("  ", s.case_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

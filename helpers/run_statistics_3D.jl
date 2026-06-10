# This helper script recomputes the 3D ensemble statistics from saved chain models.
#
# Usage:
#   julia --project=. helpers/run_statistics_3D.jl <run_dir>
#   julia --project=. helpers/run_statistics_3D.jl <run_dir> "^chain.*\.rho$"
#
# VFSA3DMT writes the per-chain best models (best_model_chainNN.rho) directly into
# the starting model's directory, so <run_dir> is that directory. The optional
# second argument is a regex pattern to match model filenames
# (default: "^best_model.*\.rho$").

using MTGeophysics

function main(args::AbstractVector{<:AbstractString} = ARGS)
    1 <= length(args) <= 2 || error(
        "usage: julia --project=. helpers/run_statistics_3D.jl <run_dir> [pattern]")

    run_dir = args[1]
    isdir(run_dir) || error("Directory not found: $run_dir")

    pattern = length(args) >= 2 ? Regex(args[2]) : r"^best_model.*\.rho$"

    println("Computing 3D ensemble statistics...")
    println("  Directory: $run_dir")
    println("  Pattern:   $pattern")

    mean_path, median_path, std_path = AnalyseEnsemble3D(run_dir; pattern = pattern)

    println()
    println("Done.")
    println("  Mean   model: $mean_path")
    println("  Median model: $median_path")
    println("  Std    model: $std_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

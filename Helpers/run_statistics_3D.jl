# This helper script recomputes the 3D ensemble statistics from saved chain models.
#
# Usage:
#   julia --project=. Helpers/run_statistics_3d.jl <run_dir>
#   julia --project=. Helpers/run_statistics_3d.jl <run_dir> "^best_model.*\.rho$"
#
# The first argument is the VFSA3DMT output directory containing the .rho files.
# The optional second argument is a regex pattern to match model filenames
# (default: "^chain.*\.rho$").

using MTGeophysics

function main(args::AbstractVector{<:AbstractString} = ARGS)
    1 <= length(args) <= 2 || error(
        "usage: julia --project=. Helpers/run_statistics_3d.jl <run_dir> [pattern]")

    run_dir = args[1]
    isdir(run_dir) || error("Directory not found: $run_dir")

    pattern = length(args) >= 2 ? Regex(args[2]) : r"^chain.*\.rho$"

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

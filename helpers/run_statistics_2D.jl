# Recompute the 2D VFSA ensemble statistics from the saved best model of each chain
# Author: @pankajkmishra
# Usage: julia --project=. helpers/run_statistics_2D.jl <run_dir>

using MTGeophysics

function main(args::AbstractVector{<:AbstractString} = ARGS)
    length(args) == 1 || error("usage: julia --project=. helpers/run_statistics_2D.jl <run_dir>")
    stats = AnalyseEnsemble2D(args[1])
    println("Chains = ", length(stats.chains))
    foreach(p -> println("  ", p), stats.paths)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# This helper script recomputes the 2D ensemble statistics from the saved best model of each chain.

using MTGeophysics

function main(args::AbstractVector{<:AbstractString} = ARGS)
    1 <= length(args) <= 2 || error("usage: julia --project=. Helpers/run_statistics_2d.jl <run_dir> [observed_data_path]")
    run_dir = args[1]
    observed_data_path = length(args) == 2 ? args[2] : nothing
    stats = AnalyseEnsemble2D(run_dir; observed_data_path = observed_data_path, copy_script = false)
    println("Summary = ", stats.summary_path)
    println("MeanModel = ", stats.mean_model_path)
    println("MedianModel = ", stats.median_model_path)
    println("StdModel = ", stats.uncertainty_table_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

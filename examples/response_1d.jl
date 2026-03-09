# Compute the 1D forward response and plot observed vs predicted data.
# Usage: julia --project=. Examples/response_1d.jl <model_path> <reference_path>

using MTGeophysics

model_path = get(ARGS, 1, normpath(@__DIR__, "0Layered1D", "Layered1D.true"))
ref_path   = get(ARGS, 2, normpath(@__DIR__, "0Layered1D", "Layered1D.ref"))

pred_path = ForwardSolve1D(model_path, ref_path)
plot_path = PlotData1D(ref_path, pred_path; output_path = joinpath(dirname(pred_path), "DataPlot1D.png"))

println("Predicted = ", pred_path)
println("Plot      = ", plot_path)

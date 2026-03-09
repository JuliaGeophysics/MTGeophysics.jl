# Compute the 2D forward response and plot observed vs predicted data.
# Usage: julia --project=. Examples/response_2d.jl <model_path> <reference_path>

using MTGeophysics

model_path = get(ARGS, 1, normpath(@__DIR__, "0COMEMI2D-I", "Comemi2D1.true"))
ref_path   = get(ARGS, 2, normpath(@__DIR__, "0COMEMI2D-I", "Comemi2D1.ref"))

pred_path = ForwardSolve2D(model_path, ref_path)
out_dir   = dirname(pred_path)
plots     = PlotData2D(pred_path;
                maps_output_path   = joinpath(out_dir, "DataMaps2D.png"),
                curves_output_path = joinpath(out_dir, "DataCurves2D.png"))

println("Predicted = ", pred_path)
println("Maps      = ", plots.maps_output_path)
println("Curves    = ", plots.curves_output_path)

# 1D MT forward modelling: a one-column model and a data file (sites, periods, errors) in, data.pred out
# Author: @pankajkmishra
# Generate the benchmark once: julia --project=. helpers/benchmarks_1D.jl
# Usage: julia --project=. examples/run_fwd1D.jl [model.rho data.dat]

using MTGeophysics

const CASE_DIR = joinpath(@__DIR__, "data", "1D-I")
const MODE = :XYYX                 # XY, YX, XYYX or DET

model_path, data_path = length(ARGS) == 2 ? ARGS : (joinpath(CASE_DIR, "model.true"), joinpath(CASE_DIR, "data.dat"))
isfile(model_path) || error("missing $model_path; run julia --project=. helpers/benchmarks_1D.jl first")

pred_path = ForwardSolve1D(model_path, data_path; mode = MODE)
plot_path = PlotData2D(data_path; predicted_path = pred_path, output_path = joinpath(dirname(pred_path), "DataFit1D.png"),
                       names = MODE == :XY ? (TE = "XY", TM = nothing) : MODE == :YX ? (TE = nothing, TM = "YX") : (TE = "XY", TM = "YX"))
println("Predicted = ", pred_path)
println("Plot      = ", plot_path)

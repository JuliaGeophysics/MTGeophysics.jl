# 2D MT forward run, ModEM style: model + data file (sites, periods, errors) + FwdCtrl -> data.pred
# Author: @pankajkmishra
# data.pred is written next to the data file, with a data-fit plot beside it
# FwdCtrl is shipped in examples/ctrl/2D; generate the models and data once: julia --project=. helpers/benchmarks_2D.jl
# Usage: julia --project=. examples/run_fwd2D.jl [model.rho data.dat FwdCtrl]

using MTGeophysics

const CASE_DIR = joinpath(@__DIR__, "data", "2D-IV")
const CTRL_DIR = joinpath(@__DIR__, "ctrl", "2D")

inputs = if isempty(ARGS)
    [joinpath(CASE_DIR, "model.true"), joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "FwdCtrl")]
elseif length(ARGS) == 3
    ARGS
else
    error("usage: julia --project=. examples/run_fwd2D.jl model.rho data.dat FwdCtrl")
end
all(isfile, inputs) || error("missing inputs $(filter(!isfile, inputs)); run julia --project=. helpers/benchmarks_2D.jl first")

pred_path = ForwardSolve2D(inputs...)
plot_path = PlotData2D(inputs[2]; predicted_path = pred_path, output_path = joinpath(dirname(pred_path), "DataFit.png"))

println("Predicted : ", pred_path)
println("Plot      : ", plot_path)

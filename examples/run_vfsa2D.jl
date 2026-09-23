# 2D MT VFSA inversion of a benchmark, through the same six-file front end as run_inv2D.jl
# Author: @pankajkmishra
# Uses examples/ctrl/2D/InvCtrl.VFSA; results and plots go to run_YYYYmmdd_HHMMSS/ next to the data, chains in its vfsa/
# Generate the models and data once: julia --project=. helpers/benchmarks_2D.jl
# Usage: julia --project=. examples/run_vfsa2D.jl

using MTGeophysics
using Printf

const CASE_DIR = joinpath(@__DIR__, "data", "2D-III")
const CTRL_DIR = joinpath(@__DIR__, "ctrl", "2D")

inputs = [joinpath(CASE_DIR, "model.start"), joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "FwdCtrl"),
          joinpath(CTRL_DIR, "InvCtrl.VFSA"), joinpath(CASE_DIR, "cov.ctrl"), joinpath(CASE_DIR, "model.prior")]
all(isfile, inputs) || error("missing inputs $(filter(!isfile, inputs)); run julia --project=. helpers/benchmarks_2D.jl first")
true_model = joinpath(CASE_DIR, "model.true")

elapsed = @elapsed run = Invert2D(inputs...)
plots = PlotInversion2D(run; true_model_path = isfile(true_model) ? true_model : nothing, maximum_depth_km = 10.0)

@printf("Final RMS : %.3f, %.1f s\n", run.rms, elapsed)
println("Outputs   : ", run.run_dir)

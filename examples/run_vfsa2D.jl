# 2D MT VFSA inversion from five files: no covariance and no prior, mask.ctrl says which cells move
# Author: @pankajkmishra
# Results and plots go to run_YYYYmmdd_HHMMSS/ next to the data, chains and ensemble in its vfsa/
# Controls are shipped in examples/ctrl/2D; generate the models and data once: julia --project=. helpers/benchmarks_2D.jl
# Usage: julia --project=. examples/run_vfsa2D.jl
#        julia --project=. examples/run_vfsa2D.jl model.start data.dat FwdCtrl InvCtrl.VFSA mask.ctrl

using MTGeophysics
using Printf

const CASE_DIR = joinpath(@__DIR__, "data", "2D-IV")
const CTRL_DIR = joinpath(@__DIR__, "ctrl", "2D")
const PLOT_DEPTH_KM = 10.0         # depth limit of the model plots
const LOG10_RHO_RANGE = (0.0, 3.5)

inputs = if isempty(ARGS)
    [joinpath(CASE_DIR, "model.start"), joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "FwdCtrl"),
     joinpath(CTRL_DIR, "InvCtrl.VFSA"), joinpath(CASE_DIR, "mask.ctrl")]
elseif length(ARGS) == 5
    ARGS
else
    error("usage: julia --project=. examples/run_vfsa2D.jl [model.start data.dat FwdCtrl InvCtrl.VFSA mask.ctrl]")
end
all(isfile, inputs) || error("missing inputs $(filter(!isfile, inputs)); run julia --project=. helpers/benchmarks_2D.jl first")
true_model = joinpath(dirname(inputs[2]), "model.true")   # plots only, skipped when absent

elapsed = @elapsed run = VFSA2D(inputs...)
plots = PlotInversion2D(run; true_model_path = isfile(true_model) ? true_model : nothing,
                        maximum_depth_km = PLOT_DEPTH_KM, resistivity_log10_range = LOG10_RHO_RANGE)

println()
@printf("Chains      : %d, best chain %d at RMS %.3f\n", length(run.vfsa.chains), run.vfsa.best_chain, run.vfsa.best_rms)
@printf("Mean RMS    : %.3f (%s), %.1f s\n", run.rms, run.reason, elapsed)
println("Outputs     : ", run.run_dir, "  (model.rho = ensemble mean, data.pred, vfsa/)")
foreach(p -> println("  ", relpath(p, run.run_dir)), plots)

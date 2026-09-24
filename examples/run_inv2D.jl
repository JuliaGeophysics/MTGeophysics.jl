# 2D MT deterministic inversion (GN or NLCG), ModEM style, from six files; VFSA runs from run_vfsa2D.jl
# Author: @pankajkmishra
# The algorithm comes from the inversion control; results and plots go to run_YYYYmmdd_HHMMSS/ next to the data
# Controls are shipped in examples/ctrl/2D; generate the models and data once: julia --project=. helpers/benchmarks_2D.jl
# Usage: julia --project=. examples/run_inv2D.jl [GN|NLCG]
#        julia --project=. examples/run_inv2D.jl model.start data.dat FwdCtrl InvCtrl cov.ctrl model.prior

using MTGeophysics
using Printf

const CASE_DIR = joinpath(@__DIR__, "data", "2D-IV")
const CTRL_DIR = joinpath(@__DIR__, "ctrl", "2D")
const PLOT_DEPTH_KM = 10.0         # depth limit of the model plots
const LOG10_RHO_RANGE = (0.0, 3.5)

inputs = if length(ARGS) <= 1
    algorithm = uppercase(isempty(ARGS) ? "GN" : ARGS[1])
    [joinpath(CASE_DIR, "model.start"), joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "FwdCtrl"),
     joinpath(CTRL_DIR, "InvCtrl.$algorithm"), joinpath(CASE_DIR, "cov.ctrl"), joinpath(CASE_DIR, "model.prior")]
elseif length(ARGS) == 6
    ARGS
else
    error("usage: julia --project=. examples/run_inv2D.jl [GN|NLCG] or model.start data.dat FwdCtrl InvCtrl cov.ctrl model.prior")
end
all(isfile, inputs) || error("missing inputs $(filter(!isfile, inputs)); run julia --project=. helpers/benchmarks_2D.jl first")
true_model = joinpath(dirname(inputs[2]), "model.true")   # plots only, skipped when absent

elapsed = @elapsed run = Invert2D(inputs...)
plots = PlotInversion2D(run; true_model_path = isfile(true_model) ? true_model : nothing,
                        maximum_depth_km = PLOT_DEPTH_KM, resistivity_log10_range = LOG10_RHO_RANGE)

println()
@printf("Algorithm   : %s\n", uppercase(string(run.algorithm)))
@printf("Termination : %s (converged = %s), %d iterations\n", run.reason, run.converged, length(run.history) - 1)
@printf("Final RMS   : %.3f, %.1f s\n", run.rms, elapsed)
println("Outputs     : ", run.run_dir, "  (model.rho, data.pred)")
foreach(p -> println("  ", relpath(p, run.run_dir)), plots)

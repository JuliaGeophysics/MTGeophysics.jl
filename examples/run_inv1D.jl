# 1D MT inversion: every site of a data file on its own, Gauss-Newton or VFSA
# Author: @pankajkmishra
# The data file and one control file (examples/ctrl/1D) are the only inputs: MakeMesh1D lays out each site's
# layers from its skin depths and starts it from its median apparent resistivity. Results and plots go to
# run_YYYYmmdd_HHMMSS/ next to the data, one folder per site
# Generate the benchmark once: julia --project=. helpers/benchmarks_1D.jl
# Usage: julia --project=. examples/run_inv1D.jl [GN|VFSA]
#        julia --project=. examples/run_inv1D.jl data.dat InvCtrl

using MTGeophysics
using Printf

const CASE_DIR = joinpath(@__DIR__, "data", "1D-I")
const CTRL_DIR = joinpath(@__DIR__, "ctrl", "1D")

#---------- mesh ----------
const FIRST_LAYER_DIV = 5.0        # first layer = δ(f_max) / this
const VERTICAL_FACTOR = 1.1        # layer growth with depth
const DEPTH_MULT      = 4.0        # model depth = this × δ(f_min)

data_path, inv_path = length(ARGS) == 2 ? ARGS :
    (joinpath(CASE_DIR, "data.dat"), joinpath(CTRL_DIR, "InvCtrl." * uppercase(get(ARGS, 1, "GN"))))
all(isfile, (data_path, inv_path)) || error("missing $(filter(!isfile, [data_path, inv_path])); run julia --project=. helpers/benchmarks_1D.jl first")
true_model = joinpath(dirname(data_path), "model.true")

data = load_data2d(data_path)
meshes = MakeMesh1D(data; mode = ReadInvCtrl1D(inv_path).mode, first_layer_div = FIRST_LAYER_DIV,
                    vertical_factor = VERTICAL_FACTOR, depth_mult = DEPTH_MULT)
foreach(m -> @printf("%-10s %d layers to %.1f km, start %.1f Ω·m\n", m.site, length(m.thicknesses),
                     sum(m.thicknesses) / 1000, m.background), meshes)

elapsed = @elapsed run = Invert1D(data_path, inv_path, meshes)
plots = PlotInversion1D(run; true_model_path = isfile(true_model) ? true_model : nothing)

@printf("Algorithm : %s\n", uppercase(string(run.algorithm)))
foreach(r -> @printf("%-10s RMS %.3f  %s\n", r.site, r.rms, r.reason), run.sites)
@printf("RMS       : %.3f, %.1f s\n", run.rms, elapsed)
println("Outputs   : ", run.run_dir)
foreach(p -> println("  ", relpath(p, run.run_dir)), plots)

# 2D MT inversion of the COMEMI-III benchmark with Gauss–Newton, NLCG, or VFSA
# Author: @pankajkmishra
# Every algorithm takes the same inputs: observed data and a start model (its mesh)
# Generate the inputs once: julia --project=. helpers/benchmarks_2D.jl
# Usage: julia --project=. examples/run_inv2D.jl [gn|nlcg|vfsa]

using MTGeophysics
using Dates
using Printf

#---------- inputs ----------

const ALGORITHM  = Symbol(get(ARGS, 1, "gn"))            # :gn, :nlcg, or :vfsa
const CASE_DIR   = joinpath(@__DIR__, "0COMEMI2D-III")
const DATA_PATH  = joinpath(CASE_DIR, "Comemi2D3.obs")   # observed impedances and errors
const START_PATH = joinpath(CASE_DIR, "Comemi2D3.ini")   # start model, defines the mesh
const TRUE_PATH  = joinpath(CASE_DIR, "Comemi2D3.true")  # plots only; "" if unknown

const PLOT_DEPTH_KM   = 10.0      # depth limit of the *_core plots
const LOG10_RHO_RANGE = (0.0, 3.5)

#---------- algorithm settings ----------

# deterministic: shared objective and stopping controls, per-algorithm settings
const OPTIONS = Dict(
    :gn   => Inv2DOptions(max_iter = 20,  beta = 1.0, log_bounds = (0.0, 4.0), max_step = 0.5, target_rms = 1.0),
    :nlcg => Inv2DOptions(max_iter = 200, beta = 1.0, log_bounds = (0.0, 4.0), max_step = 0.5, target_rms = 1.0),
)
const CONFIG = Dict(
    :gn   => GaussNewton2DConfig(damping = 1e-2),
    :nlcg => NLCG2DConfig(restart = 30, precondition = true),
)

# stochastic
const VFSA_CONFIG = VFSA2DMTConfig(
    n_chains   = 2,
    n_ctrl     = 400,
    max_iter   = 3000,
    n_trials   = 1,
    log_bounds = (0.0, 4.0),
    step_scale = 0.11,
    seed       = 20260308,
    keep_models = true,
)

ALGORITHM in (:gn, :nlcg, :vfsa) || error("algorithm must be gn, nlcg, or vfsa, got $ALGORITHM")
isfile(DATA_PATH) && isfile(START_PATH) ||
    error("missing $DATA_PATH or $START_PATH; run julia --project=. helpers/benchmarks_2D.jl first")

#---------- run ----------

run_dir   = joinpath(@__DIR__, "Results", "$(ALGORITHM)2D_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
plots_dir = joinpath(run_dir, "plots")
mkpath(plots_dir)
println("Algorithm: ", ALGORITHM, "\nRun dir  : ", run_dir)

observed = load_data2d(DATA_PATH)
start    = load_model2d(START_PATH)
mesh     = build_mesh_from_model2d(start; frequencies = observed.frequencies, receiver_positions = observed.receivers)
@printf("Mesh     : %d x %d earth cells, %d frequencies %.3g-%.3g Hz, %d sites\n",
        length(mesh.y_cell_sizes), length(mesh.z_cell_sizes) - mesh.n_air_cells, length(mesh.frequencies),
        minimum(mesh.frequencies), maximum(mesh.frequencies), length(mesh.receiver_positions))

elapsed = @elapsed if ALGORITHM == :vfsa
    vfsa = VFSA2DMT(START_PATH, DATA_PATH; run_dir = run_dir, config = VFSA_CONFIG,
                    true_model_path = isfile(TRUE_PATH) ? TRUE_PATH : nothing)
    final_model = load_model2d(vfsa.best_chain.best_model_path).resistivity
    predicted   = data_from_response2d(vfsa.best_chain.best_response; z_xy_error = observed.z_xy_error,
                                       z_yx_error = observed.z_yx_error, site_names = observed.site_names)
    final_rms   = chi2_rms2d(observed, predicted).rms
    history     = nothing
else
    result = Invert2D(START_PATH, DATA_PATH; output_dir = run_dir,
                      algorithm = CONFIG[ALGORITHM], options = OPTIONS[ALGORITHM])
    final_model = result.resistivity
    predicted   = load_data2d(joinpath(run_dir, "data_$(inv2d_tag(result.algorithm)).dat"))
    final_rms   = result.fit.rms
    history     = result.history
end

#---------- plots ----------

label = Dict(:gn => "Gauss–Newton", :nlcg => "NLCG", :vfsa => "VFSA best")[ALGORITHM]
plot_path(name) = joinpath(plots_dir, name)
core = (show_padding = false, maximum_depth_km = PLOT_DEPTH_KM, resistivity_log10_range = LOG10_RHO_RANGE)
ρ_bg = exp(sum(log, start.resistivity[mesh.n_air_cells+1:end, :]) / length(start.resistivity[mesh.n_air_cells+1:end, :]))

plots = String[
    plot_mt2d_mesh(mesh; output_path = plot_path("mesh_full.png"), region = :full, background_resistivity = ρ_bg),
    plot_mt2d_mesh(mesh; output_path = plot_path("mesh_core.png"), region = :core, background_resistivity = ρ_bg),
    plot_mt2d_model(mesh, start.resistivity; output_path = plot_path("mstart_core.png"), title = "Start model", core...),
    plot_mt2d_model(mesh, final_model; output_path = plot_path("mfinal_core.png"),
                    title = @sprintf("%s model, RMS %.2f", label, final_rms), core...),
    plot_mt2d_model(mesh, final_model; output_path = plot_path("mfinal_full.png"), title = "$label model, full mesh",
                    resistivity_log10_range = LOG10_RHO_RANGE),
    plot_mt2d_data_fit(observed, predicted; output_path = plot_path("data_fit.png")),
    plot_mt2d_data_maps(data_to_response2d(observed); output_path = plot_path("data_obs_maps.png")),
    plot_mt2d_data_maps(data_to_response2d(predicted); output_path = plot_path("data_pred_maps.png")),
]
isfile(TRUE_PATH) && push!(plots, plot_mt2d_model(mesh, load_model2d(TRUE_PATH).resistivity;
    output_path = plot_path("mtrue_core.png"), title = "True model", core...))
history === nothing || push!(plots, plot_inv2d_convergence(history; output_path = plot_path("convergence.png"),
                                                           target_rms = OPTIONS[ALGORITHM].target_rms))

#---------- summary ----------

println()
history === nothing || @printf("Termination : %s (converged = %s), %d iterations\n",
                               result.reason, result.converged, length(history) - 1)
@printf("Final RMS   : %.3f, %.1f s\n", final_rms, elapsed)
println("Outputs     : ", run_dir)
foreach(p -> println("  ", relpath(p, run_dir)), plots)

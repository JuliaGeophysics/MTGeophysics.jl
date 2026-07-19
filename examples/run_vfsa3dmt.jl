# cascadia 3D VFSA MT inversion via MTGeophysics.jl
# needs the external ModEM forward solver and an MPI runtime (OpenMPI or MPICH)

using MTGeophysics
using Dates

const BASE_DIR = @__DIR__
const MODEM_EXECUTABLE = "/projappl/project_2011796/ModEM/bin/Mod3DMT_2026"

start_model = joinpath(BASE_DIR, "cascad_half_prior.ws")
observed_data = joinpath(BASE_DIR, "cascad_errfl5.dat")

# seed from the 1st CLI arg (e.g. SLURM job ID), fixed fallback for reproducibility
seed = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1911

#---------- configure the inversion ----------
cfg = VFSA3DMTConfig(
    nchains               = 1,          # independent chains
    nprocs                = 21,         # MPI ranks for ModEM
    mpirun_cmd            = "srun",     # MPI launcher
    modem_exe             = MODEM_EXECUTABLE,
    out_root              = joinpath(BASE_DIR, "run"),
    n_ctrl                = 900,        # RBF control points
    frac_update_controls  = 1.0,        # all controls perturbed per trial 
    log_bounds            = (-0.5, 4.0), # log10(Ω·m) model bounds
    step_scale            = 0.05,       # proposal step = step_scale × bound width 
    max_iter              = 5000,       # iteration cap
    n_trials              = 4,          # greedy best-of-4, deferred metropolis on the winner
    temp_kappa            = 0.05,       # T0 ~ 3x median uphill dE -> ~75% uphill acceptance at start
    cool_ratio            = 0.02,       # T_end = 1e-3: ~50% acceptance mid-run, quench by 5000
    target_rms            = 3.0,        # early-stop only; harmless if never reached
    seed                  = seed,
    pad_tol               = 0.2,        # core/padding detection tolerance
    padding_decay_length  = 10.0,       # padding blend (cells)
    keep_models           = false,      # discard trial models (this one was added because on CSC supercomputer limited storage)
    keep_dpred            = false,      # discard predicted-data files
    model_save_every      = 100,        # keep the winning trial model every 100 iters
    sigma_scale           = 2.0,        # RBF width in cells (1σ) at the top of the core
    sigma_scale_deep      = 2.0,        # RBF width at the bottom; log-depth interpolated
    trunc_sigmas          = 3.0,        # RBF support cutoff in σ
    ctrl_depth_power      = 0.25,       # mild shallow bias; p>=0.5 starves the deep (ceiling test v5)
    water_log10           = 0.5,        # freeze start-model cells below this (ocean at -0.523)
    bathymetry_file       = "",         # "" = derive mask from the start model
)

best_model, iter_log = VFSA3DMT(start_model; dobs_path=observed_data, cfg=cfg)

println("Best model: ", best_model)
println("Log: ", iter_log)

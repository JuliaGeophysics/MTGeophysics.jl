# MT3DINV 3D VFSA MT inversion via MTGeophysics.jl
# classic vfsa (1 trial per iteration); optional per-site 2x2 distortion
# correction (variable projection) via distortion_mode = :on, logged to
# distortion_best.txt in the run directory
# one chain per job: submit several jobs with different seeds for an ensemble
# needs the external ModEM forward solver and an MPI runtime (OpenMPI or MPICH)

using MTGeophysics
using Dates

const BASE_DIR = @__DIR__
# CPU-only build of the ModEM-GPU tree: znver5 (AVX-512 + FMA) vs the SSE2-only
# Mod3DMT_2026, so the forward solves run measurably faster on the Zen5 nodes
const MODEM_EXECUTABLE = "/projappl/project_2011796/ModEM-GPU/bin/Mod3DMT_CPU"

start_model = joinpath(BASE_DIR, "cascad_half_prior.ws")
observed_data = joinpath(BASE_DIR, "cascad_errfl5.dat")

# seed from the 1st CLI arg (e.g. SLURM job ID), fixed fallback for reproducibility
seed = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1911

#---------- configure the inversion ----------
cfg = VFSA3DMTConfig(
    nprocs                = 21,         # MPI ranks for ModEM
    mpirun_cmd            = "srun",     # MPI launcher
    modem_exe             = MODEM_EXECUTABLE,
    # pins BICG: this binary's compiled-in default is QMR (the maintained
    # ModEM tree defaults to BICG), and QMR is slower here
    fwd_ctrl              = joinpath(BASE_DIR, "F.dat"),
    out_root              = joinpath(BASE_DIR, "run"),
    n_ctrl                = 1600,        # RBF control points
    frac_update_controls  = 1.0,        # all controls perturbed per trial
    log_bounds            = (1.0, 5.0), # log10(Ω·m) model bounds
    step_scale            = 0.05,       # proposal step = step_scale × bound width
    max_iter              = 12000,       # iteration cap; also sets the cooling timescale
    n_trials              = 1,          # classic vfsa: one proposal, one metropolis test
    temp_kappa            = 0.03,       # hot start: T0 = 2x prior, broad uphill acceptance early
    cool_ratio            = 1.0e-4,     # fast cool: steep exp, effectively cold (~T0/100) by iter ~1250
    target_rms            = 1.0,        # early-stop; NOTE reachable at alpha=1 (ceiling ~3.4)
    seed                  = seed,
    pad_tol               = 0.2,        # core/padding detection tolerance
    padding_decay_length  = 10.0,       # padding blend (cells)
    z_core_cells          = 27,         # perturb top 27 z layers (~100 km); deeper layers are forward-solver padding
    padding_decay_length_z = 10.0,      # below-core blend e-fold (multiples of median core dz)
    keep_models           = false,      # discard trial models (limited storage on CSC)
    keep_dpred            = false,      # discard predicted-data files
    model_save_every      = 100,        # keep the winning trial model every 100 iters
    sigma_scale           = 2.0,        # RBF width in cells (1σ) at the top of the core
    sigma_scale_deep      = 2.0,        # RBF width at the bottom; equal = uniform widths
    trunc_sigmas          = 3.0,        # RBF support cutoff in σ
    ctrl_depth_power      = 0.0,       # shallow bias in control placement
    water_log10           = 0.3,        # no water mask (land survey)
    bathymetry_file       = "cascad_bathymetry.dat",         # bathymetry file
    distortion_mode       = :off,        # per-site 2x2 galvanic C each misfit eval; :off = plain
    distortion_damping    = 1.0,        # relative damping toward C=I; Inf = C≡I
)

best_model, iter_log = VFSA3DMT(start_model; dobs_path=observed_data, cfg=cfg)

println("Best model: ", best_model)
println("Log: ", iter_log)

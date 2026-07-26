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
    n_ctrl                = 900,        # RBF control points
    frac_update_controls  = 1.0,        # all controls perturbed per trial
    log_bounds            = (0.0, 4.0), # log10(Ω·m) model bounds, symmetric about the 2.0 prior; ocean is frozen and exempt from the clamp
    step_scale            = 0.15,       # caps a single control jump at step_scale × bound width (0.6 log10); the typical step is √T smaller
    max_iter              = 12000,       # iteration cap; also sets the cooling timescale
    n_trials              = 1,          # classic vfsa: one proposal, one metropolis test
    T0                    = 0.03,       # start temperature, on the scale of the typical uphill dE_rms2
    cool_ratio            = 1.0e-3,     # geometric cool to T0/10000 at max_iter: one decade per 3000 iters, T0/100 by iter ~6000
    target_rms            = 1.0,        # early-stop threshold; unreachable here, so the chain runs to max_iter
    seed                  = seed,
    pad_tol               = 0.2,        # core/padding detection tolerance
    core_expand_cells     = 2,          # grow the perturbable core by N cells per side; 2 brings all 109 sites inside (0 leaves 8 over padding)
    padding_decay_length  = 5.0,        # lateral padding blend e-fold, in core cells of true distance
    z_core_cells          = 30,         # perturb top 28 z layers (144 km); deeper layers are forward-solver padding
    keep_models           = false,      # discard trial models (limited storage on CSC)
    keep_dpred            = false,      # discard predicted-data files
    model_save_every      = 100,        # keep the winning trial model every 100 iters
    sigma_scale           = 2.0,        # RBF kernel 1σ in cells at the core top; physical footprint = σ × cell size per axis (2 cells ≈ 48 km laterally here)
    sigma_scale_deep      = 3.0,        # kernel 1σ at the core bottom, log-depth interpolated between the two; equal = uniform widths, raise to 3-4 to widen deep kernels laterally
    trunc_sigmas          = 3.0,        # kernel zeroed beyond this many σ; keeps the control-to-cell weight map sparse, 3σ drops only ~1% tail
    ctrl_depth_power      = 0.5,       # shallow bias in control placement
    water_log10           = 0.3,        # only used for the written file header when a bathymetry file is given
    bathymetry_file       = "cascad_bathymetry.dat",         # frozen ocean columns; relative path, so submit from this directory
    distortion_mode       = :off,        # per-site 2x2 galvanic C each misfit eval; :off = plain
    distortion_damping    = 1.0,        # pull of per-site C toward identity: 0 = free fit, larger = weaker correction, Inf = correction off
)

best_model, iter_log = VFSA3DMT(start_model; dobs_path=observed_data, cfg=cfg)

println("Best model: ", best_model)
println("Log: ", iter_log)

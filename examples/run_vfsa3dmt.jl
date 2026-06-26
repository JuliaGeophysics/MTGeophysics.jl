# example: 3D VFSA MT inversion on the bundled dataset
# needs the external ModEM forward solver and an MPI runtime (OpenMPI or MPICH)

using MTGeophysics
using Dates

const BASE_DIR = @__DIR__
const MODEM_EXECUTABLE = "/projappl/project_2005537/Mod3DMT_cc"

start_model = joinpath(BASE_DIR, "model.rho")
observed_data = joinpath(BASE_DIR, "data.dat")


#---------- configure the inversion ----------
cfg = VFSA3DMTConfig(
    nchains               = 1,          # independent chains
    nprocs                = 39,         # MPI ranks for ModEM
    mpirun_cmd            = "srun",     # MPI launcher
    modem_exe             = MODEM_EXECUTABLE,
    out_root              = joinpath(BASE_DIR, "run"),
    n_ctrl                = 900,        # RBF control points
    frac_update_controls  = 0.05,       # fraction of controls perturbed per trial
    log_bounds            = (0.0, 6.0), # log10(Ω·m) model bounds
    step_scale            = 1.0,        # proposal step = step_scale × bound width
    max_iter              = 3000,       # iteration cap
    n_trials              = 4,          # trials per iteration
    temp_kappa            = 1.0,        # dimensionless start temperature T0
    explore_frac          = 1.0,        # proposal reaches cool_ratio*T0 at this fraction of max_iter
    acc_freeze_frac       = 1.0,        # acceptance reaches cool_ratio*T0 at this fraction of max_iter
    cool_ratio            = 1e-3,       # cold endpoint T = cool_ratio*T0
    ak                    = nothing,    # set Float to override proposal decay rate
    ak_acc                = nothing,    # set Float to override acceptance decay rate
    target_rms            = 3.0,        # stop once best rms ≤ this
    seed                  = Dates.value(Dates.now()), # fix to an Int for reproducibility
    pad_tol               = 0.2,        # core/padding detection tolerance
    padding_decay_length  = 10.0,       # padding blend (cells)
    keep_models           = true,       # keep trial models (ignored when model_save_every>0)
    keep_dpred            = false,      # discard predicted-data files
    model_save_every      = 100,        # keep trial models every N iters (+ best)
    sigma_scale           = 2.0,        # RBF width in cells (1σ); larger = smoother
    trunc_sigmas          = 3.0,        # RBF support cutoff in σ
)

best_model, iter_log = VFSA3DMT(start_model; dobs_path=observed_data, cfg=cfg)


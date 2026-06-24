# example: 3D VFSA MT inversion on the bundled dataset
# needs the external ModEM forward solver and an MPI runtime (OpenMPI or MPICH)

using MTGeophysics
using Dates
using Statistics

#---------- user-configurable paths ----------
# resolve relative to this script
const START_MODEL_PATH = normpath(@__DIR__, "geoenergialoikka", "model.rho")
const OBSERVED_DATA_PATH = normpath(@__DIR__, "geoenergialoikka", "data.dat") 
const MODEM_EXECUTABLE = "/projappl/project_2005537/Mod3DMT_cc"

function _print_preflight_summary(start_model::AbstractString, observed_data::AbstractString, cfg::VFSA3DMTConfig)
    m = load_ws3d_model(start_model)
    ix = core_indices(m.cx; tol = cfg.pad_tol)
    iy = core_indices(m.cy; tol = cfg.pad_tol)
    println("3D VFSA preflight summary")
    println("  Model: $start_model")
    println("  Data:  $observed_data")
    println("  Mesh cells: $(m.nx) x $(m.ny) x $(m.nz)")
    println("  Detected core x-range: $ix ($(length(ix)) cells)")
    println("  Detected core y-range: $iy ($(length(iy)) cells)")
    println("  Nominal padding estimate: $(m.npad)")
    println("  Core cell widths (x median / y median): $(median(m.dx[ix])) / $(median(m.dy[iy]))")
end

#---------- resolve input paths ----------
start_model = START_MODEL_PATH
observed_data = OBSERVED_DATA_PATH

# verify files exist
for (label, path) in [("Starting model", start_model), ("Observed data", observed_data)]
    if !isfile(path)
        error("$label not found: $path")
    end
end

#---------- configure the inversion ----------
cfg = VFSA3DMTConfig(
    nchains               = 1,          # independent chains
    nprocs                = 39,         # MPI ranks for ModEM
    mpirun_cmd            = "srun",     # MPI launcher
    modem_exe             = MODEM_EXECUTABLE,
    out_root              = "runs",     # run dir base -> runs_<seed>
    n_ctrl                = 900,        # RBF control points
    frac_update_controls  = 0.05,       # fraction of controls perturbed per trial
    log_bounds            = (0.0, 6.0), # log10(Ω·m) model bounds
    step_scale            = 1.0,        # proposal step = step_scale × bound width
    max_iter              = 3000,       # iteration cap
    n_trials              = 4,          # trials per iteration
    # decoupled temperatures: T0 = temp_kappa*rms_initial^2, cooling as T0*exp(-sqrt(ak)*(k-1)).
    # T_prop drives proposals, T_acc the metropolis test; each reaches cool_ratio*T0 at its
    # own fraction of max_iter. acc_freeze_frac >= explore_frac makes acceptance cool slower.
    temp_kappa            = 1.0,        # T0 = temp_kappa * rms_initial^2
    explore_frac          = 0.15,       # proposal reaches cool_ratio*T0 at this fraction of max_iter
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

#---------- optional preflight (inspect mesh without ModEM) ----------
if get(ENV, "MTG_PREFLIGHT_ONLY", "0") == "1"
    _print_preflight_summary(start_model, observed_data, cfg)
    exit(0)
end

best_model, iter_log = VFSA3DMT(start_model; dobs_path=observed_data, cfg=cfg)


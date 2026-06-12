# Example: 3D VFSA MT Inversion
#
# This script demonstrates how to set up and run a 3D VFSA inversion
# with MTGeophysics.jl using an example dataset bundled with the repository.
#
# IMPORTANT: This example requires the external ModEM forward solver
# and an MPI runtime (e.g. OpenMPI or MPICH).
#
# The model and data paths are defined below relative to this script.

using MTGeophysics

# --------------------------------------------------------------------------
# User-configurable paths
# --------------------------------------------------------------------------
# These resolve relative to the location of this script.
const START_MODEL_PATH = normpath(@__DIR__, "geoenergialoikka", "model.rho")
const OBSERVED_DATA_PATH = normpath(@__DIR__, "geoenergialoikka", "data.dat")
const MODEM_EXECUTABLE = "/usr/local/bin/Mod3DMT_2025"

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

# --------------------------------------------------------------------------
# Resolve input paths
# --------------------------------------------------------------------------
start_model = START_MODEL_PATH
observed_data = OBSERVED_DATA_PATH

# Verify files exist
for (label, path) in [("Starting model", start_model), ("Observed data", observed_data)]
    if !isfile(path)
        error("$label not found: $path")
    end
end

# --------------------------------------------------------------------------
# Configure the 3D VFSA inversion
# --------------------------------------------------------------------------
cfg = VFSA3DMTConfig(
    nchains               = 1,          # number of independent Markov chains
    nprocs                = 21,         # MPI processes for ModEM forward calls
    mpirun_cmd            = "mpirun",
    modem_exe             = MODEM_EXECUTABLE,
    out_root              = "runs",     # scratch dir, created next to the model file
    n_ctrl                = 900,        # RBF control points in the core
    log_bounds            = (0.0, 5.0), # log10(Ω·m) bounds
    step_scale            = 0.05,       # VFSA proposal step size
    max_iter              = 3000,       # total VFSA iterations (use 10-50 for smoke test)
    n_trials              = 4,          # trial proposals per iteration
    T0_prop               = 1.0,        # initial proposal temperature
    Tf_prop               = 1e-3,       # final proposal temperature
    T0_acc                = 1.0,        # initial acceptance temperature
    Tf_acc                = 1e-3,       # final acceptance temperature
    seed                  = 1911,       # random seed for reproducibility
    padding_decay_length  = 10.0,       # horizontal padding blend (in cell widths)
    keep_models           = true,       # keep all trial model files
    keep_dpred            = false,      # discard predicted data files to save space
)

# --------------------------------------------------------------------------
# Optional preflight mode for debugging input meshes without ModEM
# --------------------------------------------------------------------------
if get(ENV, "MTG_PREFLIGHT_ONLY", "0") == "1"
    _print_preflight_summary(start_model, observed_data, cfg)
    exit(0)
end

# --------------------------------------------------------------------------
# Check that the external ModEM solver is available
# --------------------------------------------------------------------------
modem_exe = cfg.modem_exe

modem_found = try
    success(`which $modem_exe`)
catch
    try
        success(`where $modem_exe`)
    catch
        false
    end
end

if !modem_found
    printstyled("""

    ╔══════════════════════════════════════════════════════════════════╗
    ║  ModEM forward solver ("$modem_exe") was not found.            ║
    ║                                                                ║
    ║  The 3D VFSA inversion requires the ModEM 3D MT code and an   ║
    ║  MPI runtime to compute forward responses.                    ║
    ║                                                                ║
    ║  Installation:                                                 ║
    ║    1. Download ModEM from:                                     ║
    ║       https://github.com/magnetotellurics/ModEM               ║
    ║    2. Build with MPI support (see ModEM documentation).        ║
    ║    3. Install an MPI runtime (OpenMPI, MPICH, MS-MPI, etc.)   ║
    ║                                                                ║
    ║  Current executable path:                                      ║
    ║    $modem_exe
    ╚══════════════════════════════════════════════════════════════════╝

    """; color=:yellow, bold=true)
    exit(1)
end

# --------------------------------------------------------------------------
# Run the inversion
# --------------------------------------------------------------------------
println("Starting 3D VFSA MT inversion...")
println("  Model:  $start_model")
println("  Data:   $observed_data")
println("  Config: $(cfg.nchains) chain(s), $(cfg.max_iter) iterations, $(cfg.n_trials) trials/iter")
println()

best_model, iter_log = VFSA3DMT(start_model; dobs_path=observed_data, cfg=cfg)

println()
println("Inversion complete.")
println("  Best model : $best_model")
println("  Iteration log: $iter_log")

# --------------------------------------------------------------------------
# (Optional) Run ensemble statistics if multiple chains were used
# --------------------------------------------------------------------------
if cfg.nchains > 1
    println()
    println("Computing ensemble statistics...")
    mean_path, median_path, std_path = AnalyseEnsemble3D(dirname(start_model))
    println("  Mean   model: $mean_path")
    println("  Median model: $median_path")
    println("  Std    model: $std_path")
end

# --------------------------------------------------------------------------
# (Optional) Visualize the best model using the existing 3D viewer
# --------------------------------------------------------------------------
# If GLMakie is available, you can load and view the result:
#
#   m = load_ws3d_model(best_model)
#   # Convert to linear for the ModEMModel viewer:
#   m_modem = load_model_modem(best_model)
#   # Then use plot_model_XYZ.jl or gl_modem_viewer() from examples/

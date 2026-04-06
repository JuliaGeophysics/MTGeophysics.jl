# Example: 3D VFSA MT Inversion using Cascadia data
#
# This script demonstrates how to set up and run a 3D VFSA inversion
# with MTGeophysics.jl using the included Cascadia dataset.
#
# IMPORTANT: This example requires the external ModEM forward solver
# (Mod3DMT / Mod3DMT_2025) to be installed and accessible on your PATH,
# along with an MPI runtime (e.g. OpenMPI or MPICH).
#
# If ModEM is not installed, this script will print setup instructions
# and exit. To install ModEM, please read the documentation at:
#   https://sites.google.com/site/modularem/
#
# Usage:
#   julia --project=. examples/run_vfsa3dmt.jl
#
# The Cascadia example data is from:
#   Patro & Egbert (2008), Geophys. Res. Lett., 35, L20311.

using MTGeophysics

# --------------------------------------------------------------------------
# Check that the external ModEM solver is available
# --------------------------------------------------------------------------
modem_exe = "Mod3DMT"  # change to "Mod3DMT_2025" if using the 2025 build

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
    ║  ModEM forward solver ("$modem_exe") not found on PATH.        ║
    ║                                                                ║
    ║  The 3D VFSA inversion requires the ModEM 3D MT code and an   ║
    ║  MPI runtime to compute forward responses.                    ║
    ║                                                                ║
    ║  Installation:                                                 ║
    ║    1. Download ModEM from:                                     ║
    ║       https://github.com/dong-hao/ModEM-GPU                   ║
    ║    2. Build with MPI support (see ModEM documentation).        ║
    ║    3. Ensure the executable is on your PATH.                   ║
    ║    4. Install an MPI runtime (OpenMPI, MPICH, MS-MPI, etc.)   ║
    ║                                                                ║
    ║  Once installed, re-run this script.                           ║
    ╚══════════════════════════════════════════════════════════════════╝

    """; color=:yellow, bold=true)
    exit(1)
end

# --------------------------------------------------------------------------
# Paths to Cascadia example data
# --------------------------------------------------------------------------
cascadia_dir = normpath(@__DIR__, "Cascadia")

start_model  = joinpath(cascadia_dir, "cascad_half_prior.ws")
observed_data = joinpath(cascadia_dir, "cascad_errfl5.dat")

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
    mpirun_cmd            = "mpirun",   # MPI launcher command
    modem_exe             = modem_exe,  # ModEM executable name
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
    mean_path, median_path, std_path = AnalyseEnsemble3D(cfg.out_root)
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
#   # Then use plot_model_3D.jl or gl_modem_viewer() from examples/

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
#   julia --project=. examples/run_vfsa3dmt.jl path/to/model.rho path/to/data.dat
#   julia --project=. examples/run_vfsa3dmt.jl --check path/to/model.rho path/to/data.dat
#
# The Cascadia example data is from:
#   Patro & Egbert (2008), Geophys. Res. Lett., 35, L20311.

using MTGeophysics

function _parse_cli(args::Vector{String})
    check_only = any(==("--check"), args)
    positional = filter(!=("--check"), args)
    length(positional) in (0, 2) || error("Usage: julia --project=. examples/run_vfsa3dmt.jl [--check] [model_path data_path]")
    return check_only, positional
end

function _resolve_paths(args::Vector{String})
    check_only, positional = _parse_cli(args)
    if length(positional) == 2
        start_model = abspath(positional[1])
        observed_data = abspath(positional[2])
    else
        cascadia_dir = normpath(@__DIR__, "Cascadia")
        start_model = get(ENV, "MTG_MODEL_PATH", joinpath(cascadia_dir, "cascad_half_prior.ws"))
        observed_data = get(ENV, "MTG_DATA_PATH", joinpath(cascadia_dir, "cascad_errfl5.dat"))
    end
    return check_only, start_model, observed_data
end

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
check_only, start_model, observed_data = _resolve_paths(ARGS)

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
    mpirun_cmd            = get(ENV, "MTG_MPIRUN", "mpirun"),
    modem_exe             = get(ENV, "MTG_MODEM_EXE", "Mod3DMT"),
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
if check_only || get(ENV, "MTG_PREFLIGHT_ONLY", "0") == "1"
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
    ║  ModEM forward solver ("$modem_exe") not found on PATH.        ║
    ║                                                                ║
    ║  The 3D VFSA inversion requires the ModEM 3D MT code and an   ║
    ║  MPI runtime to compute forward responses.                    ║
    ║                                                                ║
    ║  To debug inputs without ModEM, run:                           ║
    ║    julia --project=. examples/run_vfsa3dmt.jl --check         ║
    ║    julia --project=. examples/run_vfsa3dmt.jl --check         ║
    ║      /path/to/model.rho /path/to/data.dat                     ║
    ║                                                                ║
    ║  Installation:                                                 ║
    ║    1. Download ModEM from:                                     ║
    ║       https://github.com/dong-hao/ModEM-GPU                   ║
    ║    2. Build with MPI support (see ModEM documentation).        ║
    ║    3. Ensure the executable is on your PATH.                   ║
    ║    4. Install an MPI runtime (OpenMPI, MPICH, MS-MPI, etc.)   ║
    ║                                                                ║
    ║  Once installed, re-run this script without --check.           ║
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

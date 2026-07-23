# 3D VFSA Inversion

Run a 3D magnetotelluric inversion using Very Fast Simulated Annealing (VFSA)
with Gaussian-RBF parameterization and the external
[ModEM](https://github.com/magnetotellurics/ModEM) forward solver. Each run is
a single chain; for an ensemble, submit several jobs with different seeds and
output directories.

![VFSA 3D benchmark — ensemble mean and std vs deterministic ModEM inversion](assets/VFSA3DBenchmark.png)

!!! note "External dependency"
    The 3D VFSA inversion requires **ModEM** (Mod3DMT / Mod3DMT\_2025) compiled
    with MPI support and accessible on your `PATH`, along with an MPI runtime
    (OpenMPI, MPICH, or MS-MPI).

## Running the example

```bash
julia --project=. examples/run_vfsa3dmt.jl
```

The script checks for the ModEM executable before starting. If it is not found
it prints setup instructions and exits.

## From Julia

```julia
using MTGeophysics

cfg = VFSA3DMTConfig(
    nprocs    = 21,
    n_ctrl    = 900,
    log_bounds = (0.0, 5.0),
    max_iter  = 3000,
    n_trials  = 1,
    seed      = 1911,
    keep_models = true,
)

best_model, iter_log = VFSA3DMT(
    "examples/Cascadia/cascad_half_prior.ws";
    dobs_path = "examples/Cascadia/cascad_errfl5.dat",
    cfg = cfg,
)
```

## Configuration

| Parameter | Example | Description |
|:----------|:--------|:------------|
| `nprocs` | 21 | MPI processes for ModEM forward calls |
| `mpirun_cmd` | `"mpirun"` | MPI launcher command |
| `modem_exe` | `"Mod3DMT_2025"` | ModEM executable name |
| `n_ctrl` | 900 | Gaussian-RBF control points in the core |
| `log_bounds` | (0, 5) | log₁₀(Ω·m) search bounds |
| `step_scale` | 0.05 | VFSA proposal step size |
| `max_iter` | 3000 | Iteration cap |
| `n_trials` | 1 | Trial perturbations per iteration (1 = classic VFSA) |
| `temp_kappa` | 0.05 | Initial temperature (proposal and acceptance) |
| `cool_ratio` | 0.02 | Final temperature as a fraction of `temp_kappa` |
| `target_rms` | 3.0 | Early-stop RMS |
| `seed` | 1911 | Random seed for reproducibility |
| `pad_tol` | 0.2 | Tolerance for core/padding detection |
| `padding_decay_length` | 10.0 | Horizontal padding blend (cell widths) |
| `keep_models` | true | Keep all trial model files |
| `keep_dpred` | false | Keep predicted data files |
| `sigma_scale` | 2.0 | RBF width scaling factor |
| `trunc_sigmas` | 3.0 | RBF truncation distance (in σ) |

## WS3D model format

The 3D VFSA engine uses the **WS3D model format** (Siripunvaraporn et al.),
which stores resistivity in log₁₀(Ω·m). Conversion between WS3D and ModEM
linear formats is handled automatically.

```julia
# Read / write WS3D models directly
m = load_ws3d_model("model.rho")
write_ws3d_model("output.rho", m.dx, m.dy, m.dz, m.A, m.origin)
```

## Ensemble analysis

Run several jobs with different seeds, then point `AnalyseEnsemble3D` at the
parent directory of the run directories (searched recursively for
`best_model.rho`):

```julia
mean_path, median_path, std_path = AnalyseEnsemble3D("runs")
```

This writes `model.mean`, `model.median`, and `model.std` in WS3D format.

The lower-level helper `core_statistics(cores)` computes element-wise mean,
median, and standard deviation over an array of 3D resistivity cubes.

## Output

The result directory contains:

```text
run_<timestamp>/
├── 0vfsa3DMT.log              # per-iteration summary
├── 0vfsa3DMT_detailed.log     # per-trial detail
├── best_model.rho             # best model of the run
├── model_001_01.rho           # trial models (if keep_models=true)
├── dpred_001_01.dat           # predicted data (if keep_dpred=true)
├── distortion_best.txt        # per-site distortion (if distortion_mode=:on)
└── ...
```

`AnalyseEnsemble3D` writes `model.mean`, `model.median`, and `model.std` into
the parent directory it is given.

## Example data

The Cascadia 3D example is not bundled. Download it from
[ModEM-Examples](https://github.com/magnetotellurics/ModEM-Examples/tree/main/Magnetotelluric/3D_MT/Cascadia)
and place it in `examples/Cascadia/`.

# MTGeophysics.jl — Examples Guide

This guide walks through every step from generating benchmark files to running
forward modelling, plotting, inversion, and convergence GIF generation.

All commands assume you are in the `MTGeophysics.jl/` project root.

---

## 1. Prerequisites

Make sure the package compiles cleanly:

```bash
julia --project=. -e 'using MTGeophysics; println("OK")'
```

---

## 2. Generate benchmark files

Benchmark models and synthetic observed data live under `Examples/` in
directories prefixed with `0` so they sort before result directories.

### 1-D benchmarks (layered half-space)

```bash
julia --project=. Helpers/benchmarks_1d.jl
```

Writes into `Examples/0Layered1D/`:

| File | Description |
|------|-------------|
| `Layered1D.true` | True 3-layer resistivity model |
| `Layered1D.ref` | Reference data template (frequencies, receivers) |
| `Layered1D.obs` | Synthetic observed data with 5 % Gaussian noise |

### 2-D benchmarks (COMMEMI models)

```bash
julia --project=. Helpers/benchmarks_2d.jl
```

Writes into three directories:

| Directory | COMMEMI case |
|-----------|--------------|
| `Examples/0COMEMI2D-I/` | COMMEMI 2D-1 |
| `Examples/0COMEMI2D-II/` | COMMEMI 2D-2 |
| `Examples/0COMEMI2D-III/` | COMMEMI 2D-3 |

Each directory contains:

| File | Description |
|------|-------------|
| `*.true` | True resistivity model |
| `*.ini` | Uniform half-space starting model (for inversion) |
| `*.ref` | Reference data template |
| `*.obs` | Synthetic observed data with noise |

---

## 3. Forward modelling

### 1-D forward response

```bash
julia --project=. Examples/response_1d.jl                 # uses default Layered1D model
julia --project=. Examples/response_1d.jl model.true data.ref   # custom paths
```

Computes the 1-D forward response, writes a `.pred` data file, and plots
observed-vs-predicted curves.

### 2-D forward response

```bash
julia --project=. Examples/response_2d.jl                 # uses default COMEMI 2D-I model
julia --project=. Examples/response_2d.jl model.true data.ref   # custom paths
```

Computes the 2-D forward response (TE + TM modes), writes a `.pred` file,
and generates data-map and site-curve plots.

---

## 4. Plotting

### Plot a 1-D model

```bash
julia --project=. Examples/plot_model_1d.jl \
    Examples/0Layered1D/Layered1D.true \
    Examples/0Layered1D/Layered1D.obs
```

Optional: `--maximum-depth-m 5000` to limit the plotted depth.

### Plot 1-D data

```bash
julia --project=. Examples/plot_data_1d.jl \
    Examples/0Layered1D/Layered1D.obs
```

To overlay predicted data on the observed data:

```bash
julia --project=. Examples/plot_data_1d.jl \
    Examples/0Layered1D/Layered1D.obs \
    path/to/predicted.pred
```

### Plot a 2-D model

```bash
julia --project=. Examples/plot_model_2d.jl \
    Examples/0COMEMI2D-I/Comemi2D1.true
```

Options:

| Flag | Default | Description |
|------|---------|-------------|
| `--show-grid true` | `false` | Draw cell boundaries |
| `--show-padding false` | `true` | Hide padding cells |
| `--maximum-depth-km 50` | `Inf` | Limit plotted depth |
| `--resistivity-log10-range 0,4` | `0.0,4.0` | Colour-scale bounds |

### Plot 2-D data

```bash
julia --project=. Examples/plot_data_2d.jl \
    Examples/0COMEMI2D-I/Comemi2D1.obs
```

Produces two plots: `DataMaps2D.png` (pseudo-section maps) and
`SiteCurves2D.png` (per-station apparent resistivity and phase curves).

---

## 5. 2-D VFSA inversion

### Option A — script with hardcoded parameters

Edit `Examples/run_vfsa2dmt.jl` to set the model paths and inversion
parameters, then run:

```bash
julia --project=. Examples/run_vfsa2dmt.jl
```

### Option B — CLI with flags

The `VFSA2DMT` module exposes a CLI entry point:

```bash
julia --project=. -e 'using MTGeophysics; main_vfsa2dmt()' -- \
    Examples/0COMEMI2D-I/Comemi2D1.ini \
    Examples/0COMEMI2D-I/Comemi2D1.obs \
    --n-chains 3 --n-ctrl 400 --max-iter 100 \
    --log-bounds 0,4 --seed 20260308
```

Available CLI flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--n-chains N` | 4 | Number of SA chains |
| `--n-ctrl N` | 200 | Number of RBF control points |
| `--max-iter N` | 500 | Maximum iterations per chain |
| `--n-trials N` | 1 | Trials per iteration |
| `--log-bounds lo,hi` | `0,4` | log₁₀(ρ) search bounds |
| `--perturb-depth D` | `Inf` | Maximum perturbation depth (m) |
| `--seed N` | 42 | RNG seed |
| `--output-root DIR` | project root | Base directory for output |

### Inversion output structure

A run creates a timestamped directory, for example
`run_VFSA2DMT_20260309_235858/`, containing:

```
run_VFSA2DMT_20260309_235858/
├── Summary.md              # Run summary with fit statistics
├── chain_01/               # Per-chain logs and snapshots
│   ├── 0vfsa2DMT.log       # Iteration-level log
│   ├── 0vfsa2DMT.details   # Trial-level log
│   ├── best_iter_00010.rho # Best-so-far model at iter 10
│   ├── best_iter_00020.rho
│   └── itr_00001_trial_01.rho  # (if keep_models=true)
├── chain_02/
├── chain_03/
├── snapshots/              # Cross-chain averaged models
│   ├── avg_iter_00010.rho
│   ├── avg_iter_00020.rho
│   └── ...
├── data.obs                # Copy of observed data
├── model.start             # Copy of starting model
├── model.true              # True model (if auto-detected)
├── model.mean              # Ensemble mean model
├── model.median            # Ensemble median model
├── model.std               # Ensemble standard deviation
├── model.c1best            # Best model from chain 1
├── data.c1best             # Best predicted data from chain 1
├── plot_model_mean.png     # Model plots
├── plot_convergence.png    # χ² and RMS convergence curves
└── convergence.gif         # Animated convergence (if GIF generated)
```

---

## 6. Post-run analysis

### Recompute ensemble statistics

If you want to recompute the mean/median models from saved chain-best models:

```bash
julia --project=. Helpers/run_statistics_2d.jl <run_dir>
```

Optionally pass the observed data path as a second argument if it cannot be
auto-detected.

---

## 7. Generate a convergence GIF

After an inversion run that wrote snapshot models, create an animated GIF
showing how the cross-chain mean model evolves over iterations:

```bash
julia --project=. Helpers/make_gif_2d.jl <run_dir>
```

For example:

```bash
julia --project=. Helpers/make_gif_2d.jl Examples/run_VFSA2DMT_20260309_235858
```

This reads `snapshots/avg_iter_*.rho` and writes `convergence.gif` in the run
directory.

Options:

| Flag | Default | Description |
|------|---------|-------------|
| `--fps N` | `4` | Frames per second |
| `--depth_km D` | `Inf` | Maximum plotted depth (km) |
| `--rho_range lo,hi` | `0.0,4.0` | log₁₀(ρ) colour scale |

You can also specify a custom output path:

```bash
julia --project=. Helpers/make_gif_2d.jl <run_dir> my_animation.gif --fps 8
```

---

## 8. Complete workflow example

A full end-to-end run using the COMMEMI 2D-1 benchmark:

```bash
# 1. Generate benchmark files
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl

# 2. Plot the true model and observed data
julia --project=. Examples/plot_model_2d.jl Examples/0COMEMI2D-I/Comemi2D1.true
julia --project=. Examples/plot_data_2d.jl  Examples/0COMEMI2D-I/Comemi2D1.obs

# 3. Run the 2D VFSA inversion (3 chains, 100 iterations)
julia --project=. Examples/run_vfsa2dmt.jl

# 4. Generate the convergence GIF from the run directory
julia --project=. Helpers/make_gif_2d.jl Examples/run_VFSA2DMT_<timestamp>

# 5. Recompute ensemble statistics (optional)
julia --project=. Helpers/run_statistics_2d.jl Examples/run_VFSA2DMT_<timestamp>
```

Replace `<timestamp>` with the actual directory name created by the inversion.

---

## 9. Tests

Run the full test suite to verify everything works:

```bash
julia --project=. tests/runtests.jl
```

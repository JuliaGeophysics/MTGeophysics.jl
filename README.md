<h1 align="center">MTGeophysics.jl</h1>

<p align="center"><em>Magnetotelluric modelling, inversion, and visualization workflows in Julia.</em></p>

<p align="center">
	<a href="https://juliageophysics.github.io/MTGeophysics.jl/dev/">Documentation</a> |
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/issues">Issues</a> |
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/releases">Releases</a>
</p>

---


## Features

- 1D and 2D MT forward solvers (analytical and finite-difference)
- 2D/3D VFSA inversion with ensemble uncertainty quantification
- ModEM 3D model and data I/O
- Interactive 3D slice viewers (GLMakie) with GIS overlays and coordinate reprojection
- Shapefile export for GIS integration

## Installation

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick start: 2D benchmark + VFSA inversion

> Note: This is a quick demonstration workflow. For production-quality inversion, tune key settings such as the number of VFSA iterations, number of chains, cooling schedule, regularization choices, and uncertainty/ensemble controls for your survey and model size.

**1. Generate the COMMEMI 2D benchmarks.** Writes the true model, halfspace starting model, reference data template, and noisy observed data into `examples/0COMEMI2D-I/`, `-II/`, and `-III/`.

```
julia --project=. helpers/benchmarks_2d.jl
```

**2. (Optional) Inspect the 2D forward response** of the benchmark model.

```
julia --project=. examples/response_2d.jl
```

**3. Run the 2D VFSA inversion** (multi-chain ensemble, RBF parameterization) on the COMMEMI-I benchmark. Results are written to a timestamped `examples/run_VFSA2DMT_<timestamp>/` directory containing per-chain logs, best models, and ensemble mean/median/std.

```
julia --project=. examples/run_vfsa2dmt.jl
```

**4. Compute ensemble statistics** (posterior mean/median/std and misfit summaries).

```
julia --project=. helpers/run_statistics_2d.jl examples/run_VFSA2DMT_<timestamp>
```

**5. Create a convergence GIF** to visualize model evolution across iterations (requires `keep_models = true`).

```
julia --project=. helpers/make_gif_2d.jl examples/run_VFSA2DMT_<timestamp>
```

From Julia, the inversion can also be driven directly:

```julia
using MTGeophysics

result = VFSA2DMT(
    VFSA2DMTParams(
        script_path      = @__FILE__,
        start_model_path = "examples/0COMEMI2D-I/Comemi2D1.ini",
        data_path        = "examples/0COMEMI2D-I/Comemi2D1.obs",
        config = VFSA2DMTConfig(
            n_chains    = 2,
            n_ctrl      = 400,
            max_iter    = 3000,
            n_trials    = 4,
            log_bounds  = (0.0, 4.0),
            seed        = 20260308,
            keep_models = true,
        ),
    ),
)
```

See the [documentation](https://juliageophysics.github.io/MTGeophysics.jl/dev/) for the full 2D/3D workflow, configuration options, and visualization examples.

# Research using this code 

- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., and Sen, M. K.: Building uncertainty-aware subsurface models with 3D magnetotelluric inversion, EGU General Assembly 2026, Vienna, Austria, 3–8 May 2026, EGU26-4367, https://doi.org/10.5194/egusphere-egu26-4367, 2026. 

- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026. 

- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., Xiao, L., and Sen, M. K.: Stochastic model exploration in three-dimensional inversion of magnetotelluric data, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026. 

- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, JuliaCon 2026. 

- Patzer, C., Mishra, P. K., and Kamm, J.: Studying the Wiborg Rapakivi Batholith in SE Finland, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026. 




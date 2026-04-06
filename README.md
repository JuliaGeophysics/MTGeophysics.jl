<h1 align="center">MTGeophysics.jl</h1>

<p align="center"><em>Magnetotelluric modelling, inversion, and visualization workflows in Julia.</em></p>

<p align="center">
	<a href="https://juliageophysics.github.io/MTGeophysics.jl/">
		<img src="https://img.shields.io/badge/Documentation-Visit%20Site-0A7EA4?style=for-the-badge&logo=readthedocs&logoColor=white" alt="Documentation" />
	</a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/Documenter.yml">
		<img src="https://img.shields.io/github/actions/workflow/status/JuliaGeophysics/MTGeophysics.jl/Documenter.yml?branch=main&label=Docs%20Build&style=for-the-badge" alt="Docs Build" />
	</a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/releases">
		<img src="https://img.shields.io/github/v/release/JuliaGeophysics/MTGeophysics.jl?style=for-the-badge&label=Release" alt="Release" />
	</a>
</p>

<p align="center">
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/main/LICENSE.md">
		<img src="https://img.shields.io/github/license/JuliaGeophysics/MTGeophysics.jl?style=flat-square" alt="License" />
	</a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/stargazers">
		<img src="https://img.shields.io/github/stars/JuliaGeophysics/MTGeophysics.jl?style=flat-square" alt="GitHub stars" />
	</a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/issues">
		<img src="https://img.shields.io/github/issues/JuliaGeophysics/MTGeophysics.jl?style=flat-square" alt="Issues" />
	</a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/pulls">
		<img src="https://img.shields.io/github/issues-pr/JuliaGeophysics/MTGeophysics.jl?style=flat-square" alt="Pull requests" />
	</a>
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

## Repository structure

```
MTGeophysics.jl/
├── src/                        # Package source
│   ├── MTGeophysics.jl         # Main module (includes & exports)
│   ├── MTGeophysics1D.jl       # 1D forward solvers
│   ├── MTGeophysics2D.jl       # 2D forward solvers
│   ├── VFSA2DMT.jl             # 2D VFSA inversion engine
│   ├── VFSA3DMT.jl             # 3D VFSA inversion engine
│   ├── Model.jl                # ModEM 3D model I/O
│   ├── Data.jl                 # ModEM 3D data I/O
│   ├── WS3DModel.jl            # WS3D model format I/O
│   ├── CoreUtils3D.jl          # Core/padding detection utilities
│   ├── PlotModel.jl            # Visualization helpers (GLMakie)
│   └── Chi2RMS.jl              # Misfit statistics
├── examples/                   # Runnable example scripts
│   ├── response_1d.jl          # 1D forward response
│   ├── response_2d.jl          # 2D forward response
│   ├── run_vfsa2dmt.jl         # 2D VFSA inversion
│   ├── run_vfsa3dmt.jl         # 3D VFSA inversion (requires ModEM)
│   ├── plot_model_3D.jl        # Interactive 3D slice viewer
│   ├── plot_XY_slices.jl       # XY depth-slice viewer
│   ├── plot_XZ_slices.jl       # XZ cross-section viewer
│   ├── plot_YZ_slices.jl       # YZ cross-section viewer
│   ├── replace_slice_resistivity_scope.jl  # Bulk resistivity editing
│   ├── draw_and_replace_in_model.jl        # Polygon zone editing
│   └── Cascadia/               # Example data (not tracked)
├── Helpers/                    # Benchmark generation & post-processing
│   ├── benchmarks_1d.jl        # Generate 1D benchmarks
│   ├── benchmarks_2d.jl        # Generate 2D benchmarks
│   ├── run_statistics_2d.jl    # 2D ensemble statistics
│   ├── run_statistics_3D.jl    # 3D ensemble statistics
│   └── make_gif_2d.jl          # 2D convergence animation
├── docs/                       # Documenter.jl site
├── test/                       # Unit tests
├── paper/                      # JOSS paper
├── Project.toml
└── README.md
```

## Quick example: 2D VFSA inversion

> Note: This is a quick demonstration workflow. For production-quality inversion, tune key settings such as the number of VFSA iterations, number of chains, cooling schedule, regularization choices, and uncertainty/ensemble controls for your survey and model size.

1. Build a synthetic 2D benchmark model and generate inputs.

```
julia --project=. Helpers/benchmarks_2d.jl
```

2. Run the VFSA inversion workflow for the benchmark.

```
julia --project=. examples/run_vfsa2dmt.jl
```

3. Compute ensemble statistics for posterior models and misfit summaries.

```
julia --project=. Helpers/run_statistics_2d.jl
```

4. Create a GIF to visualize inversion evolution through iterations.

```
julia --project=. Helpers/make_gif_2d.jl
```

# Research using this code 

- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., and Sen, M. K.: Building uncertainty-aware subsurface models with 3D magnetotelluric inversion, EGU General Assembly 2026, Vienna, Austria, 3–8 May 2026, EGU26-4367, https://doi.org/10.5194/egusphere-egu26-4367, 2026. 




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




<h1 align="center">MTGeophysics.jl</h1>

<p align="center"><em>A software repository for magnetotelluric research and applications.</em></p>

<p align="center">
	<a href="https://juliageophysics.github.io/MTGeophysics.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation (stable)"></a>
	<a href="https://juliageophysics.github.io/MTGeophysics.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation (dev)"></a>
	<a href="https://joss.theoj.org/papers/e45b75b003b013751a4a2e1a51314103"><img src="https://joss.theoj.org/papers/e45b75b003b013751a4a2e1a51314103/status.svg" alt="JOSS status"></a>
	<a href="https://juliahub.com/ui/Packages/General/MTGeophysics"><img src="https://juliahub.com/docs/General/MTGeophysics/stable/version.svg" alt="Registry version"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/releases"><img src="https://img.shields.io/github/v/release/JuliaGeophysics/MTGeophysics.jl?label=release&color=blue" alt="Latest release"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/CI.yml/badge.svg" alt="CI"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/Documenter.yml"><img src="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/Documenter.yml/badge.svg" alt="Documentation build"></a>
	<a href="https://juliapkgstats.com/pkg/MTGeophysics"><img src="https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FMTGeophysics&query=total_requests&label=total%20downloads&color=389826" alt="Total downloads"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/stargazers"><img src="https://img.shields.io/github/stars/JuliaGeophysics/MTGeophysics.jl?label=stars&color=dfb317" alt="GitHub stars"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/forks"><img src="https://img.shields.io/github/forks/JuliaGeophysics/MTGeophysics.jl?label=forks&color=97ca00" alt="GitHub forks"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/discussions"><img src="https://img.shields.io/github/discussions/JuliaGeophysics/MTGeophysics.jl?label=discussions" alt="Discussions"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/issues"><img src="https://img.shields.io/github/issues/JuliaGeophysics/MTGeophysics.jl?label=open%20issues" alt="Open issues"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/issues?q=is%3Aissue+is%3Aclosed"><img src="https://img.shields.io/github/issues-closed/JuliaGeophysics/MTGeophysics.jl?label=closed%20issues&color=8250df" alt="Closed issues"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/main/CONTRIBUTING.md"><img src="https://img.shields.io/badge/contributing-guide-blueviolet" alt="Contributing guide"></a>
</p>

---


## Features

- 1D (exact) and 2D (finite-difference) MT forward solvers with Fréchet derivatives and their transpose
- 2D topography and water, `topo.dat`, and a mesh tool (batch or GLMakie)
- [2D and 1D inversion](docs/src/gaussnewton2d.md) from ModEM-style files: Gauss–Newton and NLCG (six files), and [VFSA](docs/src/inversion2d.md) with a mask and ensemble uncertainty (five files)
- 3D VFSA inversion with ensemble uncertainty quantification
- ModEM 3D model and data I/O
- Interactive 3D slice viewers (GLMakie) with GIS overlays and coordinate reprojection
- Shapefile export for GIS integration

## Requirements

- [Julia](https://julialang.org) 1.10 or newer
- OpenGL, for the interactive 3D viewers (GLMakie)

## Installation

**As a package** — to use MTGeophysics.jl from your own project or scripts. It is registered in the Julia General registry.

Install it into a **dedicated project environment**:

```julia
julia> ]  # press ] to enter the Pkg REPL
pkg> activate @mtgeophysics   # a named shared environment; or `activate .` for the current folder
pkg> add MTGeophysics
```

or equivalently, non-interactively:

```bash
julia --project=@mtgeophysics -e 'using Pkg; Pkg.add("MTGeophysics")'
```

> [!TIP]
> As a general Julia best practice, avoid installing packages into your default
> (global) environment. A dedicated per-project environment keeps dependencies
> isolated and reproducible, and avoids slow, unexpected version changes across
> unrelated packages you already have installed.

**From a clone** — to run the bundled examples, helper scripts, and benchmarks, or to develop the package:

```bash
git clone https://github.com/JuliaGeophysics/MTGeophysics.jl.git
cd MTGeophysics.jl/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The `julia --project=.` prefix used throughout the examples below activates that cloned environment.

## Getting started

Generate the 2D-IV benchmark. It is COMEMI 2D-III under 100–200 m of Finnish-style relief, with a
lake in the padding. Its noisy data come from a fine mesh, while the start and prior models, the
covariance and the mask sit on a coarser inversion mesh with the same topography. Then invert it with
Gauss–Newton or NLCG (six ModEM-style files), or with VFSA (five files: a mask, no covariance or prior):

```bash
julia --project=. helpers/benchmarks_2D.jl                    # examples/data/2D-IV
julia --project=. examples/run_inv2D.jl GN                     # or NLCG
julia --project=. -t 10 examples/run_vfsa2D.jl                 # threaded chains, ensemble uncertainty
```

Each run writes `run_YYYYmmdd_HHMMSS/` next to the data: `model.rho`, `data.pred`, `Summary.txt`,
the inputs, the plots and, for VFSA, the chains and the ensemble mean, median, spread and 5-95% range.

The same inversion from Julia:

```julia
using MTGeophysics

d = "examples/data/2D-IV"; c = "examples/ctrl/2D"
run = Invert2D("$d/model.start", "$d/data.dat", "$c/FwdCtrl", "$c/InvCtrl.GN", "$d/cov.ctrl", "$d/model.prior")
PlotInversion2D(run; true_model_path = "$d/model.true")
vrun = VFSA2D("$d/model.start", "$d/data.dat", "$c/FwdCtrl", "$c/InvCtrl.VFSA", "$d/mask.ctrl")
```

`helpers/benchmarks_2D.jl 2D-III` adds the same model without topography,
`helpers/benchmarks_1D.jl` a 1D sounding for `examples/run_inv1D.jl`, and
`examples/make_mesh2D.jl` builds the inputs for your own data file.

See the [documentation](https://juliageophysics.github.io/MTGeophysics.jl/dev/) for the full 1D/2D/3D workflows, ensemble statistics and convergence animations, ModEM I/O, configuration options, and the interactive 3D viewers.

# Research using this code 

- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., and Sen, M. K.: Building uncertainty-aware subsurface models with 3D magnetotelluric inversion, EGU General Assembly 2026, Vienna, Austria, 3–8 May 2026, EGU26-4367, https://doi.org/10.5194/egusphere-egu26-4367, 2026. 

- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026. 

- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., Xiao, L., and Sen, M. K.: Stochastic model exploration in three-dimensional inversion of magnetotelluric data, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026. 

- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, JuliaCon 2026. 

- Patzer, C., Mishra, P. K., and Kamm, J.: Studying the Wiborg Rapakivi Batholith in SE Finland, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026. 




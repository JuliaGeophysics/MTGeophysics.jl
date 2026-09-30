<h1 align="center">MTGeophysics.jl</h1>

<p align="center"><em>A software repository for magnetotelluric research and applications.</em></p>

<p align="center">
	<a href="https://juliageophysics.github.io/MTGeophysics.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/CI.yml/badge.svg" alt="CI"></a>
	<a href="https://joss.theoj.org/papers/e45b75b003b013751a4a2e1a51314103"><img src="https://joss.theoj.org/papers/e45b75b003b013751a4a2e1a51314103/status.svg" alt="JOSS status"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/releases"><img src="https://img.shields.io/github/v/release/JuliaGeophysics/MTGeophysics.jl?label=release&color=blue" alt="Latest release"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
</p>

## Requirements

- [Julia](https://julialang.org) 1.12 or newer (tested on 1.12 and 1.13)
- OpenGL, for the interactive 3D viewers (GLMakie)

## Installation

```julia
julia> ]
pkg> activate @mtgeophysics
pkg> add MTGeophysics
```

To get started, see the [**Getting Started**](https://juliageophysics.github.io/MTGeophysics.jl/dev/getting_started) guide. It covers installing from a clone, running the examples, and a first forward model and inversion. The [documentation](https://juliageophysics.github.io/MTGeophysics.jl/dev/) covers the full 1D/2D/3D workflows.

Contributions are welcome, see [CONTRIBUTING.md](CONTRIBUTING.md).

## Research using this code

- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., and Sen, M. K.: Building uncertainty-aware subsurface models with 3D magnetotelluric inversion, EGU General Assembly 2026, Vienna, Austria, 3–8 May 2026, EGU26-4367, https://doi.org/10.5194/egusphere-egu26-4367, 2026.
- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026.
- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., Xiao, L., and Sen, M. K.: Stochastic model exploration in three-dimensional inversion of magnetotelluric data, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026.
- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, JuliaCon 2026.
- Patzer, C., Mishra, P. K., and Kamm, J.: Studying the Wiborg Rapakivi Batholith in SE Finland, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026.

<h1 align="center">MTGeophysics.jl</h1>

<p align="center"><em>A software repository for magnetotelluric research and applications.</em></p>

<p align="center">
	<a href="https://juliageophysics.github.io/MTGeophysics.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaGeophysics/MTGeophysics.jl/actions/workflows/CI.yml/badge.svg" alt="CI"></a>
	<a href="https://joss.theoj.org/papers/e45b75b003b013751a4a2e1a51314103"><img src="https://joss.theoj.org/papers/e45b75b003b013751a4a2e1a51314103/status.svg" alt="JOSS status"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/releases"><img src="https://img.shields.io/github/v/release/JuliaGeophysics/MTGeophysics.jl?label=release&color=blue" alt="Latest release"></a>
	<a href="https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
</p>

## Installation

> [!TIP]
> **New to Julia?** Install it first from [julialang.org/install](https://julialang.org/install/), which sets up Julia and its version manager, [juliaup](https://github.com/JuliaLang/juliaup), in one step. Then open a terminal, type `julia`, and follow the steps below.

Install MTGeophysics.jl into its own environment (press `]` to enter the package manager):

```julia
julia> ]
pkg> activate @mtgeophysics
pkg> add MTGeophysics
```

Then load it in a Julia session or script:

```julia
julia> using MTGeophysics
```

In a new session, run `pkg> activate @mtgeophysics` again before `using MTGeophysics`.

To get started, see the [**Getting Started**](https://juliageophysics.github.io/MTGeophysics.jl/dev/getting_started) guide. It covers installing from a clone, running the examples, and a first forward model and inversion. The [documentation](https://juliageophysics.github.io/MTGeophysics.jl/dev/) covers the full 1D/2D/3D workflows.

Contributions are welcome, see [CONTRIBUTING.md](CONTRIBUTING.md).

> [!NOTE]
> **Also check out** other JuliaGeophysics packages: [Timekeepers.jl](https://github.com/JuliaGeophysics/Timekeepers.jl) for time-series analysis in electromagnetic geophysics, and [TerraScope.jl](https://github.com/JuliaGeophysics/TerraScope.jl) for looking at diverse geoscientific data together.

## Publications using or citing this package

- Rulff, P., Castillo-Reyes, O., Deleersnyder, W., Carrizo Mascarell, M., Minsley, B. J., and King, J.: Computational electromagnetic geophysics for groundwater system studies: A review of established practices and recent advances, Journal of Hydrology, 677, 135542, https://doi.org/10.1016/j.jhydrol.2026.135542, 2026.
- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, JuliaCon 2026.
- Mishra, P. K.: MTGeophysics.jl: A software repository for magnetotelluric research and application, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026.
- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., Xiao, L., and Sen, M. K.: Stochastic model exploration in three-dimensional inversion of magnetotelluric data, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026.
- Patzer, C., Mishra, P. K., and Kamm, J.: Studying the Wiborg Rapakivi Batholith in SE Finland, 27th Electromagnetic Induction Workshop (EMIW 2026), St. John's, Newfoundland and Labrador, Canada, 2026.
- Mishra, P. K., Kamm, J., Patzer, C., Autio, U., and Sen, M. K.: Building uncertainty-aware subsurface models with 3D magnetotelluric inversion, EGU General Assembly 2026, Vienna, Austria, 3–8 May 2026, EGU26-4367, https://doi.org/10.5194/egusphere-egu26-4367, 2026.

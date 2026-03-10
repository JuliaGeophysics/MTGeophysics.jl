# MTGeophysics.jl

*A software repository for magnetotelluric geophysics research and applications.*

## Features

- **1D forward modelling** — layered-earth analytical and finite-difference solvers
- **2D forward modelling** — TE/TM mode solvers on tensor meshes
- **2D VFSA inversion** — stochastic inversion with multi-chain ensemble sampling
- **3D model viewer** — interactive GLMakie slice viewer with GIS overlays

## Quick start

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl
julia --project=. Examples/response_1d.jl
julia --project=. Examples/response_2d.jl
julia --project=. Examples/run_vfsa2dmt.jl
```

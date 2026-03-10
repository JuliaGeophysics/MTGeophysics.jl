# MTGeophysics.jl

*A software repository for magnetotelluric research and applications*



**[Documentation](https://pankajkmishra.github.io/MTGeophysics.jl)**

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

```
julia --project=. Helpers/benchmarks_2d.jl
julia --project=. examples/run_vfsa2dmt.jl
```

![2D VFSA convergence](images/vfsa2d_convergence.gif)

## License

[MIT](LICENSE.md)


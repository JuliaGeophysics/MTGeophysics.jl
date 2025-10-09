# MTGeophysics

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeophysics.github.io/MTGeophysics.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeophysics.github.io/MTGeophysics.jl/dev)

Julia package for magnetotelluric (MT) geophysics research and applications.

## Features

- MT 1D forward modeling with analytical and finite-difference methods
- ModEM data loading and processing
- 3D model visualization with GLMakie
- Chi-squared and RMS calculations for inversion analysis
- UBC format mesh and model file generation

## Installation

```julia
using Pkg
Pkg.add("MTGeophysics")
```

## Quick Start

```julia
using MTGeophysics

# Define a 3-layer earth model
frequencies = 10.0 .^ range(-2, 2, length=50)
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]
resistivities = [100.0, 10.0, 1000.0]

# Compute MT response
rho_a, phi = MT1D_response(frequencies, layer_boundaries, resistivities)
```

## Documentation

For detailed documentation, tutorials, and API reference, please visit the [documentation](https://JuliaGeophysics.github.io/MTGeophysics.jl/stable). 


![Alt text](images/MT1D_comparison.png) 
![Alt text](images/Model3D.png)
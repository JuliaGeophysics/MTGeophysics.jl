# MTGeophysics.jl

Julia package for magnetotelluric (MT) geophysics research and applications.

## Features

- **MT 1D Forward Modeling**: Analytical and finite-difference methods for computing magnetotelluric responses
- **ModEM Data Loading**: Read and process ModEM data files
- **3D Model Visualization**: Interactive plotting of 3D resistivity models
- **Inversion Analysis**: Chi-squared and RMS calculations for inversion quality assessment
- **UBC Format Support**: Generate mesh and model files in UBC format

## Package Overview

MTGeophysics.jl provides a comprehensive suite of tools for magnetotelluric (MT) geophysics research. The package implements both analytical and numerical solutions for 1D MT forward modeling, supports popular data formats like ModEM and UBC, and includes visualization capabilities for 3D models.

### Key Components

- **Forward Modeling**: Multiple implementations for MT 1D response calculations
  - Analytical solution using recursive algorithms
  - Finite-difference methods for numerical stability
  - Diagonal finite-difference for efficiency

- **Data Structures**: Native support for ModEM data format
  - Load and manipulate impedance and tipper data
  - Calculate apparent resistivity and phase
  - Error handling and data validation

- **Visualization**: Interactive 3D model viewer (requires GLMakie)
  - Slice-based visualization
  - Configurable color scales and ranges
  - Toggle between linear and logarithmic scales

- **Inversion Tools**: Quality metrics for inversion results
  - Chi-squared calculation
  - RMS computation
  - Component-wise analysis

## Quick Start

```julia
using MTGeophysics

# Define a 3-layer earth model
frequencies = 10.0 .^ range(-2, 2, length=50)  # 0.01 to 100 Hz
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]  # depths in meters
resistivities = [100.0, 10.0, 1000.0]  # resistivity in Ω⋅m

# Compute MT response
rho_a, phi = MT1D_response(frequencies, layer_boundaries, resistivities)
```

See the [Tutorials](@ref) section for more examples.

## Citation

If you use MTGeophysics.jl in your research, please cite:

```
@software{mtgeophysics2025,
  author = {JuliaGeophysics community and Pankaj K Mishra},
  title = {MTGeophysics.jl: Julia package for magnetotelluric geophysics},
  year = {2025},
  url = {https://github.com/JuliaGeophysics/MTGeophysics.jl}
}
```

## Contents

```@contents
Pages = ["installation.md", "tutorials.md", "api.md"]
Depth = 2
```

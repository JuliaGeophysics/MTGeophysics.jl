# Installation

## Prerequisites

MTGeophysics.jl requires Julia 1.6 or later. If you don't have Julia installed, download it from [julialang.org](https://julialang.org/downloads/).

## Installing MTGeophysics.jl

### From the Julia REPL

The package can be installed using Julia's package manager:

```julia
using Pkg
Pkg.add("MTGeophysics")
```

Or in the Pkg REPL mode (press `]` in the Julia REPL):

```
pkg> add MTGeophysics
```

### Development Version

To install the latest development version from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaGeophysics/MTGeophysics.jl")
```

## Optional Dependencies

MTGeophysics.jl has some optional dependencies for additional functionality:

### Visualization (GLMakie)

For 3D model visualization using the `gl_modem_viewer` function:

```julia
using Pkg
Pkg.add("GLMakie")
```

### Advanced Forward Modeling (LinearSolve)

For finite-difference forward modeling with advanced linear solvers:

```julia
using Pkg
Pkg.add("LinearSolve")
```

## Verifying the Installation

To verify that MTGeophysics.jl is installed correctly:

```julia
using MTGeophysics

# Create a simple test
frequencies = [1.0, 10.0, 100.0]
mesh = [0.0, 1000.0, 2000.0]
resistivities = [100.0, 1000.0]

# Run a simple forward model
rho_a, phi = MT1D_response(frequencies, mesh, resistivities)
println("Installation successful! Computed $(length(rho_a)) frequency points.")
```

If this runs without errors, your installation is working correctly.

## Troubleshooting

### Package Not Found

If you get a "package not found" error, make sure you have the General registry added:

```julia
using Pkg
Pkg.Registry.add("General")
```

### Dependency Issues

If you encounter issues with dependencies, try updating your packages:

```julia
using Pkg
Pkg.update()
```

### Platform-Specific Issues

#### Linux

On some Linux distributions, you may need to install additional graphics libraries for GLMakie visualization. For Ubuntu/Debian:

```bash
sudo apt-get install libglfw3 libglfw3-dev
```

#### macOS

GLMakie should work out of the box on macOS with XQuartz installed.

#### Windows

GLMakie should work out of the box on Windows.

## Getting Help

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/JuliaGeophysics/MTGeophysics.jl/issues) page
2. Open a new issue with a minimal reproducible example
3. Ask on the [Julia Discourse](https://discourse.julialang.org/) with the `geophysics` tag

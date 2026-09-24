# Getting Started

## Requirements

- Julia 1.10+ (developed on 1.12.4)
- OpenGL for interactive 3D viewers (GLMakie)

## Installation

There are two ways to install, depending on how you intend to use the package.

### 1. From the registry

MTGeophysics.jl is registered in the Julia General registry. Use this if you want to call the
package from your own project or scripts. Install it into a dedicated project environment:

```julia
julia> ]  # press ] to enter the Pkg REPL
pkg> activate @mtgeophysics   # a named shared environment; or `activate .` for the current folder
pkg> add MTGeophysics
```

or, non-interactively:

```bash
julia --project=@mtgeophysics -e 'using Pkg; Pkg.add("MTGeophysics")'
```

To add it to a specific project directory instead:

```bash
julia --project=/path/to/your/project -e 'using Pkg; Pkg.add("MTGeophysics")'
```

!!! tip
    As a general Julia best practice, avoid installing packages into your default (global)
    environment. A dedicated per-project environment keeps dependencies isolated and
    reproducible, and avoids slow, unexpected version changes across unrelated packages you
    already have installed.

### 2. From a clone

Use this if you want the bundled `examples/` and `helpers/` scripts, the benchmark generators,
and the test suite, or if you plan to develop the package:

```bash
git clone https://github.com/JuliaGeophysics/MTGeophysics.jl.git
cd MTGeophysics.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` installs the exact dependency versions recorded in `Manifest.toml`. All
example commands in this documentation assume you are in the clone root and pass
`--project=.` so that this environment is active.

To develop the clone while using it from another environment, use `Pkg.develop`:

```bash
julia --project=/path/to/your/project -e 'using Pkg; Pkg.develop(path="/path/to/MTGeophysics.jl")'
```

## Verify

```bash
julia --project=. -e 'using MTGeophysics; println("OK")'
julia --project=. test/runtests.jl
```

## Generate benchmarks

Before running the examples, generate the synthetic benchmark data:

```bash
julia --project=. helpers/benchmarks_1D.jl
julia --project=. helpers/benchmarks_2D.jl
```

This creates `examples/data/1D-I` (a five-layer sounding) and `examples/data/2D-IV` (COMEMI
2D-III under topography with a lake). `helpers/benchmarks_2D.jl 2D-I 2D-II 2D-III` adds the
flat COMEMI cases. Each 2D folder holds `model.true`, `data.dat`, `model.start`,
`model.prior`, `cov.ctrl` (GN, NLCG) and `mask.ctrl` (VFSA), and 2D-IV also `topo.dat`.
1D-I holds only `data.dat` and `model.true`, since the 1D inversion lays out its own mesh.
The control files ship in `examples/ctrl/1D` and `examples/ctrl/2D`.

## First session

```bash
julia --project=. helpers/benchmarks_1D.jl
julia --project=. helpers/benchmarks_2D.jl
julia --project=. examples/run_fwd1D.jl
julia --project=. examples/run_inv1D.jl GN
julia --project=. examples/run_inv1D.jl VFSA
julia --project=. examples/run_fwd2D.jl
julia --project=. examples/run_inv2D.jl GN
julia --project=. -t 10 examples/run_vfsa2D.jl
```

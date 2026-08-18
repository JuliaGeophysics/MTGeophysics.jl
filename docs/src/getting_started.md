# Getting Started

## Requirements

- Julia 1.10+ (developed on 1.12.4)
- OpenGL for interactive 3D viewers (GLMakie)

## Installation

There are two ways to install, depending on how you intend to use the package.

### 1. From the registry

MTGeophysics.jl is registered in the Julia General registry. Use this if you want to call the
package from your own project or scripts:

```julia
julia> ]  # press ] to enter the Pkg REPL
pkg> add MTGeophysics
```

or, non-interactively:

```bash
julia -e 'using Pkg; Pkg.add("MTGeophysics")'
```

To add it to a specific project environment rather than your default one:

```bash
julia --project=/path/to/your/project -e 'using Pkg; Pkg.add("MTGeophysics")'
```

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

This creates:

- `examples/0Layered1D/` — 1D layered benchmark
- `examples/0COMEMI2D-I/`, `0COMEMI2D-II/`, `0COMEMI2D-III/` — 2D COMEMI benchmarks

## First session

```bash
julia --project=. helpers/benchmarks_1D.jl
julia --project=. helpers/benchmarks_2D.jl
julia --project=. examples/run_fwd1D.jl
julia --project=. examples/run_fwd2D.jl
julia --project=. examples/run_vfsa2D.jl
```

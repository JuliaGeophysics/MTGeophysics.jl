# Getting Started

## Requirements

- Julia 1.10+ (developed on 1.12.4)
- OpenGL for interactive 3D viewers (GLMakie)

## Installation

```bash
git clone https://github.com/pankajkmishra/MTGeophysics.jl.git
cd MTGeophysics.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Verify

```bash
julia --project=. -e 'using MTGeophysics; println("OK")'
julia --project=. test/runtests.jl
```

## Generate benchmarks

Before running the examples, generate the synthetic benchmark data:

```bash
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl
```

This creates:

- `Examples/0Layered1D/` — 1D layered benchmark
- `Examples/0COMEMI2D-I/`, `0COMEMI2D-II/`, `0COMEMI2D-III/` — 2D COMEMI benchmarks

## First session

```bash
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl
julia --project=. Examples/response_1d.jl
julia --project=. Examples/response_2d.jl
julia --project=. Examples/run_vfsa2dmt.jl
```

# Getting Started

In this section we install MTGeophysics.jl, check that it works, and run a first forward model and
inversion on a small synthetic example. It takes about ten minutes.

## What you need

- Julia 1.10 or newer
- OpenGL, only for the interactive 3D viewers (GLMakie)
- ModEM compiled with MPI, only for 3D inversion (see [3D VFSA](inversion/3d_vfsa.md))

## Install the package

There are two ways to install MTGeophysics.jl. Pick the one that matches how you want to use it.

### Option 1: as a package

Use this if you want to call MTGeophysics from your own scripts. We install it into its own
environment:

```julia
julia> ]                      # press ] to enter the package manager
pkg> activate @mtgeophysics   # a named environment; or `activate .` for the current folder
pkg> add MTGeophysics
```

!!! tip
    Avoid installing packages into your default (global) Julia environment. A separate environment
    per project keeps versions stable and your work reproducible.

### Option 2: from a clone

Use this if you want the example scripts, the benchmark generators and the tests, or if you plan to
change the code:

```bash
git clone https://github.com/JuliaGeophysics/MTGeophysics.jl.git
cd MTGeophysics.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` installs the exact versions listed in `Manifest.toml`. All commands in this
documentation assume you are in the clone folder and pass `--project=.`.

## Check the installation

```bash
julia --project=. -e 'using MTGeophysics; println("OK")'
```

If this prints `OK`, you are ready to go. The first `using` compiles the package and can take a few
minutes.

## Create the example data

The examples use small synthetic data sets. Let's create them:

```bash
julia --project=. helpers/benchmarks_1D.jl
julia --project=. helpers/benchmarks_2D.jl
```

This writes two folders:

- `examples/data/1D-I`: one sounding over a five-layer earth.
- `examples/data/2D-IV`: a 2D profile over hills with a lake, based on the COMEMI 2D-III model.

Each folder holds the true model and noisy data computed from it. The 2D folder also holds a start
model and the other inversion inputs. The control files the examples use are in `examples/ctrl/`.

## Your first run

Now we can compute forward responses and run two inversions:

```bash
julia --project=. examples/run_fwd1D.jl        # 1D forward response
julia --project=. examples/run_inv1D.jl GN     # 1D Gauss–Newton inversion
julia --project=. examples/run_fwd2D.jl        # 2D forward response
julia --project=. examples/run_inv2D.jl GN     # 2D Gauss–Newton inversion
```

Each inversion writes a new folder, `run_YYYYmmdd_HHMMSS/`, next to the data. Inside you will find the
recovered model (`model.rho`), its predicted data (`data.pred`), a short `Summary.txt` and a `plots/`
folder. Open the plots to see how well the model fits the data and how it compares with the true
model.

## Where to go next

- Have your own data? Start with [Data and model files](data/files.md) and
  [1D and 2D meshes](data/mesh2d.md).
- Want to know what a forward run computes? See [Forward](forward/2d.md).
- Ready to invert? See [Inversion](inversion/2d_deterministic.md).
- Working in 3D? See [3D meshes](data/mesh3d.md), [3D VFSA](inversion/3d_vfsa.md) and the
  [3D viewers](visualisation/3d_models.md).

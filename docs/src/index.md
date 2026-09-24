# MTGeophysics.jl

*A software repository for magnetotelluric geophysics research and applications.*

MTGeophysics.jl is part of the [JuliaGeophysics ecosystem](https://github.com/JuliaGeophysics) and is intended for both research and real-world applications. It provides reusable forward-modelling, inversion, and visualization components so you can prototype new machine-learning methods, inversion strategies, and data-analysis ideas quickly without rebuilding core MT tooling from scratch. More broadly, JuliaGeophysics aims to build a tightly integrated yet modular ecosystem for multiphysics workflows, multisource data integration, and uncertainty quantification.

## Features

- **1-D modelling and inversion** — exact layered-earth responses and Fréchet derivatives; each site of a data file inverted on its own (Gauss–Newton or VFSA) on a skin-depth mesh from `MakeMesh1D`, with one small control file. See [1D](forward1d.md).
- **2-D forward modelling** — TE/TM finite-difference (box integration) solver with topography and water, ModEM model and data files, `fwd.ctrl`, Fréchet derivatives and their transpose. See [2D forward](forward2d.md) and [topography](topography2d.md).
- **2-D deterministic inversion** — Gauss–Newton and NLCG from six ModEM-style files, with covariance masks, fixed air and water, and a standard run folder. See [deterministic inversion](gaussnewton2d.md).
- **2-D VFSA inversion** — very fast simulated annealing from five files (a mask, no covariance or prior), threaded chains and ensemble uncertainty (mean, median, std, 5–95 %). See [VFSA](inversion2d.md).
- **2-D mesh tool** — inversion inputs from a data file and `topo.dat`, batch or GLMakie window. See [mesh tool](mesh2d.md).
- **3-D ModEM data & model I/O** — full impedance-tensor + tipper reader/writer for ModEM and WS3D formats, apparent-resistivity/phase derivation, and χ²/RMS misfit evaluation.
- **3-D VFSA inversion** — 3-D VFSA engine using ModEM as the external MPI-parallel forward solver, with RBF control-point parameterisation, padding decay, and multi-chain ensemble analysis.
- **3-D model utilities** — headless core/padding detection, depth truncation, and core sub-array extraction (no GLMakie required).
- **3-D interactive visualisation** (GLMakie) — volume-slice viewer, XY/XZ/YZ cross-section browsers, polygon-based model editor, bulk deep-layer editor; all with CRS reprojection, shapefile overlays, and high-resolution export.

## Quick start

Requires Julia 1.10 or newer. To use the package from your own project, `pkg> add MTGeophysics`;
the commands below assume a clone of the repository, which also provides the examples and helper
scripts — see [Getting Started](getting_started.md) for both routes.

```bash
# One-time setup
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Generate synthetic benchmark models
julia --project=. helpers/benchmarks_1D.jl
julia --project=. helpers/benchmarks_2D.jl

# 1-D forward response and inversion
julia --project=. examples/run_fwd1D.jl
julia --project=. examples/run_inv1D.jl GN
julia --project=. examples/run_inv1D.jl VFSA

# 2-D forward response, deterministic and VFSA inversion
julia --project=. examples/run_fwd2D.jl
julia --project=. examples/run_inv2D.jl GN
julia --project=. -t 10 examples/run_vfsa2D.jl

# 3-D VFSA inversion (requires ModEM + MPI on PATH)
julia --project=. examples/run_vfsa3D.jl

# Interactive 3-D viewers (requires GLMakie)
julia --project=. examples/plot_model_XYZ.jl <model.ws> <data.dat>
julia --project=. examples/plot_model_XY_slices.jl <model.ws> <data.dat> EPSG:32610
julia --project=. examples/plot_model_XY_with_shapefiles.jl <model.ws> <data.dat> EPSG:32610
julia --project=. examples/plot_model_XZ_slices.jl <model.ws> <data.dat>
julia --project=. examples/plot_model_YZ_slices.jl <model.ws> <data.dat>
```

The 3-D viewers and model editors are exported by the package, so they can also
be called directly — `PlotModelXYZ(model, data)`, `EditModelByLayers(model)` —
without a repository checkout. The Cascadia model used in the
[3D Visualization](visualization3d.md) examples is not distributed with the
package; see [Example data](visualization3d.md#Example-data) for the download
link.

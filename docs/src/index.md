# MTGeophysics.jl

*A software repository for magnetotelluric geophysics research and applications.*

MTGeophysics.jl is part of the [JuliaGeophysics ecosystem](https://github.com/JuliaGeophysics) and is intended for both research and real-world applications. It provides reusable forward-modelling, inversion, and visualization components so you can prototype new machine-learning methods, inversion strategies, and data-analysis ideas quickly without rebuilding core MT tooling from scratch. More broadly, JuliaGeophysics aims to build a tightly integrated yet modular ecosystem for multiphysics workflows, multisource data integration, and uncertainty quantification.

## Features

- **1-D forward modelling** — analytical (recursive impedance) and finite-difference solvers for layered-earth models, with mesh building, model/data I/O, and plotting.
- **2-D forward modelling** — TE/TM finite-volume solver on tensor meshes, file-driven workflows, model/data I/O, misfit metrics, and pseudo-section/curve plotting.
- **2-D VFSA inversion** — Very Fast Simulated Annealing with Gaussian-RBF parameterisation, multi-chain ensemble sampling, and automatic ensemble statistics (mean, median, std).
- **3-D ModEM data & model I/O** — full impedance-tensor + tipper reader/writer for ModEM and WS3D formats, apparent-resistivity/phase derivation, and χ²/RMS misfit evaluation.
- **3-D VFSA inversion** — 3-D VFSA engine using ModEM as the external MPI-parallel forward solver, with RBF control-point parameterisation, padding decay, and multi-chain ensemble analysis.
- **3-D model utilities** — headless core/padding detection, depth truncation, and core sub-array extraction (no GLMakie required).
- **3-D interactive visualisation** (GLMakie) — volume-slice viewer, XY/XZ/YZ cross-section browsers, polygon-based model editor, bulk deep-layer editor; all with CRS reprojection, shapefile overlays, and high-resolution export.

## Quick start

```bash
# One-time setup
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Generate synthetic benchmark models
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl

# 1-D forward response
julia --project=. examples/response_1d.jl

# 2-D forward response
julia --project=. examples/response_2d.jl

# 2-D VFSA inversion
julia --project=. examples/run_vfsa2dmt.jl

# 3-D VFSA inversion (requires ModEM + MPI on PATH)
julia --project=. examples/run_vfsa3dmt.jl

# Interactive 3-D viewers (requires GLMakie)
julia --project=. examples/plot_model_3D.jl
julia --project=. examples/plot_XY_slices.jl
julia --project=. examples/plot_XZ_slices.jl
julia --project=. examples/plot_YZ_slices.jl
```

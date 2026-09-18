# API Reference

## Types

| Type | Description |
|:-----|:------------|
| `MT1DMesh` | 1D layered-earth model and finite-difference discretisation |
| `MT1DDataSpec` | 1D survey specification (frequencies and errors) |
| `MT1DResponse` | 1D forward response (impedance, apparent resistivity, phase) |
| `MT2DMesh` | 2D tensor mesh (nodes, cells, receivers, frequencies) |
| `MT2DResponse` | 2D forward response (TE/TM apparent resistivity, phase, impedance) |
| `ModelFile2D` | Parsed 2D model file |
| `DataFile2D` | Parsed 2D data file |
| `FitSummary2D` | Misfit summary (χ², RMS, count) |
| `Inv2DOptions` | Algorithm-independent 2D inversion controls |
| `Inv2DResult` | Recovered model, fit, history, and stopping reason |
| `AbstractInversion2D` | Supertype of 2D deterministic inversion algorithms |
| `GaussNewton2DConfig` | Damped Gauss-Newton algorithm settings |
| `GaussNewton2DResult` | `Inv2DResult` produced by Gauss-Newton |
| `NLCG2DConfig` | Preconditioned NLCG algorithm settings |
| `NLCG2DResult` | `Inv2DResult` produced by NLCG |
| `VFSA2DMTConfig` | VFSA 2D inversion parameters |
| `VFSA2DMTParams` | VFSA 2D run paths and configuration |
| `VFSA3DMTConfig` | VFSA 3D inversion parameters |
| `WS3DModel` | 3D resistivity model in WS3D format |
| `RBFMap` | Gaussian-RBF mapping for 3D parameterization |
| `Data` / `Model` | 3D data and model containers |

## 1D Functions

| Function | Description |
|:---------|:------------|
| `BuildMesh1D(t, ρ)` | Build a 1D mesh from layer thicknesses and resistivities |
| `solve_mt1d_analytical(f, ρ, t)` | Analytical recursive impedance solver |
| `solve_mt1d_fd(f, mesh)` | Finite-difference 1D solver |
| `Forward1D(f, ρ, t)` | Run both solvers |
| `ForwardSolve1D(model, data)` | File-based forward solve |
| `PlotData1D(ref, pred)` | Plot observed vs predicted |
| `PlotModel1D(model)` | Plot the layered model |
| `plot_mt1d_data(responses)` | Plot response curves |
| `plot_mt1d_model(mesh)` | Plot model structure |

## 2D Functions

| Function | Description |
|:---------|:------------|
| `BuildMesh2D(; ...)` | Build a 2D tensor mesh |
| `BuildMesh2D(; ...)` | Profile mesh; vertical core uniform to one skin depth of the lowest frequency |
| `mt2d_skin_depth(ρ, f)` | Skin depth in metres |
| `mt2d_skin_depth_layers(f; ...)` | Uniform-core plus geometric-padding ground layers |
| `build_default_mt2d_mesh()` | Default COMEMI mesh |
| `build_mt2d_halfspace_model(mesh)` | Uniform resistivity model |
| `build_mt2d_layered_model(mesh)` | Layered model |
| `build_mt2d_block_model(mesh)` | Model with rectangular anomalies |
| `run_mt2d_forward(mesh, ρ)` | TE/TM forward solver |
| `ForwardSolve2D(model, data)` | File-based forward solve |
| `PlotData2D(data)` | Plot data maps and site curves |
| `PlotModel2D(model)` | Plot the 2D model |
| `plot_mt2d_model(mesh, ρ)` | Plot model cross-section |
| `plot_mt2d_data_maps(response)` | Plot TE/TM maps |
| `plot_mt2d_site_curves(response)` | Plot per-site curves |
| `plot_mt2d_mesh(mesh; region)` | Plot mesh edges, air, uniform core, and skin depths |
| `plot_mt2d_data_fit(obs, pred)` | Observed vs predicted curves at selected sites |
| `plot_inv2d_convergence(history)` | RMS and objective terms per iteration |
| `chi2_rms2d(obs, pred)` | Compute misfit |
| `Invert2D(mesh, initial, observed; algorithm, options, ...)` | Regularized impedance inversion with any algorithm |
| `Invert2D(model_path, data_path; output_dir, ...)` | File-driven inversion workflow |
| `GaussNewton2D(...)` | `Invert2D` with `algorithm = GaussNewton2DConfig()` |
| `NLCG2D(...)` | `Invert2D` with `algorithm = NLCG2DConfig()` |
| `inv2d_frechet(problem, state)` | Data-weighted Fréchet derivative C_D^{-1/2} G in log10 resistivity |
| `inv2d_gradient(problem, state[, J])` | Objective gradient (explicit or adjoint) |
| `FrechetDerivative2D(mesh, ρ; ...)` | Explicit Fréchet derivative G = ∂g/∂m |
| `ApplyFrechet2D(mesh, ρ, δm; ...)` | Tangent linear application δd = G δm |
| `ApplyFrechetTranspose2D(mesh, ρ, δd̂; ...)` | Transpose application δm̂ = Gᵗ δd̂ |

## VFSA Inversion

| Function | Description |
|:---------|:------------|
| `VFSA2DMT(params)` | Run the 2D VFSA inversion workflow |
| `AnalyseEnsemble2D(chains)` | Compute 2D ensemble statistics |
| `VFSA3DMT(model; dobs_path, cfg)` | Run the 3D VFSA inversion workflow |
| `AnalyseEnsemble3D(dir)` | Compute 3D ensemble mean/median/std |
| `core_statistics(cores)` | Element-wise mean, median, std over 3D cubes |

## 3D Functions

| Function | Description |
|:---------|:------------|
| `load_data_modem(path)` | Load 3D data file |
| `load_model_modem(path)` | Load 3D model file |
| `write_model_modem(path, model)` | Write 3D model file |
| `chi2_and_rms(obs, pred)` | Compute 3D misfit |
| `phase_tensor(Zxx, Zxy, Zyx, Zyy)` | Phase tensor invariants for one site and period |
| `induction_vector(Tzx, Tzy; convention)` | Induction arrow components from the tipper |
| `phase_tensors_from_data(d)` | Phase tensors for every period and site |
| `induction_vectors_from_data(d; convention)` | Induction vectors for every period and site |
| `has_tipper_data(d)` | Whether a data set carries usable tippers |
| `write_ptiv_gis(data_file; ...)` | Write phase tensor / induction vector shapefiles per period |
| `PlotPTIVMap(data_file; ...)` | Interactive phase tensor and induction vector map |

## WS3D Model I/O

| Function / Type | Description |
|:---------|:------------|
| `WS3DModel` | 3D resistivity model in WS3D format (log₁₀ internal) |
| `load_ws3d_model(path)` | Load a WS3D model file |
| `read_ws3d_model(path)` | Alias for `load_ws3d_model` |
| `write_ws3d_model(path, ...)` | Write a WS3D model file |

## Core Utilities

| Function | Description |
|:---------|:------------|
| `edges_from_centers(c)` | Cell-edge coordinates from cell centres |
| `core_indices(c; tol)` | Index range of the unpadded core cells |
| `z_indices_for_max_depth(zc, d)` | Depth-limited vertical index range |
| `lateral_core_ranges(m; tol)` | `(ix, iy)` core ranges for a model |
| `core_view(m; tol)` | View into the core resistivity block |
| `RBFMap` | Gaussian-RBF mapping structure |
| `build_rbf_map(m, ix, iy, n, rng)` | Build a 3D RBF control-point map |
| `apply_rbf_map!(delta, rbf, params)` | Apply RBF perturbations to a 3D field |

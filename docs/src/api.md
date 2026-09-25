# API Reference

The exported functions and types, grouped by topic. In the Julia REPL, `?name` shows the docstring
of a function.

## Types

| Type | Description |
|:-----|:------------|
| `MT2DMesh` | 2D tensor mesh (nodes, cells, receivers, frequencies, topographic air) |
| `FwdCtrl2D`, `InvCtrl2D`, `Cov2D` | `fwd.ctrl`, `inv.ctrl` (GN, NLCG) and covariance file contents |
| `VFSACtrl2D` | VFSA control file contents |
| `InvCtrl1D` | 1D inversion control (GN or VFSA) |
| `Topo2D` | Topography points (WGS84 lat, lon, elevation) |
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
| `VFSA2DConfig` | VFSA 2D (and 1D) inversion parameters; `VFSA2DConfig(ctrl::VFSACtrl2D)` from a control file |
| `VFSA3DMTConfig` | VFSA 3D inversion parameters |
| `WS3DModel` | 3D resistivity model in WS3D format |
| `RBFMap` | Gaussian-RBF mapping for 3D parameterization |
| `Data` / `Model` | 3D data and model containers |

## 1D Functions

| Function | Description |
|:---------|:------------|
| `mt1d_impedance(f, ρ, h)` | Layered-earth surface impedance by the recursion |
| `mt1d_layers(f; ...)`, `mt1d_skin_depth(ρ, f)` | Skin-depth layering, skin depth |
| `mt1d_frechet(f, ρ, h)`, `WriteFrechet1D(path, h, ρ, site)` | G = ∂Z/∂log10 ρ, and its file |
| `MakeMesh1D(data; ...)` | Skin-depth layering and background resistivity of every site |
| `mt1d_site_data(data, i; mode)` | One site as a 1D survey; `:DET` gives √det Z |
| `ReadInvCtrl1D`, `WriteInvCtrl1D` | 1D inversion control |
| `ForwardSolve1D(model, data; mode)` | File forward run, writes `data.pred` |
| `Invert1D(data, inv, meshes)` | 1D GN or VFSA inversion, every site on its own |
| `PlotInversion1D(run)` | Model, data fit and convergence plots per site |
| `PlotModel1D(model)`, `plot_mt1d_model(h, models)` | Resistivity-depth steps |
| `plot_mt1d_convergence(histories)` | GN or VFSA convergence |

## 2D Functions

| Function | Description |
|:---------|:------------|
| `ReadModel2D`, `WriteModel2D` | ModEM-layout model files (earth cells, air tagged 1e17) |
| `load_data2d`, `write_data2d` | ModEM Full_Impedance data (ZXY = TE, ZYX = TM), `[mV/km]/[nT]` |
| `EstimateStrike2D(data)` | Phase tensor strike of a survey, with consistency and skew |
| `StrikeData2D(data, strike)`, `RotateData2D(data, strike)` | Data rotated to the strike frame (`nothing` = auto) |
| `RotateToStrike2D(path; strike)` | Write the rotated data as `<stem>-r<ext>` |
| `ReadFwdCtrl2D`, `WriteFwdCtrl2D`, `ReadInvCtrl2D`, `WriteInvCtrl2D`, `ReadCov2D`, `WriteCov2D` | Control and covariance files |
| `ReadVFSACtrl2D`, `WriteVFSACtrl2D`, `ReadMask2D`, `WriteMask2D` | VFSA control and `mask.ctrl` |
| `Mesh2DFromInputs(model, data, fwd)` | Solver mesh and resistivity, air from `fwd.ctrl`, stations snapped to the ground |
| `BuildMesh2D(; ...)` | Padded profile mesh |
| `build_default_mt2d_mesh()` | Small benchmark mesh of the tests |
| `build_mt2d_halfspace_model(mesh)` | Uniform earth with the mesh's air |
| `mt2d_geometric_layers(f; ...)` | MakeMesh3D-style layers from skin depths |
| `mt2d_skin_depth_layers(f; ...)` | Skin-depth core layers with geometric padding below |
| `mt2d_air_layers(n, thickness, growth)` | Air layer thicknesses of `fwd.ctrl` |
| `mt2d_skin_depth(ρ, f)` | Skin depth in metres |
| `mt2d_air_mask(mesh)`, `mt2d_topo_air(mesh)` | Air cells, topographic air per column |
| `mt2d_receiver_depths(mesh)`, `mt2d_receiver_columns(mesh)`, `mt2d_station_offsets(mesh, data)` | Station ground depths, columns and snapping |
| `ReadTopo2D`, `WriteTopo2D` | `topo.dat` |
| `mt2d_profile_topography(topo, data)` | Topography projected onto the profile |
| `Topography2D(model, data, topo; water)` | Model with topography and water, mask, data Z |
| `Mask2D(model; water, fixed_below_m, fixed)` | Mask of `cov.ctrl` and `mask.ctrl` |
| `mt2d_ground(model)` | Ground depth of each column |
| `MakeMesh2D(data; ...)` | Inversion inputs from a data file (batch or GUI) |
| `run_mt2d_forward(mesh, ρ)`, `Forward2D(mesh, ρ)` | TE/TM forward solve |
| `data_from_response2d(response; ...)` | `DataFile2D` of a forward response, with errors |
| `build_mt2d_data_template(mesh; ...)`, `write_mt2d_data_template(path, mesh; ...)` | Empty data file of a mesh's stations and frequencies |
| `ForwardSolve2D(model, data, fwd)` | File forward run, writes `data.pred` |
| `chi2_rms2d(obs, pred)` | Misfit |
| `Invert2D(start, data, fwd, inv, cov, prior)` | Six-file GN or NLCG inversion |
| `Invert2D(mesh, initial, observed; algorithm, options, ...)` | In-memory deterministic inversion |
| `GaussNewton2D(...)`, `NLCG2D(...)` | `Invert2D` shorthands |
| `inv2d_frechet(problem, state)` | Data-weighted Fréchet derivative C_D^{-1/2} G in log10 resistivity |
| `inv2d_gradient(problem, state[, G])` | Objective gradient (explicit or adjoint) |
| `FrechetDerivative2D(mesh, ρ; ...)` | Explicit Fréchet derivative G = ∂g/∂m |
| `ApplyFrechet2D(mesh, ρ, δm; ...)` | δd = G δm |
| `ApplyFrechetTranspose2D(mesh, ρ, δd̂; ...)` | δm̂ = Gᵗ δd̂ |
| `WriteFrechet2D(path, mesh, ρ, data)` | G of a data file's impedances |
| `PlotInversion2D(run)` | Standard plots of a run |
| `PlotData2D(data)`, `PlotModel2D(model)` | Plots of files |
| `plot_mt2d_model`, `plot_mt2d_mesh`, `plot_mt2d_data_fit`, `plot_mt2d_site_curves`, `plot_inv2d_convergence`, `plot_vfsa2d_convergence` | Plot building blocks |

## VFSA Inversion

| Function | Description |
|:---------|:------------|
| `VFSA2D(start, data, fwd, vfsa, mask)` | Five-file 2D VFSA: no covariance, no prior |
| `VFSA2D(mesh, ρ0, observed; config, ...)` | In-memory 2D (and 1D) VFSA with threaded chains and ensemble |
| `mt2d_ensemble(models)` | Cell-wise log10 mean, median, std, p05, p95 |
| `AnalyseEnsemble2D(run_dir)` | Recompute the ensemble of a run |
| `VFSA3DMT(model; dobs_path, cfg)` | Run the 3D VFSA inversion workflow |
| `AnalyseEnsemble3D(dir)` | Compute 3D ensemble mean/median/std |
| `core_statistics(cores)` | Element-wise mean, median, std over 3D cubes |

## 3D Functions

| Function | Description |
|:---------|:------------|
| `load_data_modem(path)` | Load 3D data file, with its rotation history |
| `rotate_data(path_or_data, angle; kind)` | Rotate Z and tipper for a mesh, a strike or declination, see [rotation](data/rotation.md) |
| `RotationStep` | One step of a data file's rotation history |
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

## 3D Meshes

| Function | Description |
|:---------|:------------|
| `MakeMesh3D(data; ...)` | Start model and ModEM covariance from a data file (batch or GUI), see [3D meshes](data/mesh3d.md) |
| `MeshToMesh(model, data, target)` | Resample a model onto another mesh, with matching data and covariance |
| `extract_topography`, `write_topography`, `read_topography` | Ground surface of a 3D model |
| `extract_bathymetry`, `write_bathymetry`, `read_bathymetry` | Sea floor of a 3D model |
| `air_mask_from_topography`, `air_mask_from_model` | Air cells of a 3D model |
| `water_mask_from_bathymetry`, `water_mask_from_model` | Water cells of a 3D model |
| `chi2_and_rms_distorted(obs, pred; damping)` | 3D misfit after per-site galvanic distortion correction |

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

# Code layout

This section is for people who want to change or extend MTGeophysics.jl. It collects the details
that the user pages leave out: every control file key, how the 2D solver treats topography, the
equations of the 2D inversions, and how to add a new inversion algorithm.

For how to contribute (branches, tests, pull requests) see
[CONTRIBUTING.md](https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/main/CONTRIBUTING.md).

## 1D

| File | Contents |
|:-----|:---------|
| `src/Fwd1D.jl` | layer recursion, `mt1d_layers`, `MakeMesh1D`, Fréchet derivatives, `ForwardSolve1D` |
| `src/Inv1D.jl` | `InvCtrl1D`, `Invert1D` (its own GN and VFSA), 1D plots |

1D is self-contained: it shares only the data and model file formats and the plot style with 2D, so
changes to the 2D solver or inversions cannot change 1D results.

## 2D

| File | Contents |
|:-----|:---------|
| `src/Control2D.jl` | `fwd.ctrl`, `inv.ctrl`, VFSA control, `cov.ctrl` and `mask.ctrl` readers and writers |
| `src/Mesh2D.jl` | mesh geometry, topography helpers, model builders and model I/O |
| `src/Fwd2D.jl` | TE/TM solver, data I/O, misfit, Fréchet derivatives (G, G δm, Gᵗ δd̂) |
| `src/Strike2D.jl` | phase tensor strike estimate and rotation to the strike frame |
| `src/Topo2D.jl` | `topo.dat`, topography and water cut into models |
| `src/Mask2D.jl` | masks for `cov.ctrl` and `mask.ctrl` |
| `src/MakeMesh2D.jl`, `src/MakeMesh2DGUI.jl` | the 2D mesh tool, batch and GLMakie |
| `src/Inv2D.jl` | problem setup, objective, regularisation, line search, driver, six-file `Invert2D` |
| `src/Inv2D_GN.jl` | damped Gauss–Newton |
| `src/Inv2D_NLCG.jl` | preconditioned nonlinear conjugate gradients |
| `src/Inv2D_VFSA.jl` | very fast simulated annealing, its ensemble and the five-file `VFSA2D` |
| `src/PlotModel2D.jl` | resistivity sections and the mesh |
| `src/PlotData2D.jl` | data curves and fit, convergence, `PlotInversion2D` |

The continuous problem the 2D solver discretises is stated at the top of `src/Fwd2D.jl`.

## 3D and shared

| File | Contents |
|:-----|:---------|
| `src/Data.jl`, `src/Model.jl` | ModEM data and model I/O |
| `src/WS3DModel.jl` | WS3D model I/O |
| `src/Rotate.jl` | `rotate_data` and the rotation history |
| `src/Chi2RMS.jl`, `src/Distortion.jl` | 3D misfit, with and without galvanic distortion |
| `src/MakeMesh3D.jl`, `src/MakeMesh3DGUI.jl` | the 3D mesh tool |
| `src/MeshToMesh.jl` | resampling a model onto another mesh |
| `src/Mask3D.jl` | topography and bathymetry of 3D models |
| `src/VFSA3DMT.jl` | 3D VFSA with ModEM as the forward solver, and ensemble statistics |
| `src/CoreUtils3D.jl` | core and padding detection |
| `src/PhaseTensor.jl` | phase tensors, induction vectors and their GIS export |
| `src/PlotModel3D.jl`, `src/PlotData3D.jl`, `src/EditModel3D.jl` | GLMakie viewers and editors |

## Tests

```bash
julia --project=. test/runtests.jl
```

`test/TestFrechet2D.jl` checks the 2D Fréchet derivatives against finite differences and with the
dot-product test; `test/TestTopography2D.jl` checks the solver against a published topography model
(see [2D solver notes](solver2d.md)).

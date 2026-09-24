# 2D Mesh Tool

`MakeMesh2D` writes the inversion inputs for a 2D data file, as `MakeMesh3D` does in 3D:
`model.start`, `model.prior`, `cov.ctrl` and `inv.ctrl` (GN, NLCG), `mask.ctrl` and
`vfsa.ctrl` (VFSA), `fwd.ctrl`, `Mesh.png`, and, with topography, `data.dat` with the
station depths.

```bash
julia --project=. examples/make_mesh2D.jl data.dat topo.dat out_dir          # GLMakie window
julia --project=. examples/make_mesh2D.jl data.dat topo.dat out_dir batch    # no window
```

```julia
MakeMesh2D("data.dat"; out_dir = "mesh", topo_path = "topo.dat",
           water = [(y_range = (-30e3, -18e3), level = 95.0)], mode = :batch)
```

| Setting | Default | Meaning |
|:--------|:--------|:--------|
| `cell_width_frac` | 0.5 | core cell width as a fraction of the median station spacing |
| `core_margin_cells` | 4 | core cells beyond the outer stations |
| `n_pad`, `pad_factor` | 12, 1.3 | lateral padding cells each side and their growth |
| `first_layer_div` | 5 | first layer = δ(f_max) / this |
| `vertical_factor` | 1.1 | layer growth with depth |
| `depth_mult` | 4 | model depth = this × δ(f_min) |
| `background_resistivity` | 0 | 0 = median apparent resistivity of the data |
| `air_layers`, `air_thickness`, `air_growth`, `air_resistivity` | 10, 50 km, 2, 1e9 | the air in `fwd.ctrl` |
| `dipole_length` | 100 | TM dipole averaging next to topographic steps |
| `cov_smoothing`, `n_smooth` | 0.3, 1 | smoothing values kept in the ModEM covariance layout |
| `fixed_below_m` | Inf | fix cells below this depth (mask 0) |
| `water_resistivity` | 100 | lakes and sea |
| `inv_ctrl` | `examples/ctrl/2D/InvCtrl.GN` | copied as `inv.ctrl` |
| `vfsa_ctrl` | `examples/ctrl/2D/InvCtrl.VFSA` | copied as `vfsa.ctrl` |

The core is centred on y = 0 (the data's origin), as ModEM expects. The tool prints the
mesh summary and advice: more than one station per cell, padding shorter than δ(f_min),
a first layer thicker than δ(f_max)/3, a model shallower than δ(f_min), stations snapping
by more than half a cell.

In `mode = :gui` (GLMakie with a display) each setting is a slider with a live preview of
the mesh, the ground, the stations and δ(f_min); **Save inputs** writes the files.

## Masks

```julia
mask = Mask2D(model; water = t.mask .== 9, fixed_below_m = 30e3,
              fixed = [(y_range = (-1e3, 1e3), z_range = (0.0, 500.0))])
WriteCov2D("cov.ctrl", Cov2D(sy = fill(0.3, nz), sz = 0.3, n_smooth = 1, mask))   # GN, NLCG
WriteMask2D("mask.ctrl", mask)                                                  # VFSA
y, ground = mt2d_ground(model)                     # ground depth of each column
```

Mask values: 0 = air or fixed, 9 = water (fixed, no regularization), others free. The
two files hold the same mask. VFSA reads `mask.ctrl`, because it has no covariance.

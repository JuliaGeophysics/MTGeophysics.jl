# 3D meshes

In 3D the forward problem is solved by ModEM, but the mesh still has to come from somewhere. In this
section we build a 3D mesh from a data file, move a model from one mesh to another, and pull the
topography and water out of a model.

## Build a 3D mesh

`MakeMesh3D` reads a ModEM data file and writes a start model and a matching ModEM covariance file:

```bash
julia --project=. examples/make_mesh3D.jl data.dat start.rho dem.tif
```

or, from Julia:

```julia
MakeMesh3D("data.dat"; out_model = "start.rho", topo_file = "dem.tif", topo_crs = "EPSG:26918", mode = :nogui)
```

The model (`mesh_start_model.rho` unless you name it) and the covariance (`C3.dat`) are written
next to each other. `mode = :gui` opens a window where you can change the settings and see the mesh
update.

`topo_file` is optional. It is a digital elevation model (DEM) in the coordinate system `topo_crs`.
With it, the ground surface is cut into the mesh as air cells and the station depths in the data are
set from the DEM.

| Setting | Default | Meaning |
|:--------|:--------|:--------|
| `cell_width_frac` | 0.5 | core cell width, as a fraction of the station spacing |
| `n_pad`, `pad_factor` | 12, 1.5 | number of padding cells on each side, and their growth |
| `first_layer_div` | 4 | first layer = skin depth at the highest frequency / this |
| `vertical_factor` | 1.2 | growth of the layers with depth |
| `depth_mult` | 3 | model depth = this × skin depth at the lowest frequency |
| `air_layers` | 6 | number of air layers above the highest ground |
| `cov_value` | 0.3 | smoothing written to the covariance file |

## Move a model to another mesh

Often you want to start one inversion from the result of another, on a finer or coarser mesh.
`MeshToMesh` resamples a model onto a target mesh:

```julia
MeshToMesh("vfsa_best.rho", "data.dat", "fine_mesh.rho")
```

It writes three files into a new timestamped folder next to the target mesh: the resampled model, a
copy of the data with the station depths moved onto the new ground surface, and a matching
covariance file.

The target mesh keeps its own topography and water. Only its earth cells are filled from the
source model. `method = :nearest` (the default) keeps sharp contrasts, and `method = :linear`
smooths them.

## Topography and water from a model

A model with topography stores the air as very high resistivity. You can extract the ground surface,
or the sea floor, as a table:

```julia
m = load_ws3d_model("model.rho")
topo  = extract_topography(m)
bathy = extract_bathymetry(m; water_log10 = 0.3)   # cells below 2 Ω·m at the surface are water
write_bathymetry("bathymetry.dat", bathy)
```

A bathymetry file like this tells [3D VFSA](../inversion/3d_vfsa.md) which cells are sea and should
stay fixed.

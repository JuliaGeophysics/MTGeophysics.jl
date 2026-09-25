# 1D and 2D meshes

Before we can model or invert data we need a mesh: the cells in which the resistivity is defined.
A good mesh is fine near the stations and at shallow depth, where the data see the most detail,
and coarse far away. In this section we build meshes directly from a data file.

## 1D: nothing to do

A 1D inversion builds its own layering for every site, so usually you do not need this step.
`MakeMesh1D` makes layers that start thin at the surface (a fifth of the skin depth at the highest
frequency) and grow with depth. You can call it yourself to change the settings:

```julia
meshes = MakeMesh1D(load_data2d("data.dat"); first_layer_div = 5, vertical_factor = 1.1, depth_mult = 4)
```

## Build a 2D mesh

`MakeMesh2D` reads a data file and writes everything a 2D inversion needs:

```bash
julia --project=. examples/make_mesh2D.jl data.dat topo.dat out_dir          # opens a window
julia --project=. examples/make_mesh2D.jl data.dat topo.dat out_dir batch    # no window
```

or, from Julia:

```julia
MakeMesh2D("data.dat"; out_dir = "mesh", topo_path = "topo.dat", mode = :batch)
```

In `out_dir` you will find:

- `model.start` and `model.prior`: the start and reference models
- `fwd.ctrl`, `inv.ctrl`, `cov.ctrl`: the inputs of [Gauss–Newton and NLCG](../inversion/2d_deterministic.md)
- `vfsa.ctrl`, `mask.ctrl`: the inputs of [VFSA](../inversion/2d_vfsa.md)
- the data rotated to the strike, `data-r.dat`
- `Mesh.png` and `MeshCore.png`: pictures of the whole mesh and of its core

The tool also prints a summary with advice, for example when two stations share a cell or when the
padding is shorter than the skin depth at the lowest frequency.

With `mode = :gui` every setting is a slider with a live preview. **Save inputs** writes the files.

## The main settings

| Setting | Default | Meaning |
|:--------|:--------|:--------|
| `cell_width_frac` | 0.5 | core cell width, as a fraction of the station spacing |
| `n_pad`, `pad_factor` | 12, 1.3 | number of padding cells on each side, and their growth |
| `first_layer_div` | 5 | first layer = skin depth at the highest frequency / this |
| `vertical_factor` | 1.1 | growth of the layers with depth |
| `depth_mult` | 4 | model depth = this × skin depth at the lowest frequency |
| `strike` | `nothing` | `nothing` estimates it from the data; or degrees clockwise from north |

All settings are listed in the [control file reference](../developers/control_files.md).

## Topography and water

Give `MakeMesh2D` a `topo.dat` and it cuts the ground surface into the mesh. `topo.dat` holds WGS84
`lat lon elevation` points along the line, in metres above sea level. Where there is a lake or sea,
give the elevation of its bottom and mark the water:

```julia
MakeMesh2D("data.dat"; out_dir = "mesh", topo_path = "topo.dat",
           water = [(y_range = (-5500.0, -2500.0), level = 100.0)], water_resistivity = 200.0)
```

Here a lake lies between −5.5 and −2.5 km along the profile, with its surface at 100 m. The tool
places every station on the ground of its column, and the depths in the data file (Z) are set to
match.

## Masks

The mask says which cells the inversion may change. Gauss–Newton and NLCG read it from `cov.ctrl`,
VFSA from `mask.ctrl`. Both hold the same values:

| Value | Cells |
|:------|:------|
| 0 | air, or cells you want to keep fixed |
| 9 | water: fixed, and not smoothed across |
| any other | free |

You can build a mask yourself, for example to fix everything below 30 km:

```julia
mask = Mask2D(ReadModel2D("model.start"); fixed_below_m = 30e3)
WriteMask2D("mask.ctrl", mask)
```

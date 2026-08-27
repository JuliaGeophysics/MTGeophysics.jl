# 3D Visualization

Interactive GLMakie viewers for exploring 3D resistivity models with depth slices, cross-sections, and GIS overlays.

!!! note
    Requires GLMakie and a working OpenGL environment.

Every viewer is exported by the package, so the examples below work from any
Julia session after `using MTGeophysics` — a `Pkg.add("MTGeophysics")` install is
enough, no repository checkout needed. The scripts under `examples/` are thin
command-line wrappers around the same functions and are available if you cloned
the repository.

## Example data

The model used below is the **Cascadia** example from the ModEM example
collection, which is not bundled with this package. Download it from
[ModEM-Examples](https://github.com/magnetotellurics/ModEM-Examples/tree/main/Magnetotelluric/3D_MT/Cascadia)
and place it in `examples/Cascadia/`:

```bash
git clone https://github.com/magnetotellurics/ModEM-Examples.git
cp -r ModEM-Examples/Magnetotelluric/3D_MT/Cascadia examples/Cascadia
```

The viewers need two files from that directory:

| File | Role |
|:-----|:-----|
| `cascad_half_inverse.ws` | inverted resistivity model (WS3D format) |
| `cascad_errfl5.dat` | ModEM data file, used for station positions and georeferencing |

Set the paths once and reuse them in the examples that follow:

```julia
using MTGeophysics

model_file = "examples/Cascadia/cascad_half_inverse.ws"
data_file  = "examples/Cascadia/cascad_errfl5.dat"
```

Shapefiles for the GIS overlays are not distributed either; supply your own
coastline or boundary files (for example from [Natural Earth](https://www.naturalearthdata.com/)).
Any shapefile with a `.prj` sidecar is reprojected automatically into the
chosen coordinate system.

## Full 3D slice viewer

```julia
PlotModelXYZ(model_file, data_file)
```

Combined XY/XZ/YZ slices with depth controls, padding toggle, shapefile overlays, and figure export.
The data file is required here — it supplies the georeferencing.

From a repository checkout:

```bash
julia --project=. examples/plot_model_XYZ.jl <model.ws> <data.dat>
```

![3D model viewer](assets/plot_model_3d.png)

## XY depth slices

```julia
PlotModelXY(model_file, data_file; crs = "EPSG:32610")
```

Horizontal map-view slices at selectable depths.

From a repository checkout:

```bash
julia --project=. examples/plot_model_XY_slices.jl <model.ws> <data.dat> EPSG:32610
```

![XY slices](assets/plot_xy_slices.png)

## XY depth slices with shapefile overlays

The same viewer accepts any number of shapefiles overlaid on the map view. Each
one is reprojected from its native CRS (`.prj` sidecar) into the chosen
coordinate system:

```julia
shapefiles = [
    (path = "gis/coastline.shp",        enabled = true, color = :black, alpha = 0.95, point_size = 7, line_width = 1.2),
    (path = "gis/state_boundaries.shp", enabled = true, color = :black, alpha = 0.85, point_size = 7, line_width = 1.0),
]

PlotModelXY(model_file, data_file; crs = "EPSG:32610", shapefiles = shapefiles)
```

This works in WGS 84 lat/lon or **any** projected EPSG CRS — e.g. `EPSG:3067`
(Finland), `EPSG:32610` / `EPSG:26910` (Cascadia UTM 10N), `EPSG:3005` (BC
Albers). Set `gis_output_dir` to export the slices as georeferenced rasters;
exports carry a `.prj` generated for the chosen CRS.

From a repository checkout, edit the `shapefiles` list at the top of the script:

```bash
julia --project=. examples/plot_model_XY_with_shapefiles.jl <model.ws> <data.dat> EPSG:32610
```

## XZ cross-sections

```julia
PlotModelXZ(model_file, data_file)
```

North-South vertical sections at selectable East-West positions.

From a repository checkout:

```bash
julia --project=. examples/plot_model_XZ_slices.jl <model.ws> <data.dat>
```

![XZ slices](assets/plot_xz_slices.png)

## YZ cross-sections

```julia
PlotModelYZ(model_file, data_file)
```

East-West vertical sections at selectable North-South positions.

From a repository checkout:

```bash
julia --project=. examples/plot_model_YZ_slices.jl <model.ws> <data.dat>
```

![YZ slices](assets/plot_yz_slices.png)

## Phase tensor and induction vector map

```julia
PlotPTIVMap(data_file)
```

Period-stepped map of phase tensor ellipses and induction arrows, read straight
from a ModEM data file. Ellipses are normalised so every site is the same size
and coloured by beta skew, Phimin or Phi2; arrows use the Parkinson convention
by default and so point towards conductors. Shapefile overlays can be imported
from the window, and both the figure and the shapefiles can be exported.

From a repository checkout, edit the settings at the top of the script:

```bash
julia --project=. examples/plot_PTIV_map.jl <data.dat> [EPSG:XXXX]
```

Every period can also be written to GIS without opening a window, which is what
`--gis-only` does and what to call on a headless machine:

```julia
write_ptiv_gis(data_file; crs = "EPSG:4326", output_dir = "ptiv-gis")
```

Each period becomes `*_PT_T<period>.shp` (ellipse polygons, with the invariants
as attributes) and `*_IV_T<period>.shp` (arrow polylines), plus `.prj`/`.qpj`
sidecars and a `README.txt` recording the conventions.

The invariants themselves are available for scripting:

```julia
d  = load_data_modem(data_file)
PT = phase_tensors_from_data(d)              # nf x ns, `nothing` where undefined
IV = induction_vectors_from_data(d; convention = :parkinson)
PT[1, 1].beta                                 # skew angle in degrees
```

## Coordinate systems

Pass `crs` to the 2-D viewers (or as the third command-line argument to the
scripts):

| Mode | Description |
|:-----|:------------|
| `"model"` | Local model coordinates (metres) |
| `"EPSG:3067"` | Finnish national grid |
| `"EPSG:4326"` | WGS84 latitude/longitude |

`PlotModelXY` also accepts any projected `EPSG:XXXX` code, such as
`EPSG:32610`, `EPSG:26910`, or `EPSG:3005`.

## GIS overlays in the 3-D viewer

`PlotModelXYZ` takes a single overlay path:

```julia
PlotModelXYZ(model_file, data_file;
    shapefile_path = "gis/coastline.shp",
    crs            = "EPSG:3067",
)
```

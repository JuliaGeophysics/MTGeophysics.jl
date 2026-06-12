# 3D Visualization

Interactive GLMakie viewers for exploring 3D resistivity models with depth slices, cross-sections, and GIS overlays.

!!! note
    Requires GLMakie and a working OpenGL environment.

## Full 3D slice viewer

```bash
julia --project=. examples/plot_model_XYZ.jl
```

Combined XY/XZ/YZ slices with depth controls, padding toggle, shapefile overlays, and figure export.

![3D model viewer](assets/plot_model_3d.png)

## XY depth slices

```bash
julia --project=. examples/plot_model_XY_slices.jl
```

Horizontal map-view slices at selectable depths.

![XY slices](assets/plot_xy_slices.png)

## XY depth slices with shapefile overlays

```bash
julia --project=. examples/plot_model_XY_with_shapefiles.jl
```

Same viewer with any number of shapefiles overlaid on the map view; configure
the `shapefiles` list at the top of the script. Each shapefile is reprojected
from its native CRS (`.prj` sidecar) into the chosen coordinate system.

## Lat/lon map view (any projected CRS)

```bash
julia --project=. examples/plot_model_LL_with_shapefiles.jl
```

Generalized map view with shapefile overlays that works in WGS 84 lat/lon
(default) or **any** projected EPSG CRS — e.g. `EPSG:3067` (Finland),
`EPSG:32610` / `EPSG:26910` (Cascadia UTM 10N), `EPSG:3005` (BC Albers).
GIS exports carry a `.prj` generated for the chosen CRS.

## XZ cross-sections

```bash
julia --project=. examples/plot_model_XZ_slices.jl
```

North-South vertical sections at selectable East-West positions.

![XZ slices](assets/plot_xz_slices.png)

## YZ cross-sections

```bash
julia --project=. examples/plot_model_YZ_slices.jl
```

East-West vertical sections at selectable North-South positions.

![YZ slices](assets/plot_yz_slices.png)

## Coordinate systems

Switch the coordinate mode at the top of each viewer script:

| Mode | Description |
|:-----|:------------|
| `"model"` | Local model coordinates (metres) |
| `"EPSG:3067"` | Finnish national grid |
| `"EPSG:4326"` | WGS84 latitude/longitude |

## GIS overlays

Add shapefile overlays by setting `shapefile_path` in the viewer script:

```julia
shapefile_path = "path/to/coastline.shp"
target_crs     = "EPSG:3067"
```

## Example data

The Cascadia 3D example is not bundled. Download it from [ModEM-Examples](https://github.com/magnetotellurics/ModEM-Examples/tree/main/Magnetotelluric/3D_MT/Cascadia) and place it in `examples/Cascadia/`.

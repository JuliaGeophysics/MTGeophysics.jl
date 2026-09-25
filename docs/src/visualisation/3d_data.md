# 3D data maps

Before inverting 3D data it helps to look at them on a map. Two quantities are especially useful:

- The **phase tensor** is drawn as an ellipse at each site. It is not affected by galvanic
  distortion, and its shape and orientation show the dimensionality and strike of the structure
  under the site.
- The **induction vector** is drawn as an arrow computed from the tipper. By default the arrows
  point towards conductors.

## Plot a phase tensor and induction vector map

```julia
using MTGeophysics

PlotPTIVMap("data.dat")
```

The window steps through the periods. Ellipses have the same size at every site and are coloured by
β skew, Φmin or Φ2. You can add shapefiles from the window and export the figure.

From the command line:

```bash
julia --project=. examples/plot_PTIV_map.jl <data.dat> [EPSG:XXXX]
```

## Export to GIS

To write every period as shapefiles, without opening a window (for example on a cluster):

```julia
write_ptiv_gis("data.dat"; crs = "EPSG:4326", output_dir = "ptiv-gis")
```

Each period gives one shapefile of ellipses (`*_PT_T<period>.shp`) and one of arrows
(`*_IV_T<period>.shp`). A `README.txt` in the folder records the conventions used.

## Use the values in a script

```julia
d  = load_data_modem("data.dat")
PT = phase_tensors_from_data(d)                       # periods × sites
IV = induction_vectors_from_data(d; convention = :parkinson)
PT[1, 1].beta                                         # skew angle, in degrees
```

Entries are `nothing` where a value cannot be computed, for example when a component is missing.

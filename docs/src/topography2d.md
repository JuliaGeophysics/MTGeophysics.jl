# Topography

Topography follows the 3D code (`MakeMesh3D`, `load_ws3d_model`).

- **Model**: cells above the ground are air, written as 1e17 Ω·m; any cell above
  1e15 Ω·m reads as air. The model top is the datum. The solver gives air cells the air
  resistivity of `fwd.ctrl`.
- **Data**: Z is the depth below the model top, as in ModEM 3D. The datum is the highest
  station, so that station has Z = 0.
- **Stations** sit on the ground of their model column. A station whose Z is more than half
  a surface cell away from that ground is moved there with a warning, and the offsets are
  listed in `Summary.txt`.
- **Masks** (`cov.ctrl` for GN and NLCG, `mask.ctrl` for VFSA): air has mask 0 and water mask 9 (fixed, e.g. lakes). Neither is
  inverted nor regularized. Water may lie between stations (lakes) or in the padding (sea),
  never under a station.
- **`topo.dat`**: WGS84 `lat lon elevation` points along the line (m a.s.l.), with the lake
  and sea bottom (bathymetry) where there is water. They are projected onto the station
  polyline, so curved lines stay curved.

## Building a model with topography

```julia
topo = ReadTopo2D("topo.dat")
t = Topography2D(ReadModel2D("model.rho"), load_data2d("data.dat"), topo;
                 water = [(y_range = (-5500.0, -2500.0), level = 100.0)], water_resistivity = 200.0)
# t.model (air 1e17, water), t.mask (0 air, 9 water, 1 free), t.data (Z set), t.datum, t.ground
```

In station columns the ground moves to the cell boundary nearest the station.
[`MakeMesh2D`](mesh2d.md) does all of this from a data file and `topo.dat`.

## Stations on a staircase

Finite differences describe the ground as a staircase. Right at a step corner the TM
electric field is singular: it vanishes at a convex corner and diverges at a concave one.
A station 1 m from a corner can read anything. Next to a step, TM Ey is therefore averaged
over the tread centres of the columns within half the `Dipole length (m)` each side
(fwd.ctrl, default 100 m; the station's own column alone for a dipole shorter than a cell),
as a real electric dipole integrates E. On flat ground TM Ey is the tread-centre Ey of the
cells, each from its own ρ, interpolated between cell centres: a nodal Ey averages both
cells beside the node, and a station next to a lateral contrast such as a lake shore would
read part of its neighbour. H and TE Ex are sampled at the station itself.

Against a 4 × 4 refined mesh, the 2D-IV inversion mesh (500 m × 32 m surface cells) is
within TE 0.34 / TM 2.1 of the reference (rms in 5 % errors), and within TM 0.7 of the
2 × 2 refined data mesh. The largest misfits are at the lake-shore station and next to a
step; longer dipoles (250–1000 m) were tried and made it worse.

## Trapezoidal hill

`test/TestTopography2D.jl` solves the hill of Wannamaker, Stodt & Rijo (1986) on a 50 m
staircase mesh. The hill is 450 m high, 2 km wide at its base and 450 m at its flat top,
over a 100 Ω·m halfspace, at 2 Hz. The test checks the published behaviour: TM ρa drops
below 40 Ω·m on the hilltop and peaks near the foot, TE is mildly raised on top, and both
return to 100 Ω·m away from the hill.

## Fréchet derivatives

G and Gᵗ include the topography: the receiver sampling reads every surface row it needs,
and air and water cells have zero columns. `test/TestTopography2D.jl` checks them against
finite differences and with the dot-product test.

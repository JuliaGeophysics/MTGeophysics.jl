# 1D Modelling and Inversion

1D reads the same data file as 2D and 3D, and every site is a sounding inverted on its
own. The only other input is one small control file. `MakeMesh1D` builds each site's
layering from its skin depths, and the site starts from its median apparent resistivity.
A 1D model is solved exactly by the layered-earth impedance recursion, not by finite
differences.

| File | Contents |
|:-----|:---------|
| `src/Fwd1D.jl` | layer recursion, `Mesh1D`, `MakeMesh1D`, Fréchet derivatives, `ForwardSolve1D` |
| `src/Inv1D.jl` | `InvCtrl1D`, `Invert1D` (GN or VFSA), 1D model plots |

## Files

- **Data**: the ModEM Full_Impedance file of 2D and 3D.
- **Model** (forward runs and results): the 2D model layout with one column: `1 nz LOGE`,
  a dummy width, the layer thicknesses from the surface down, then ln ρ per layer. The last
  layer is the halfspace.
- **Control** (`examples/ctrl/1D/InvCtrl.GN`, `InvCtrl.VFSA`):

  ```text
  Algorithm                         : GN
  Mode                              : XYYX
  Exit search when rms is less than : 1
  Maximum number of iterations      : 20
  Log10 resistivity bounds          : 0 4
  Initial damping factor lambda     : 1
  ```

  `Mode` picks the impedance that is fitted: `XY` (Zxy), `YX` (Zyx = −Z), `XYYX` (both),
  or `DET`, the determinant impedance √det Z. A VFSA file replaces the lambda line with
  `Number of chains` and `Control points`. Everything else is fixed inside the code:
  - GN: vertical smoothing with smallness 0.01, damping 0.01, and a step of at most 0.5
    in log10 ρ;
  - VFSA: temperature 1, cooling ratio 0.001, step scale 0.11, and an RBF width of 2.5
    layers.

## Benchmark and examples

```bash
julia --project=. helpers/benchmarks_1D.jl         # writes examples/data/1D-I (data.dat, model.true)
julia --project=. examples/run_fwd1D.jl            # model.true -> data.pred and a data plot
julia --project=. examples/run_inv1D.jl GN         # or VFSA
julia --project=. examples/run_inv1D.jl data.dat InvCtrl
```

1D-I is a five-layer earth (100, 20, 350, 40, 800 Ω·m; 120, 280, 650, 1400 m) under one
site near Jyväskylä, 23 frequencies from 0.003 to 1000 Hz, 5 % noise. The mesh settings
(first layer δ(f_max)/5, growth 1.1, depth 4 δ(f_min)) are constants at the top of
`examples/run_inv1D.jl`.

## From Julia

```julia
using MTGeophysics

pred = ForwardSolve1D("model.rho", "data.dat"; mode = :XYYX)       # writes data.pred

meshes = MakeMesh1D(load_data2d("data.dat"); first_layer_div = 5, vertical_factor = 1.1, depth_mult = 4)
run = Invert1D("data.dat", "InvCtrl.GN", meshes)
PlotInversion1D(run; true_model_path = "model.true")
```

`Invert1D` inverts every site (or those named in `sites`) on its own. It writes one folder
per site (`model.start`, `model.rho`, `data.pred`, `History.csv` or `vfsa/`, `plots/`),
plus `data.pred` of all sites, `Summary.txt` and the inputs. Without `meshes` it calls
`MakeMesh1D` with its defaults.

In memory, a 1D mesh is an `MT2DMesh` with one column:

```julia
mesh = Mesh1D(thicknesses, frequencies)            # last layer = halfspace
ρ = fill(100.0, length(thicknesses), 1)
response = run_mt2d_forward(mesh, ρ)               # exact, z_yx = -z_xy
G = FrechetDerivative2D(mesh, ρ; parameterization = :log10_resistivity)
result = Invert2D(mesh, ρ, observed_site; algorithm = GaussNewton2DConfig())
Z = mt1d_impedance(frequencies, resistivities, thicknesses)
```

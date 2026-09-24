# 1D Modelling and Inversion

1D reads the same data file as 2D and 3D, and every site is a sounding inverted on its
own. The only other input is one small control file. `MakeMesh1D` builds each site's
layering from its skin depths, and the site starts from its median apparent resistivity.
A 1D model is solved exactly by the layered-earth impedance recursion, not by finite
differences. 1D is self-contained: it shares only the data and model file formats and the
plot style with 2D, so changes to the 2D solver or inversions cannot change 1D results.

| File | Contents |
|:-----|:---------|
| `src/Fwd1D.jl` | layer recursion, `mt1d_layers`, `MakeMesh1D`, Fréchet derivatives, `ForwardSolve1D` |
| `src/Inv1D.jl` | `InvCtrl1D`, `Invert1D` (its own GN and VFSA), 1D plots |

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
  Initial damping factor lambda     : 1
  ```

  `Mode` picks the impedance that is fitted: `XY` (Zxy), `YX` (Zyx = −Z), `XYYX` (both),
  or `DET`, the determinant impedance √det Z. A VFSA file replaces the lambda line with
  `Log10 resistivity bounds` (the search box) and `Number of chains`; GN is unbounded. Everything else is fixed inside the code:
  - GN: vertical smoothing with smallness 0.01, Levenberg–Marquardt damping 0.01, and a
    step of at most 0.5 in log10 ρ;
  - VFSA: every layer is a parameter (no sparse parameterisation), a fifth of them moved
    per proposal; temperature 0.03, cooling ratio 0.001, step scale 0.11. The chains' best
    models form the ensemble (mean, median, 5–95 %). `model.rho` is the ensemble mean, as
    in 2D. Unregularised per-layer models are rough and equivalent layers trade off, so the
    mean fits worse than the chains do: `Summary.txt` gives both, and `vfsa/data.best.pred`
    holds the best chain's response.

## Benchmark and examples

```bash
julia --project=. helpers/benchmarks_1D.jl         # writes examples/data/1D-I (data.dat, model.true)
julia --project=. examples/run_fwd1D.jl            # model.true -> data.pred and a data plot
julia --project=. examples/run_inv1D.jl GN         # or VFSA
julia --project=. examples/run_inv1D.jl data.dat InvCtrl
```

1D-I is a five-layer earth (100, 20, 350, 40, 800 Ω·m; 120, 280, 650, 1400 m) under one
site (`Fin001`) near Jyväskylä, 23 frequencies from 0.003 to 1000 Hz, 5 % noise. The mesh settings
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

In memory a 1D model is its thicknesses and resistivities (the last layer the halfspace):

```julia
h = mt1d_layers(frequencies; background_resistivity = 100, first_layer_div = 5)
ρ = fill(100.0, length(h))
Z = mt1d_impedance(frequencies, ρ, h)              # Zxy; Zyx = -Z
G = mt1d_frechet(frequencies, ρ, h)                # ∂Z/∂log10 ρ, nf × layers
```

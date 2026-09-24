# 2D Forward Modelling

TE and TM responses of a 2D profile by finite differences (box integration), with the
continuous problem stated at the top of `src/Fwd2D.jl`. Inputs and outputs follow ModEM:
a model, a data file and a forward control file in, `data.pred` out.

## Files

- **Model** (`model.rho`): the ModEM-style layout, earth cells only, after one `#` line:
  `ny nz LOGE`, the cell widths, the layer thicknesses, `0`, then ln ρ, one layer per
  block. The grid is centred on the data's y = 0. Topographic air cells hold 1e17 Ω·m
  (see [Topography](topography2d.md)).
- **Data** (`data.dat`): ModEM Full_Impedance, ZXY = TE and ZYX = TM, `exp(+iωt)`,
  `[mV/km]/[nT]`, the same file 1D and 3D read. Y is the position along the profile,
  Z the depth below the model top (0 on flat ground).
- **`fwd.ctrl`** (`examples/ctrl/2D/FwdCtrl`), required, it defines the air:

  ```text
  Mode                     : TETM
  Air layers               : 10
  Air thickness (m)        : 50000
  Air growth factor        : 2
  Air resistivity (ohm m)  : 1e+09
  Write Frechet derivative : no
  Dipole length (m)        : 100
  ```

  `Write Frechet derivative : yes` also writes G = ∂d/∂m as `data.frechet`.
  `Dipole length` only matters next to topographic steps.

## Running

```bash
julia --project=. helpers/benchmarks_2D.jl                     # examples/data/2D-IV
julia --project=. examples/run_fwd2D.jl                        # model.true -> data.pred
julia --project=. examples/run_fwd2D.jl model.rho data.dat FwdCtrl
```

```julia
using MTGeophysics
pred = ForwardSolve2D("model.rho", "data.dat", "FwdCtrl")     # data.pred next to data.dat
PlotData2D("data.dat"; predicted_path = pred, output_path = "DataFit.png")
PlotModel2D("model.rho"; output_path = "Model.png")
```

A data file whose impedances are all zero is a template: `ForwardSolve2D` fills it and
reads its errors as fractions of |Z|.

## In memory

```julia
mesh, ρ = Mesh2DFromInputs(ReadModel2D("model.rho"), load_data2d("data.dat"), ReadFwdCtrl2D("FwdCtrl"))
response = run_mt2d_forward(mesh, ρ; mode = :TETM)            # MT2DResponse, (frequency, site)
plot_mt2d_model(mesh, ρ; output_path = "model.png")
```

`BuildMesh2D(; ...)` builds a padded mesh directly; `mt2d_geometric_layers(f; ...)` gives
MakeMesh3D-style layers (first layer δ(f_max)/`first_layer_div`, growth `vertical_factor`,
down to `depth_mult` δ(f_min)).

## Benchmarks

`helpers/benchmarks_2D.jl` writes each case into `examples/data/<case>`. The data come from a
fine mesh, while the start and prior models, `cov.ctrl` and `mask.ctrl` sit on a coarser
inversion mesh, so there is no inverse crime. Stations `JYV001…` run E–W near Jyväskylä
(62.25°N, 25.75°E). 2D-IV is the default.

For 2D-IV the data mesh is the inversion mesh with every cell split 2 × 2, and each fine
cell takes its parent's air or water. The two meshes therefore share one topography
staircase and differ only below the ground. Cutting the relief into each mesh on its own
gives two different staircases, and the TM response next to a step depends strongly on
where the step lies: the two meshes then disagree at a level of rms ~100, which no
inversion can fit.

| Case | Description |
|:-----|:------------|
| 2D-I | Thin vertical conductive dyke in a two-layer background |
| 2D-II | Two resistive blocks in a two-layer background |
| 2D-III | Mixed conductors and resistors in a three-layer background |
| 2D-IV | 2D-III under 100–200 m of smooth relief, with a lake in the padding and `topo.dat` |

```bash
julia --project=. helpers/benchmarks_2D.jl 2D-I 2D-II 2D-III 2D-IV
```


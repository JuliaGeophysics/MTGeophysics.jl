# Data and model files

MTGeophysics.jl reads and writes the same files as ModEM. This means you can prepare data with your
usual tools, run a model here, and continue in ModEM, or the other way around. In this section we
look at the files you will meet most often.

## Data files

1D, 2D and 3D all use one data format: the ModEM `Full_Impedance` file. Impedances are in
`[mV/km]/[nT]` with the `exp(+iωt)` sign convention. 3D files may also hold the tipper.

```julia
using MTGeophysics

d3 = load_data_modem("data.dat")   # the full file: impedance tensor and tipper
d2 = load_data2d("data.dat")       # the ZXY and ZYX parts, as 1D and 2D use them
```

In 2D, ZXY is the TE mode and ZYX the TM mode. This only holds when the data are in the strike frame,
so every 2D workflow first rotates the data to the strike (see [Rotation and strike](rotation.md)).

!!! note
    A data file whose impedances are all zero works as a template. A forward run fills it with the
    predicted response for those stations and frequencies.

## Model files

| Dimension | Format | Read | Write |
|:----------|:-------|:-----|:------|
| 1D | the 2D layout with a single column | `ReadModel2D` | `WriteModel2D` |
| 2D | ModEM 2D layout, ln ρ | `ReadModel2D` | `WriteModel2D` |
| 3D | WS3D / ModEM 3D | `load_ws3d_model`, `load_model_modem` | `write_ws3d_model`, `write_model_modem` |

Cells above the ground are air. They are written as 10¹⁷ Ω·m, and any cell above 10¹⁵ Ω·m is read
back as air. This is how topography is stored in a model file.

## Control files

The 1D and 2D workflows are driven by short text files with one `key : value` per line, for example:

```text
Mode                     : TETM
Strike (deg)             : auto
```

| File | Used by | What it sets |
|:-----|:--------|:-------------|
| `fwd.ctrl` | 2D forward and all 2D inversions | mode, strike and the air layers |
| `inv.ctrl` | 2D Gauss–Newton and NLCG | algorithm, regularisation, stopping rules |
| `vfsa.ctrl` | 2D VFSA | search bounds, chains, cooling |
| `cov.ctrl` | 2D Gauss–Newton and NLCG | which cells are free, fixed, air or water |
| `mask.ctrl` | 2D VFSA | the same mask as `cov.ctrl`, without a covariance |
| `InvCtrl.GN`, `InvCtrl.VFSA` | 1D inversion | algorithm and stopping rules |
| `topo.dat` | 2D mesh tool | ground elevation along the profile |

Ready-to-use examples are in `examples/ctrl/1D` and `examples/ctrl/2D`, and the
[2D mesh tool](mesh2d.md) writes a full set for new data. Every key is listed in the
[control file reference](../developers/control_files.md).

## Misfit

The misfit between observed and predicted data is reported as RMS: the square root of χ² divided by
the number of real data. The real and imaginary parts of each impedance count as two data.

```julia
fit = chi2_rms2d(load_data2d("data.dat"), load_data2d("data.pred"))   # 1D and 2D
fit = chi2_and_rms("data.dat", "pred.dat")                            # 3D
fit.rms
```

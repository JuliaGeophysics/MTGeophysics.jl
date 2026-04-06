# MTGeophysics.jl — Examples and workflows

This file is the command-oriented companion to the main repository README.
Use it when you already know what you want to run and just need the exact
commands.

All commands assume you are in the repository root:

```powershell
cd MTGeophysics.jl
```

## 0. One-time setup

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using MTGeophysics; println("OK")'
```

Optional validation:

```powershell
julia --project=. test/runtests.jl
```

## 1. External Cascadia example data

The 3-D Cascadia example is not bundled in this repository. Download the
`Cascadia` directory from:

https://github.com/magnetotellurics/ModEM-Examples/tree/main/Magnetotelluric/3D_MT/Cascadia

Place it here:

```text
examples/Cascadia/
```

The 3-D viewer scripts currently default to:

- `examples/Cascadia/cascad_half_inverse.ws`
- `examples/Cascadia/cascad_errfl5.dat`

## 2. Generate synthetic benchmarks

```powershell
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl
```

Generated output directories:

- `examples/0Layered1D/`
- `examples/0COMEMI2D-I/`
- `examples/0COMEMI2D-II/`
- `examples/0COMEMI2D-III/`

## 3. Run the forward examples

### 1-D

```powershell
julia --project=. .\examples\response_1d.jl
```

Custom input files:

```powershell
julia --project=. .\examples\response_1d.jl examples/0Layered1D/Layered1D.true examples/0Layered1D/Layered1D.ref
```

### 2-D

```powershell
julia --project=. .\examples\response_2d.jl
```

Custom input files:

```powershell
julia --project=. .\examples\response_2d.jl examples/0COMEMI2D-I/Comemi2D1.true examples/0COMEMI2D-I/Comemi2D1.ref
```

## 4. Write standard plots from the library API

### 1-D model plot

```powershell
julia --project=. -e 'using MTGeophysics; PlotModel1D("examples/0Layered1D/Layered1D.true", "examples/0Layered1D/Layered1D.ref"; output_path="examples/0Layered1D/ModelPlot1D.png")'
```

### 1-D data plot

```powershell
julia --project=. -e 'using MTGeophysics; PlotData1D("examples/0Layered1D/Layered1D.obs"; output_path="examples/0Layered1D/DataPlot1D.png")'
```

### 2-D model plot

```powershell
julia --project=. -e 'using MTGeophysics; PlotModel2D("examples/0COMEMI2D-I/Comemi2D1.true"; output_path="examples/0COMEMI2D-I/ModelPlot2D.png")'
```

### 2-D data plots

```powershell
julia --project=. -e 'using MTGeophysics; PlotData2D("examples/0COMEMI2D-I/Comemi2D1.obs"; maps_output_path="examples/0COMEMI2D-I/DataMaps2D.png", curves_output_path="examples/0COMEMI2D-I/DataCurves2D.png")'
```

## 5. Run the 2-D VFSA inversion

### Bundled example script

```powershell
julia --project=. .\examples\run_vfsa2dmt.jl
```

### CLI entry point

```powershell
julia --project=. -e 'using MTGeophysics; MTGeophysics.main_vfsa2dmt()' -- examples/0COMEMI2D-I/Comemi2D1.ini examples/0COMEMI2D-I/Comemi2D1.obs --n-chains 3 --n-ctrl 400 --max-iter 300 --log-bounds 0,4 --seed 20260308
```

Useful flags:

- `--n-chains`
- `--n-trials`
- `--max-iter`
- `--n-ctrl`
- `--log-bounds lo,hi`
- `--perturb-depth D`
- `--output-root DIR`
- `--seed N`

## 6. Post-process a 2-D inversion run

### Recompute ensemble statistics

```powershell
julia --project=. Helpers/run_statistics_2d.jl examples/run_VFSA2DMT_<timestamp>
```

### Build a convergence GIF

```powershell
julia --project=. Helpers/make_gif_2d.jl examples/run_VFSA2DMT_<timestamp>
```

With explicit options:

```powershell
julia --project=. Helpers/make_gif_2d.jl examples/run_VFSA2DMT_<timestamp> my_animation.gif --fps 8 --depth_km 50 --rho_range 0,4
```

## 7. Explore 3-D models interactively

### Full 3-D slice viewer

```powershell
julia --project=. .\examples\plot_model_3D.jl
```

### XY depth-slice viewer

```powershell
julia --project=. .\examples\plot_XY_slices.jl
```

### XZ cross-section viewer

```powershell
julia --project=. .\examples\plot_XZ_slices.jl
```

### YZ cross-section viewer

```powershell
julia --project=. .\examples\plot_YZ_slices.jl
```

In the slice-viewer scripts, switch the coordinate mode at the top of the file
to one of:

- `"EPSG:3067"`
- `"EPSG:4326"`
- `"model"`

## 8. Replace slice resistivity

Interactively replace resistivity below a chosen depth layer, with scope
control (core-only or full model including padding):

```powershell
julia --project=. .\examples\replace_slice_resistivity_scope.jl
```

Use the depth slider to pick a cutoff layer, enter a target resistivity and
optional blend percentage, then click **Apply Changes**. Toggle between
core-only and full-model scope with the **Show Core Model / Show Full Model**
button. Click **Save Model** to write the edited model to disk.

## 9. Experimental editing scripts

The repository also contains interactive editing scripts for manual model
modification:

- `examples/edit_model_by_slice.jl`
- `examples/edit_model_by_drawing.jl`

These are currently more development-oriented than the main workflows and may
need local path adjustment before use.

## 10. Recommended first session

If you are using the repository for the first time, this is the safest path:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl
julia --project=. .\examples\response_1d.jl
julia --project=. .\examples\response_2d.jl
julia --project=. .\examples\run_vfsa2dmt.jl
```

Then, if you have downloaded the Cascadia dataset:

```powershell
julia --project=. .\examples\plot_XY_slices.jl
```

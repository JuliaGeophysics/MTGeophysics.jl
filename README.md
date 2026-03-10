# MTGeophysics.jl

MTGeophysics.jl is a Julia repository for magnetotelluric forward modelling,
data I/O, inversion workflows, and interactive visualisation.

The repository currently combines three layers of functionality:

1. 3-D ModEM model and data I/O, plus misfit utilities.
2. 1-D and 2-D MT forward modelling, plotting, and inversion helpers.
3. Interactive 3-D slice viewers for inspecting ModEM models, including GIS
	 overlays and coordinate reprojection.

![2D ensemble mean model](images/vfsa2d_plot_model_mean.png)

![2D convergence history](images/vfsa2d_plot_convergence.png)

![2D convergence animation](images/vfsa2d_convergence.gif)

## What the repository provides

### Core package functionality

- 3-D ModEM data loading and model loading/writing.
- Chi-square and RMS evaluation for ModEM-style datasets.
- 1-D layered MT forward modelling and analytical / finite-difference solves.
- 2-D MT forward modelling, data-template generation, observed-data I/O, and
	response plotting.
- 2-D VFSA inversion with run logging, snapshots, and ensemble statistics.
- Interactive 3-D XY, XZ, and YZ slice viewers using GLMakie.
- Optional shapefile overlay support for the 3-D viewers through `Shapefile.jl`,
	`GeoInterface.jl`, and `Proj.jl`.

### Repository extras

- Benchmark generators for 1-D and 2-D synthetic examples.
- Example scripts for forward modelling and interactive viewing.
- Post-processing helpers for inversion statistics and GIF generation.
- A lightweight test suite covering core I/O and 1-D / 2-D forward modelling.

## Repository layout

```text
MTGeophysics.jl/
├── src/
│   ├── MTGeophysics.jl      # top-level module and exported API
│   ├── Data.jl              # 3D ModEM data I/O
│   ├── Model.jl             # 3D ModEM model I/O
│   ├── Chi2RMS.jl           # misfit utilities
│   ├── MTGeophysics1D.jl    # 1D modelling and plotting
│   ├── MTGeophysics2D.jl    # 2D modelling and plotting
│   ├── VFSA2DMT.jl          # 2D VFSA inversion workflow
│   └── PlotModel.jl         # 3D interactive plotting utilities
├── Helpers/
│   ├── benchmarks_1d.jl
│   ├── benchmarks_2d.jl
│   ├── run_statistics_2d.jl
│   └── make_gif_2d.jl
├── examples/
│   ├── response_1d.jl
│   ├── response_2d.jl
│   ├── run_vfsa2dmt.jl
│   ├── plot_model_3D.jl
│   ├── plot_XY_slices.jl
│   ├── plot_XZ_slices.jl
│   ├── plot_YZ_slices.jl
│   ├── Cascadia/           # external example data, downloaded separately
│   └── 0COMEMI2D-*/        # generated synthetic benchmarks
├── images/                 # documentation figures and GIFs
├── test/
└── README.md
```

## Requirements

- Julia `1.10` or newer. The current working environment in this repository is
	using Julia `1.12.4`.
- A working OpenGL environment for GLMakie-based interactive viewers.
- Internet access on first install to resolve Julia packages and registries.

On Windows, if your network uses TLS interception or an enterprise proxy, make
sure Julia package resolution is working before trying to instantiate the
project.

## Installation

From the repository root:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

To validate the installation:

```powershell
julia --project=. -e 'using MTGeophysics; println("MTGeophysics loaded")'
julia --project=. test/runtests.jl
```

## Quick start

### 1. Generate the bundled synthetic benchmarks

```powershell
julia --project=. Helpers/benchmarks_1d.jl
julia --project=. Helpers/benchmarks_2d.jl
```

This creates:

- `examples/0Layered1D/` for the 1-D layered benchmark.
- `examples/0COMEMI2D-I/`, `examples/0COMEMI2D-II/`, and
	`examples/0COMEMI2D-III/` for the 2-D synthetic examples.

### 2. Run a 1-D forward example

```powershell
julia --project=. .\examples\response_1d.jl
```

This uses the default generated benchmark files and writes:

- a predicted 1-D response file next to the benchmark,
- `DataPlot1D.png` for observed-versus-predicted comparison.

You can also drive the plotting helpers directly from Julia:

```powershell
julia --project=. -e 'using MTGeophysics; PlotModel1D("examples/0Layered1D/Layered1D.true", "examples/0Layered1D/Layered1D.ref"; output_path="examples/0Layered1D/ModelPlot1D.png")'
```

### 3. Run a 2-D forward example

```powershell
julia --project=. .\examples\response_2d.jl
```

This writes a predicted response plus:

- `DataMaps2D.png`
- `DataCurves2D.png`

To write a standalone 2-D model figure from the library API:

```powershell
julia --project=. -e 'using MTGeophysics; PlotModel2D("examples/0COMEMI2D-I/Comemi2D1.true"; output_path="examples/0COMEMI2D-I/ModelPlot2D.png")'
```

### 4. Run the 2-D VFSA inversion workflow

The bundled script uses the `0COMEMI2D-I` benchmark:

```powershell
julia --project=. .\examples\run_vfsa2dmt.jl
```

This creates a timestamped result directory under `examples/`, containing:

- run summaries,
- per-chain logs,
- snapshot models,
- ensemble mean / median / standard deviation models,
- data and model plots,
- optional convergence GIF input files.

You can also use the CLI entry point defined in `VFSA2DMT.jl`:

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

### 5. Post-process an inversion run

Recompute ensemble statistics:

```powershell
julia --project=. Helpers/run_statistics_2d.jl examples/run_VFSA2DMT_<timestamp>
```

Generate a convergence GIF from saved snapshots:

```powershell
julia --project=. Helpers/make_gif_2d.jl examples/run_VFSA2DMT_<timestamp>
```

## 3-D ModEM viewers

The repository includes several GLMakie-based interactive viewers for 3-D
models.

### Full 3-D slice viewer

```powershell
julia --project=. .\examples\plot_model_3D.jl
```

This viewer supports:

- XY / XZ / YZ slice inspection,
- depth and padding control,
- shapefile overlays,
- coordinate reprojection with `Proj.jl`,
- figure export.

### XY depth slices

```powershell
julia --project=. .\examples\plot_XY_slices.jl
```

The current default paths are:

- `examples/Cascadia/cascad_half_inverse.ws`
- `examples/Cascadia/cascad_errfl5.dat`

### XZ cross-sections

```powershell
julia --project=. .\examples\plot_XZ_slices.jl
```

### YZ cross-sections

```powershell
julia --project=. .\examples\plot_YZ_slices.jl
```

The slice viewers allow three coordinate-system modes at the top of the script:

- `"EPSG:3067"`
- `"EPSG:4326"`
- `"model"`

If you use a projected or geographic CRS, the corresponding data file is used
for georeferencing the model and station locations.

## External 3-D example data

The Cascadia 3-D example data are not bundled in this repository. Download the
`Cascadia` directory from:

https://github.com/magnetotellurics/ModEM-Examples/tree/main/Magnetotelluric/3D_MT/Cascadia

After downloading, place it here:

```text
examples/Cascadia/
```

The current 3-D viewer scripts are already configured to use:

- `examples/Cascadia/cascad_half_inverse.ws`
- `examples/Cascadia/cascad_errfl5.dat`

## Public API summary

The top-level module exports functions and types from several layers of the
codebase.

### 3-D ModEM utilities

- `load_data_modem`
- `load_model_modem`
- `write_model_modem`
- `chi2_and_rms`

### 1-D workflow

- `BuildMesh1D`, `MakeMesh1D`
- `Forward1D`, `ForwardSolve1D`
- `load_mt1d_model`, `write_mt1d_model`
- `load_mt1d_data_spec`, `load_mt1d_observed_data`
- `PlotData1D`, `PlotModel1D`

### 2-D workflow

- `BuildMesh2D`, `MakeMesh2D`
- `Forward2D`, `ForwardSolve2D`
- `load_model2d`, `write_model2d`
- `build_mt2d_data_template`, `write_mt2d_data_template`
- `load_data2d`, `write_data2d`
- `PlotData2D`, `PlotModel2D`
- `VFSA2DMT`, `VFSA2DMTConfig`, `VFSA2DMTParams`

## Example library usage

```julia
using MTGeophysics

pred_path = ForwardSolve1D("examples/0Layered1D/Layered1D.true",
													 "examples/0Layered1D/Layered1D.ref")

PlotData1D("examples/0Layered1D/Layered1D.ref", pred_path;
					 output_path = "examples/0Layered1D/DataPlot1D.png")

PlotModel2D("examples/0COMEMI2D-I/Comemi2D1.true";
						output_path = "examples/0COMEMI2D-I/ModelPlot2D.png")
```

## Experimental scripts

The repository also contains interactive model-editing scripts:

- `examples/edit_model_by_slice.jl`
- `examples/edit_model_by_drawing.jl`

These are useful development tools, but they currently contain hard-coded paths
and references that are more project-specific than the main forward / inversion
examples. Treat them as advanced or experimental utilities rather than the
recommended starting point.

## Notes and troubleshooting

### Package setup on managed networks

If Julia package downloads fail behind a proxy or TLS-inspecting network,
ensure your Julia environment can resolve registries and packages before first
use.

### GLMakie viewer issues

If an interactive viewer fails to open:

- confirm that GLMakie precompiled successfully,
- update graphics drivers,
- test a simple `using GLMakie` session first,
- try one of the non-interactive plotting helpers to isolate whether the issue
	is OpenGL-specific.

### Case-sensitive paths in the docs

The repository uses lowercase `examples/`, `Helpers/`, and `test/` directory
names on disk. Commands in this README follow the actual current layout.

## Additional guide

For a shorter command-oriented walkthrough, see `examples/0Example.md`.


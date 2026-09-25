# 2D deterministic inversion

A deterministic inversion starts from a model and improves it step by step until the predicted data
fit the observed data. Many models can fit MT data, so we also ask for the model to be smooth and
close to a reference model. The inversion balances these two goals.

MTGeophysics.jl has two deterministic algorithms:

- **Gauss–Newton (GN)** takes few, large steps. Each step is expensive, so it suits small and
  medium models.
- **Nonlinear conjugate gradients (NLCG)** takes many cheap steps and uses little memory, so it
  suits large models.

Both read the same six files, as ModEM does.

## What you need

| File | Contents |
|:-----|:---------|
| `model.start` | the model to start from |
| `data.dat` | the observed data |
| `fwd.ctrl` | mode, strike and air (see [2D forward](../forward/2d.md)) |
| `inv.ctrl` | the algorithm and its settings |
| `cov.ctrl` | the mask: which cells are free, fixed, air or water |
| `model.prior` | the reference model for the regularisation |

You do not have to write these by hand. The [2D mesh tool](../data/mesh2d.md) writes all of them
for a new data file, and the benchmark folders already contain them.

## Run the example

```bash
julia --project=. helpers/benchmarks_2D.jl        # writes examples/data/2D-IV
julia --project=. examples/run_inv2D.jl GN        # or NLCG
```

To use your own files, pass all six:

```bash
julia --project=. examples/run_inv2D.jl model.start data.dat FwdCtrl InvCtrl cov.ctrl model.prior
```

## The inversion control file

The example `examples/ctrl/2D/InvCtrl.GN`:

```text
Algorithm                         : GN
Initial damping factor lambda     : 1
Exit search when rms is less than : 1
Maximum number of iterations      : 20
Mode                              : TETM
Max log10 step                    : 0.5
Max line search steps             : 12
Smallness weight                  : 0.01
Smoothing weight y                : 1
Smoothing weight z                : 1
GN damping                        : 0.01
```

The lines you will change most often are:

- `Algorithm`: `GN` or `NLCG`.
- `lambda`: the weight of the regularisation. Larger values give smoother models that fit the data
  less closely.
- `Exit search when rms is less than`: the target RMS. An RMS of 1 means the data are fitted to
  within their errors.
- `Mode`: invert `TE`, `TM` or both (`TETM`).

The other keys are explained in the [control file reference](../developers/control_files.md).

## From Julia

```julia
using MTGeophysics

run = Invert2D("model.start", "data.dat", "FwdCtrl", "InvCtrl.GN", "cov.ctrl", "model.prior")
PlotInversion2D(run; true_model_path = "model.true")
```

## What you get

Everything goes into a new folder, `run_YYYYmmdd_HHMMSS/`, next to the data:

```text
run_YYYYmmdd_HHMMSS/
├── model.rho       final model (can be used as a new start model)
├── data.pred       its predicted data
├── History.csv     rms and objective at every iteration
├── Summary.txt     final rms, stopping reason, strike
├── inputs/         copies of the six input files
└── plots/          mesh, start, final and true models, data fit, convergence
```

If the data were rotated to the strike, the rotated data are also saved, as `data-r.dat`.

!!! tip "Is the inversion done?"
    Look at `plots/` first. A good result fits the data (RMS close to your target) with a model
    that is no rougher than it needs to be. If the RMS stays high, check the data fit plot for
    single stations or periods that are fitted badly.

The equations behind both algorithms are described in [2D inversion theory](../developers/inversion2d.md).

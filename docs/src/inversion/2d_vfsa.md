# 2D VFSA inversion

Very fast simulated annealing (VFSA) is a global search. Instead of following the gradient downhill,
it tries random changes to the model. A change that improves the fit is kept. A change that makes
the fit worse is sometimes kept too, so the search can climb out of a poor local solution. As the
run goes on, the "temperature" drops: the changes become smaller and bad changes are accepted less
often.

Because each run is random, we run several independent chains. Their best models form an
**ensemble**, and the spread of the ensemble shows how well each part of the model is resolved.

To keep the search small, VFSA does not change every cell. It moves a few hundred control points,
and the model between them follows smoothly (Gaussian radial basis functions).

## What you need

VFSA reads five files. There is no covariance and no prior model, because VFSA has no
regularisation: the start model is the centre of the search, and the mask says which cells may
change.

| File | Contents |
|:-----|:---------|
| `model.start` | the model to start from |
| `data.dat` | the observed data |
| `fwd.ctrl` | mode, strike and air (see [2D forward](../forward/2d.md)) |
| `vfsa.ctrl` | the VFSA settings |
| `mask.ctrl` | which cells are free, fixed, air or water (see [Masks](../data/mesh2d.md)) |

The [2D mesh tool](../data/mesh2d.md) writes all five.

## Run the example

Each chain runs on its own Julia thread, so start Julia with as many threads as chains (`-t`):

```bash
julia --project=. -t 10 examples/run_vfsa2D.jl        # examples/data/2D-IV
julia --project=. -t 10 examples/run_vfsa2D.jl model.start data.dat FwdCtrl InvCtrl.VFSA mask.ctrl
```

## The VFSA control file

The example `examples/ctrl/2D/InvCtrl.VFSA` starts like this:

```text
Exit search when rms is less than : 1
Maximum number of iterations      : 3000
Mode                              : TETM
Log10 resistivity bounds          : 0 4
Number of chains                  : 2
Control points                    : 400
...
```

The settings that matter most are:

- `Maximum number of iterations`: iterations per chain. More iterations give a slower, more careful
  search.
- `Log10 resistivity bounds`: the range the search may explore, here 1 to 10 000 Ω·m.
- `Number of chains`: the size of the ensemble.
- `Control points`: how many points the search moves. More points allow more detail but make the
  search harder.

The remaining keys (cooling, control-point widths, core depth, padding) have sensible defaults and are
explained in the [control file reference](../developers/control_files.md).

## From Julia

```julia
using MTGeophysics

run = VFSA2D("model.start", "data.dat", "FwdCtrl", "InvCtrl.VFSA", "mask.ctrl")
PlotInversion2D(run; true_model_path = "model.true")
```

## What you get

```text
run_YYYYmmdd_HHMMSS/
├── model.rho        ensemble mean, the result of the run
├── data.pred        its predicted data
├── Summary.txt      rms of the mean and of each chain
├── inputs/          the five input files
├── plots/           mean, median, best, 5 % and 95 % models, standard deviation, data fit, convergence
└── vfsa/            each chain's best model and history, and the ensemble statistics
```

The standard deviation and the 5–95 % models are the ones to look at for uncertainty. Where they
are wide, the data do not constrain the model well.

To recompute the ensemble statistics from the chains later, for example after removing a chain:

```bash
julia --project=. helpers/run_statistics_2D.jl run_dir
```

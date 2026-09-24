# 2D VFSA Inversion

Very fast simulated annealing over Gaussian RBF control points (`src/Inv2D_VFSA.jl`).
It runs from five files: the start model, the observed data, `fwd.ctrl`, the VFSA
control and `mask.ctrl`. There is no model covariance and no prior. VFSA has no
regularization term, so the start model is the centre of the search, and the mask
alone says which cells move. The mesh, topography, water, data rows, misfit and run
folder are those of Gauss–Newton and NLCG ([deterministic inversion](gaussnewton2d.md)).
Independent chains run on the Julia threads. Their best models form the ensemble, and
its mean, median, spread and 5–95 % range estimate the uncertainty.

```bash
julia --project=. -t 10 examples/run_vfsa2D.jl                  # examples/data/2D-IV
julia --project=. -t 10 examples/run_vfsa2D.jl model.start data.dat FwdCtrl InvCtrl.VFSA mask.ctrl
```

```julia
run = VFSA2D("model.start", "data.dat", "FwdCtrl", "InvCtrl.VFSA", "mask.ctrl")
PlotInversion2D(run; true_model_path = "model.true")
```

Start Julia with as many threads as chains (`-t`). The sparse solves then run single
threaded.

## Mask

`mask.ctrl` is our own format. Its first line is `ny nz`, followed by `nz` rows of `ny`
integers, top row first. It uses the same codes as the covariance mask:

| Value | Cells |
|:------|:------|
| 0 | air (topographic air must be 0) or fixed |
| 9 | water: fixed at the model's value; no station may stand over it |
| other | free |

`MakeMesh2D` and the benchmark helpers write `mask.ctrl` next to `cov.ctrl`, holding
the same mask. `ReadMask2D` and `WriteMask2D` read and write it, and `Mask2D` builds one
from a model.

## Controls

`InvCtrl.VFSA` (`examples/ctrl/2D/InvCtrl.VFSA`) holds VFSA keys only; there is no
`Algorithm` and no lambda:

```text
Exit search when rms is less than : 1
Maximum number of iterations      : 3000
Mode                              : TETM
Log10 resistivity bounds          : 0 4
Number of chains                  : 2
Control points                    : 400
Trials per iteration              : 1
Share of controls moved           : 0.2
Step scale                        : 0.2
Starting temperature              : 0.01
Cooling ratio                     : 0.001
RBF width top (cells)             : 2
RBF width bottom (cells)          : 3
Control depth power               : 0.2
Core depth (skin depths)          : 2
Core depth (layers)               : 0
Core expansion (cells)            : 4
Padding decay (core cells)        : 8
Random seed                       : 20260308
Snapshot interval                 : 0
```

| Key | Meaning |
|:----|:--------|
| `Maximum number of iterations` | iterations per chain |
| `Exit search when rms is less than` | a chain stops once its best rms reaches it |
| `Log10 resistivity bounds` | the search box |
| `Number of chains` | independent chains (and ensemble size) |
| `Control points` | RBF controls per chain, drawn among the free core cells |
| `Trials per iteration` | proposals per iteration, the best takes one Metropolis test |
| `Share of controls moved` | controls perturbed per proposal (3D `frac_update_controls`) |
| `Step scale` | proposal width as a share of the box |
| `Starting temperature`, `Cooling ratio` | start temperature, on the scale of a typical relative uphill change of rms², and its ratio at the last iteration |
| `RBF width top (cells)`, `RBF width bottom (cells)` | kernel widths at the top and bottom of the core, linear in depth between (3D `sigma_scale`, `sigma_scale_deep`) |
| `Control depth power` | control placement weight (depth + z₁)^(−p), 0 = uniform (3D `ctrl_depth_power`) |
| `Core depth (skin depths)` | core depth in skin depths of the data: median off-diagonal ρa, longest period |
| `Core depth (layers)` | core depth as the top N layers instead, when N > 0 |
| `Core expansion (cells)` | cells added to each side of the uniform lateral core (3D `core_expand_cells`) |
| `Padding decay (core cells)` | e-fold of the blend from the core edge back to the start model |
| `Random seed` | chain k uses seed + 1000(k−1) |
| `Snapshot interval` | write each chain's best model every N iterations (0 = off) |

`ReadVFSACtrl2D` and `WriteVFSACtrl2D` read and write the file. `VFSA2DConfig(ctrl)`
turns it into the in-memory configuration.

The parameterisation is the one of 3D VFSA (`VFSA3DMT`). Controls sit in the core only:
the uniform lateral block (plus `Core expansion`), down to the core depth. Outside it, as
in 3D, the lateral padding is blended row by row from the median of the edge core columns
back to the start model, and the cells below the core carry its bottom value down,
keeping a third per layer. Air, water and mask-0 cells never change. The energy is the
chi2 of Gauss–Newton, and the temperature schedule is the 3D one.

## Output

```text
run_YYYYmmdd_HHMMSS/
├── model.rho            ensemble mean (log10), the run's model
├── data.pred            its prediction
├── Summary.txt          rms, best chain, per-chain acceptance and rms
├── inputs/              the five input files
├── plots/               ModelFinal (mean), ModelMedian, ModelBest, ModelP05, ModelP95,
│                        ModelStd, DataFit, DataFitBest, ConvergenceVFSA, Mesh, ...
└── vfsa/
    ├── chain_XX/best.rho, History.csv, best_iter_NNNNN.rho
    ├── model.mean.rho, model.median.rho, model.p05.rho, model.p95.rho, model.best.rho
    ├── Uncertainty.csv  per cell: log10 mean, median, std, p05, p95
    ├── Chains.csv
    └── data.best.pred
```

```bash
julia --project=. helpers/run_statistics_2D.jl run_dir     # recompute the ensemble from chain_XX/best.rho
julia --project=. helpers/make_gif_2D.jl run_dir           # convergence GIF from the snapshots
```

## From Julia

```julia
result = VFSA2D(mesh, ρ0, observed; config = VFSA2DConfig(n_chains = 10, n_ctrl = 400, max_iter = 4000),
                active_cells = active, water_cells = water, run_dir = "run")
result.resistivity          # ensemble mean
result.ensemble.std         # log10 standard deviation per cell
ens = mt2d_ensemble(models) # statistics of any set of models
```

The in-memory method writes only the chains and the ensemble (`run_dir/vfsa/`). The
five-file method adds `model.rho`, `data.pred`, `Summary.txt` and `inputs/`.

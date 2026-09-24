# 2D Deterministic Inversion

Gauss–Newton and NLCG share one driver and one six-file front end, as ModEM does. VFSA
has its own five-file front end, with a mask but no covariance or prior
([2D VFSA Inversion](inversion2d.md)). 1D reuses both ([1D](forward1d.md)).

| File | Contents |
|:-----|:---------|
| `src/Control2D.jl` | `fwd.ctrl`, `inv.ctrl`, VFSA control, `cov.ctrl` and `mask.ctrl` readers and writers |
| `src/Mesh2D.jl` | mesh geometry, topography helpers, model builders and model I/O |
| `src/Fwd2D.jl` | TE/TM solver, data I/O, misfit, Fréchet derivatives (G, G δm, Gᵗ δd̂) |
| `src/Topo2D.jl` | `topo.dat`, topography and water cut into models |
| `src/Inv2D.jl` | problem setup, objective, regularization, line search, driver, six-file `Invert2D` |
| `src/Inv2D_GN.jl` | damped Gauss–Newton |
| `src/Inv2D_NLCG.jl` | preconditioned nonlinear conjugate gradients |
| `src/Inv2D_VFSA.jl` | very fast simulated annealing, its ensemble and the five-file `VFSA2D` |
| `src/PlotModel2D.jl` | resistivity sections and the mesh |
| `src/PlotData2D.jl` | data curves and fit, convergence, `PlotInversion2D` |

## Quick start

```bash
julia --project=. helpers/benchmarks_2D.jl                 # examples/data/2D-IV
julia --project=. examples/run_inv2D.jl GN                  # or NLCG
julia --project=. examples/run_inv2D.jl model.start data.dat FwdCtrl InvCtrl cov.ctrl model.prior
```

```julia
run = Invert2D("model.start", "data.dat", "FwdCtrl", "InvCtrl.GN", "cov.ctrl", "model.prior")
PlotInversion2D(run; true_model_path = "model.true")
```

The inputs are the start model, the observed data, `fwd.ctrl` (mode and air), `inv.ctrl`
(`Algorithm : GN` or `NLCG`, and settings), the covariance file and the prior (the
reference model of the regularization). The regularizer is the gradient one below, with
its weights in `inv.ctrl`. The covariance file keeps the ModEM layout, but only its mask
is read (0 = air or fixed, 9 = water, others free).
The controls ship in `examples/ctrl/2D`; `MakeMesh2D` ([mesh tool](mesh2d.md)) writes
all of them for a new data file.

Everything goes to `run_YYYYmmdd_HHMMSS/` next to the data: `model.rho` (restartable),
`data.pred`, `History.csv`, `Summary.txt`, the inputs in `inputs/`, and the plots of
`PlotInversion2D` in `plots/` (mesh, start, final and true models, data fit, convergence).

`inv.ctrl` (`examples/ctrl/2D/InvCtrl.GN`):

```text
Algorithm                         : GN
Initial damping factor lambda     : 1
Exit search when rms is less than : 1
Maximum number of iterations      : 20
Mode                              : TETM
Log10 resistivity bounds          : 0 4
Max log10 step                    : 0.5
Max line search steps             : 12
Smallness weight                  : 0.01
Smoothing weight y                : 1
Smoothing weight z                : 1
GN damping                        : 0.01
```

lambda is β, fixed through the run. NLCG replaces `GN damping` with `NLCG restart` and
`NLCG precondition`.

The in-memory form is `Invert2D(mesh, initial, observed; algorithm, options,
active_cells, reference_resistivity, water_cells)`; `GaussNewton2D(...; config)` and
`NLCG2D(...; config)` are shorthands.

## Objective and options

For `m = log10.(resistivity)` on the active cells:

```math
\Phi(m) = \tfrac12\|W_d(F(m)-d)\|_2^2 + \tfrac{\beta}{2}\|R(m-m_{ref})\|_2^2.
```

Each valid complex datum contributes a real and an imaginary residual, both
divided by the supplied impedance error. RMS is `sqrt(chi2 / n_real_data)`,
which matches `chi2_rms2d`. `R` combines area-weighted smallness with first
differences on the nonuniform mesh. These terms act on the departure from
the reference model, which defaults to the start model.

`Inv2DOptions` holds everything that does not depend on the algorithm:

- `mode`
- `max_iter`
- `beta`, `smallness`, `smooth_y`, `smooth_z`
- `log_bounds`
- `max_step`
- `max_linesearch`
- `target_rms`
- the gradient, step, and objective tolerances
- `verbose`

The driver caps each step at `max_step` (in log₁₀ units), projects it onto
`log_bounds`, and accepts it by Armijo backtracking on the total objective.
Termination reasons are:

- `:target_rms`
- `:gradient_tolerance`
- `:step_tolerance`
- `:objective_tolerance`
- `:max_iter`
- `:line_search_failed`

Without an active mask, all earth cells, including the padding, are inverted.
Air, topographic air included, is fixed. Active cells can be a model-shaped Boolean mask
or a vector of Cartesian `(z, y)` or linear indices. Air and water cells take no
regularization term, and no smoothing pair crosses into them.

## Gauss–Newton

`GaussNewton2DConfig(damping, max_damping_trials)` solves the following
system in each iteration:

```math
(G^tG + \beta R^tR + \lambda S^2)\,\delta m = -g,
```

where `S` is the diagonal column scaling `sqrt(diag(GᵗG + βRᵗR))`. This gives
the same step as the damped augmented least-squares problem.

- λ is divided by 3 after an accepted step.
- λ is multiplied by 10 after a failed line search.
- After `max_damping_trials` failures, the run stops with `:line_search_failed`.

`inv2d_frechet` builds G from one transpose (adjoint) solve per real datum
when there are fewer data than active cells. Otherwise it uses one implicit
tangent linear solve per active cell. The step is a dense `n_active²` Cholesky factorization,
so this algorithm suits small to medium models.

## NLCG

`NLCG2DConfig(restart, precondition, epsilon)` uses preconditioned
Polak–Ribière+ conjugate gradients:

- **Gradient:** one adjoint solve per frequency and mode, via
  `inv2d_gradient(problem, state)`. G is never built, so memory grows
  linearly with the number of cells.
- **Preconditioner:** `P = (βRᵀR + εI)⁻¹`, a sparse Cholesky factorization
  computed once. It smooths updates the same way the regularization does.
- **Restarts:** the direction resets to steepest descent every `restart`
  iterations, and whenever the conjugate direction is not a descent direction.
- **Step length:** the first step is `max_step`. Later steps start from the
  previous one, scaled by `2(g_prev·d_prev)/(g·d)`. The factor of 2 lets a step
  that was cut short grow back. The driver then caps each step and backtracks.
- **Failed search:** a failed conjugate search is retried once as steepest
  descent. If that also fails, the run stops.

NLCG needs more iterations than Gauss–Newton, but each one is much cheaper.

## Adding an algorithm

Create `src/Inv2D_<Name>.jl` and include it after `Inv2D.jl` in
`src/MTGeophysics.jl`. `Inv2D_NLCG.jl` is the smallest complete example. Define
a config type `<: AbstractInversion2D` and these methods:

- `inv2d_tag`
- `inv2d_init`
- `inv2d_prepare!`, which returns the gradient
- `inv2d_direction`
- `inv2d_accept!`

Optional methods have defaults:

- `inv2d_reject!` returns false, so a failed search ends the run.
- `inv2d_info` returns no extra history columns.
- `inv2d_validate` performs no checks.

The driver handles the stopping tests, bounds, step capping, line search,
history, and file output.

Shared building blocks in `Inv2D.jl`:

- `inv2d_frechet(problem, state)`: the weighted `2n_data × n_active` Fréchet derivative
- `inv2d_gradient(problem, state, G)`: `Gᵗr + βRᵗ reg`
- `inv2d_gradient(problem, state)`: the same gradient from adjoint solves, for
  gradient-only methods

## Fréchet derivatives

For the forward relation d = g(m), `Fwd2D.jl` provides the Fréchet derivative
G = ∂g/∂m and its transpose in Tarantola (2005) notation:

```julia
G = FrechetDerivative2D(mesh, ρ;
    active_cells=active, parameterization=:log10_resistivity)
# G.z_xy and G.z_yx: complex (nf*nr) × n_active matrices
# rows follow vec(response.z_xy), frequency fastest
# columns follow G.cells, Cartesian indices in model (z, y) order

δd = ApplyFrechet2D(mesh, ρ, δm; parameterization=:log10_resistivity)             # δd = G δm
δm̂ = ApplyFrechetTranspose2D(mesh, ρ, (z_xy=δẑ,); parameterization=:log10_resistivity) # δm̂ = Gᵗ δd̂
```

G is available for:

- complex impedances
- apparent resistivity
- unfolded phase, in degrees

They are derivatives of the discrete forward equations, computed with reused
LU factors. No model perturbations or finite-difference reruns are involved.
The small boundary-field and receiver-sampling kernels use ForwardDiff dual
numbers. For complex output weights, the adjoint pairing is
`real(dot(weight, delta))`.

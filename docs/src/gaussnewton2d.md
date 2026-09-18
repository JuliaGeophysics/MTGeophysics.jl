# 2D deterministic inversion

The 2D inversion code is split into a framework and separate algorithm files.

| File | Contents |
|:-----|:---------|
| `src/Mesh2D.jl` | Mesh geometry, skin-depth vertical core, model builders, model I/O |
| `src/Fwd2D.jl` | TE/TM forward solver, data I/O, misfit, Fréchet derivatives (G, G δm, Gᵗ δd̂) |
| `src/Inv2D.jl` | Problem setup, objective, regularization, line search, driver `Invert2D` |
| `src/Inv2D_GN.jl` | Damped Gauss–Newton algorithm |
| `src/Inv2D_NLCG.jl` | Preconditioned nonlinear conjugate gradients |
| `src/PlotModel2D.jl` | Model sections, mesh layout, convergence |
| `src/PlotData2D.jl` | Data pseudo-sections, site curves, observed-vs-predicted fit |

`Inv2D.jl` does not depend on any particular algorithm. Future algorithms, such as
L-BFGS or Occam, go in their own `Inv2D_<Name>.jl`.

## Quick start

All inversions take the same inputs: an observed data file and a start model,
whose grid defines the mesh. Generate the COMEMI-III benchmark once. It uses a
skin-depth mesh, 17 frequencies from 0.1 to 1000 Hz, 11 sites, and 5% noise,
all set in `BENCHMARK_MESH` in the helper:

```bash
julia --project=. helpers/benchmarks_2D.jl        # writes examples/0COMEMI2D-III/
```

Then invert it with any algorithm:

```bash
julia --project=. examples/run_inv2D.jl gn        # Gauss–Newton
julia --project=. examples/run_inv2D.jl nlcg      # nonlinear conjugate gradients
julia --project=. examples/run_inv2D.jl vfsa      # very fast simulated annealing
```

Each run writes its results to `examples/Results/<alg>2D_<timestamp>/`. The
`plots/` subfolder gets the same plots for every algorithm:

- `mesh_full`, `mesh_core`
- `mstart_core`, `mtrue_core`
- `mfinal_core`, `mfinal_full`
- `data_fit`
- `data_obs_maps`, `data_pred_maps`
- `convergence` (Gauss–Newton and NLCG only)

In Julia:

```julia
using MTGeophysics

result = Invert2D("Comemi2D3.ini", "Comemi2D3.obs"; output_dir = "Results/run",
    algorithm = NLCG2DConfig(),                 # or GaussNewton2DConfig()
    options = Inv2DOptions(max_iter = 200, beta = 1.0, target_rms = 1.0))
```

The in-memory form is `Invert2D(mesh, initial, observed; algorithm, options)`.
`GaussNewton2D(...; config)` and `NLCG2D(...; config)` are shorthands.

With `output_dir`, the file workflow writes `model_<tag>.rho`, `data_<tag>.dat`,
`history_<tag>.csv`, and `summary_<tag>.txt`. The tag is `gn` or `nlcg`,
depending on the algorithm.
Existing result files are never overwritten. Input data must contain observed
impedances and absolute errors, not a `.ref` template of relative errors.

## Mesh: skin-depth vertical core

`BuildMesh2D` designs the ground layers from the survey unless you pass
explicit `ground_layers`:

- the core is **uniform** (constant dz) from the surface to at least
  `z_core_skin_depths` (default 1) skin depths of the lowest frequency;
- the skin depth is `δ = sqrt(2ρ/(ωμ₀)) ≈ 503·sqrt(ρ/f)`, using `background_resistivity`;
- dz defaults to `max(δ_min/3, δ_max/max_core_layers)`, or you can set `z_core_cell`;
- below the core, layers grow by `pad_factor` until the mesh bottom reaches
  `z_bottom_skin_depths` (default 4) skin depths.

`plot_mt2d_mesh(mesh; region=:full)` shows the whole mesh, with the uniform core
outlined and the skin depths of `f_min` and `f_max` marked. `region=:core`
zooms in on the core.

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
Air is fixed. Active cells can be a model-shaped Boolean mask or a vector of
Cartesian `(z, y)` or linear indices.

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

# 2D inversion theory

The objective, the two deterministic algorithms, the VFSA parameterisation and the Fréchet
derivatives behind the 2D inversions.

## Objective

For `m = log10.(resistivity)` on the active cells:

```math
\Phi(m) = \tfrac12\|W_d(F(m)-d)\|_2^2 + \tfrac{\beta}{2}\|R(m-m_{ref})\|_2^2.
```

Each valid complex datum contributes a real and an imaginary residual, both divided by the supplied
impedance error. RMS is `sqrt(chi2 / n_real_data)`, which matches `chi2_rms2d`. `R` combines
area-weighted smallness with first differences on the nonuniform mesh. These terms act on the
departure from the reference model, which defaults to the start model.

Without an active mask, all earth cells, including the padding, are inverted. Air, topographic air
included, is fixed. Active cells can be a model-shaped Boolean mask or a vector of Cartesian `(z, y)`
or linear indices. Air and water cells take no regularisation term, and no smoothing pair crosses
into them.

## Options and stopping

`Inv2DOptions` holds everything that does not depend on the algorithm: `mode`, `max_iter`, `beta`,
`smallness`, `smooth_y`, `smooth_z`, `max_step`, `max_linesearch`, `target_rms`, the gradient, step
and objective tolerances, and `verbose`.

The driver caps each step at `max_step` (in log₁₀ units) and accepts it by Armijo backtracking on
the total objective. There are no bounds; the regularisation and the step cap keep the model in
range. A run stops with one of `:target_rms`, `:gradient_tolerance`, `:step_tolerance`,
`:objective_tolerance`, `:max_iter` or `:line_search_failed`.

The in-memory entry point is `Invert2D(mesh, initial, observed; algorithm, options, active_cells,
reference_resistivity, water_cells)`; `GaussNewton2D(...; config)` and `NLCG2D(...; config)` are
shorthands.

## Gauss–Newton

`GaussNewton2DConfig(damping, max_damping_trials)` solves in each iteration

```math
(G^tG + \beta R^tR + \lambda S^2)\,\delta m = -g,
```

where `S` is the diagonal column scaling `sqrt(diag(GᵗG + βRᵗR))`. This gives the same step as the
damped augmented least-squares problem.

- λ is divided by 3 after an accepted step.
- λ is multiplied by 10 after a failed line search.
- After `max_damping_trials` failures, the run stops with `:line_search_failed`.

`inv2d_frechet` builds G from one transpose (adjoint) solve per real datum when there are fewer data
than active cells, and otherwise from one tangent linear solve per active cell. The step is a dense
`n_active²` Cholesky factorisation, so the algorithm suits small to medium models.

## NLCG

`NLCG2DConfig(restart, precondition, epsilon)` uses preconditioned Polak–Ribière+ conjugate
gradients:

- **Gradient:** one adjoint solve per frequency and mode, via `inv2d_gradient(problem, state)`.
  G is never built, so memory grows linearly with the number of cells.
- **Preconditioner:** `P = (βRᵀR + εI)⁻¹`, a sparse Cholesky factorisation computed once. It smooths
  updates the same way the regularisation does.
- **Restarts:** the direction resets to steepest descent every `restart` iterations, and whenever
  the conjugate direction is not a descent direction.
- **Step length:** the first step is `max_step`. Later steps start from the previous one, scaled by
  `2(g_prev·d_prev)/(g·d)`; the factor of 2 lets a step that was cut short grow back. The driver then
  caps each step and backtracks.
- **Failed search:** a failed conjugate search is retried once as steepest descent. If that also
  fails, the run stops.

## VFSA parameterisation

The 2D parameterisation is the one of 3D VFSA (`VFSA3DMT`). Controls sit in the core only: the
uniform lateral block (plus `Core expansion`), down to the core depth. Outside it the lateral padding
is blended row by row from the median of the edge core columns back to the start model, and the
cells below the core carry its bottom value down, keeping a third per layer. Air, water and mask-0
cells never change. The energy is the χ² of Gauss–Newton, and the temperature schedule is the 3D one.

In memory:

```julia
result = VFSA2D(mesh, ρ0, observed; config = VFSA2DConfig(n_chains = 10, n_ctrl = 400, max_iter = 4000),
                active_cells = active, water_cells = water, run_dir = "run")
result.resistivity          # ensemble mean
result.ensemble.std         # log10 standard deviation per cell
ens = mt2d_ensemble(models) # statistics of any set of models
```

The in-memory method writes only the chains and the ensemble (`run_dir/vfsa/`). The five-file method
adds `model.rho`, `data.pred`, `Summary.txt` and `inputs/`.

## Fréchet derivatives

For the forward relation d = g(m), `Fwd2D.jl` provides the Fréchet derivative G = ∂g/∂m and its
transpose, in Tarantola (2005) notation:

```julia
G = FrechetDerivative2D(mesh, ρ; active_cells = active, parameterization = :log10_resistivity)
# G.z_xy and G.z_yx: complex (nf*nr) × n_active matrices
# rows follow vec(response.z_xy), frequency fastest
# columns follow G.cells, Cartesian indices in model (z, y) order

δd = ApplyFrechet2D(mesh, ρ, δm; parameterization = :log10_resistivity)                  # δd = G δm
δm̂ = ApplyFrechetTranspose2D(mesh, ρ, (z_xy = δẑ,); parameterization = :log10_resistivity) # δm̂ = Gᵗ δd̂
```

G is available for complex impedances, apparent resistivity and unfolded phase (in degrees). These
are derivatives of the discrete forward equations, computed with reused LU factors, without model
perturbations or finite-difference reruns. The small boundary-field and receiver-sampling kernels
use ForwardDiff dual numbers. For complex output weights, the adjoint pairing is
`real(dot(weight, delta))`.

# Adding an algorithm

The 2D deterministic inversion is built so that a new algorithm only has to say how it chooses a
direction. The driver in `Inv2D.jl` does the rest: stopping tests, step capping, line search,
history and file output.

## Steps

1. Create `src/Inv2D_<Name>.jl` and include it after `Inv2D.jl` in `src/MTGeophysics.jl`.
   `Inv2D_NLCG.jl` is the smallest complete example.
2. Define a config type `<: AbstractInversion2D`.
3. Define these methods for it:
   - `inv2d_tag`
   - `inv2d_init`
   - `inv2d_prepare!`, which returns the gradient
   - `inv2d_direction`
   - `inv2d_accept!`

Optional methods have defaults:

- `inv2d_reject!` returns false, so a failed search ends the run.
- `inv2d_info` returns no extra history columns.
- `inv2d_validate` performs no checks.

## Building blocks

`Inv2D.jl` provides:

- `inv2d_frechet(problem, state)`: the weighted `2n_data × n_active` Fréchet derivative
- `inv2d_gradient(problem, state, G)`: `Gᵗr + βRᵗ reg`
- `inv2d_gradient(problem, state)`: the same gradient from adjoint solves, for gradient-only methods

See [2D inversion theory](inversion2d.md) for the objective these act on.

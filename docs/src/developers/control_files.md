# Control file reference

Every key of every control file, with the settings of the mesh tools and of the 3D VFSA
configuration. The user pages show only the keys that are changed most often.

## 1D inversion (`InvCtrl.GN`, `InvCtrl.VFSA`)

| Key | Meaning |
|:----|:--------|
| `Algorithm` | `GN` or `VFSA` |
| `Mode` | fitted impedance: `XY` (Zxy), `YX` (Zyx = −Z), `XYYX` (both) or `DET` (√det Z) |
| `Exit search when rms is less than` | target rms |
| `Maximum number of iterations` | iterations (per chain for VFSA) |
| `Initial damping factor lambda` | GN only: regularisation weight |
| `Log10 resistivity bounds` | VFSA only: the search box; GN is unbounded |
| `Number of chains` | VFSA only: independent chains, and ensemble size |

Everything else is fixed in the code:

- GN: vertical smoothing with smallness 0.01, Levenberg–Marquardt damping 0.01, and a step of at
  most 0.5 in log10 ρ.
- VFSA: every layer is a parameter (no sparse parameterisation), a fifth of them moved per proposal;
  temperature 0.03, cooling ratio 0.001, step scale 0.11. `model.rho` is the ensemble mean, and
  `vfsa/data.best.pred` holds the best chain's response.

The 1D mesh settings (first layer δ(f_max)/5, growth 1.1, depth 4 δ(f_min)) are constants at the top
of `examples/run_inv1D.jl`.

## 2D forward (`fwd.ctrl`)

| Key | Default | Meaning |
|:----|:--------|:--------|
| `Mode` | `TETM` | `TE`, `TM` or `TETM` |
| `Strike (deg)` | `auto` | `auto` or degrees clockwise from north (see [2D solver notes](solver2d.md)) |
| `Air layers` | 10 | number of air layers |
| `Air thickness (m)` | 50000 | total air thickness |
| `Air growth factor` | 2 | growth of the air layers upwards |
| `Air resistivity (ohm m)` | 1e9 | resistivity given to all air cells, topographic air included |
| `Write Frechet derivative` | `no` | `yes` also writes G = ∂d/∂m as `data.frechet` |
| `Dipole length (m)` | 100 | TM electric dipole length, only used next to topographic steps |

`fwd.ctrl` is required: it defines the air.

## 2D Gauss–Newton and NLCG (`inv.ctrl`)

| Key | Meaning |
|:----|:--------|
| `Algorithm` | `GN` or `NLCG` |
| `Initial damping factor lambda` | β, the regularisation weight, fixed through the run |
| `Exit search when rms is less than` | target rms |
| `Maximum number of iterations` | iteration cap |
| `Mode` | `TE`, `TM` or `TETM` |
| `Max log10 step` | largest change of any cell in one step, in log10 ρ |
| `Max line search steps` | backtracking steps before a step is rejected |
| `Smallness weight` | weight of the smallness term in the regulariser |
| `Smoothing weight y`, `Smoothing weight z` | weights of the horizontal and vertical first differences |
| `GN damping` | GN only: initial Levenberg–Marquardt λ |
| `NLCG restart` | NLCG only: restart to steepest descent every N iterations |
| `NLCG precondition` | NLCG only: `yes` uses the regulariser as preconditioner |

GN and NLCG are unbounded: `Log10 resistivity bounds` belongs to VFSA and is an error in `inv.ctrl`.

## 2D VFSA (`vfsa.ctrl`)

| Key | Meaning |
|:----|:--------|
| `Maximum number of iterations` | iterations per chain |
| `Exit search when rms is less than` | a chain stops once its best rms reaches it |
| `Mode` | `TE`, `TM` or `TETM` |
| `Log10 resistivity bounds` | the search box |
| `Number of chains` | independent chains (and ensemble size) |
| `Control points` | RBF controls per chain, drawn among the free core cells |
| `Trials per iteration` | proposals per iteration; the best takes one Metropolis test |
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

`ReadVFSACtrl2D` and `WriteVFSACtrl2D` read and write the file, and `VFSA2DConfig(ctrl)` turns it into
the in-memory configuration.

## `mask.ctrl`

Our own format: a first line `ny nz`, then `nz` rows of `ny` integers, top row first. The values are
those of the covariance mask: 0 = air or fixed, 9 = water, others free. `cov.ctrl` keeps the ModEM
covariance layout, but only its mask is read.

## `MakeMesh2D`

| Setting | Default | Meaning |
|:--------|:--------|:--------|
| `cell_width_frac` | 0.5 | core cell width as a fraction of the median station spacing |
| `core_margin_cells` | 4 | core cells beyond the outer stations |
| `n_pad`, `pad_factor` | 12, 1.3 | lateral padding cells each side and their growth |
| `first_layer_div` | 5 | first layer = δ(f_max) / this |
| `vertical_factor` | 1.1 | layer growth with depth |
| `depth_mult` | 4 | model depth = this × δ(f_min) |
| `background_resistivity` | 0 | 0 = median apparent resistivity of the data |
| `air_layers`, `air_thickness`, `air_growth`, `air_resistivity` | 10, 50 km, 2, 1e9 | the air in `fwd.ctrl` |
| `dipole_length` | 100 | TM dipole averaging next to topographic steps |
| `cov_smoothing`, `n_smooth` | 0.3, 1 | smoothing values kept in the ModEM covariance layout |
| `fixed_below_m` | Inf | fix cells below this depth (mask 0) |
| `water_resistivity` | 100 | lakes and sea |
| `strike` | `nothing` | `nothing` = auto from the data, or degrees clockwise from north; written to `fwd.ctrl` |
| `inv_ctrl` | `examples/ctrl/2D/InvCtrl.GN` | copied as `inv.ctrl` |
| `vfsa_ctrl` | `examples/ctrl/2D/InvCtrl.VFSA` | copied as `vfsa.ctrl` |

The core is centred on y = 0 (the data's origin), as ModEM expects.

## `VFSA3DMTConfig`

Fields without a default must be given; `examples/run_vfsa3D.jl` sets all of them.

| Field | Default | Meaning |
|:------|:--------|:--------|
| `nprocs` | | MPI processes for each ModEM forward call |
| `mpirun_cmd` | | MPI launcher, e.g. `"mpirun"` or `"srun"` |
| `modem_exe` | | ModEM executable |
| `fwd_ctrl` | `""` | ModEM forward control file; `""` uses the binary's compiled-in solver defaults |
| `out_root` | | run folder base name; `_<timestamp>` is appended, relative paths sit next to the start model |
| `n_ctrl` | | Gaussian RBF control points in the core |
| `frac_update_controls` | | share of controls perturbed per trial |
| `log_bounds` | | search box in log10 Ω·m |
| `step_scale` | | largest control jump, as a share of the box width |
| `max_iter` | | iteration cap; also sets the cooling timescale |
| `n_trials` | | trials per iteration; the best takes one Metropolis test (1 = classic VFSA) |
| `T0`, `cool_ratio` | | start temperature, and the final temperature as a share of it |
| `target_rms` | | early-stop rms |
| `seed` | | random seed |
| `pad_tol` | | tolerance of the core/padding detection |
| `core_expand_cells` | 0 | grow the lateral core by N cells per side |
| `padding_decay_length` | | lateral padding blend e-fold, in core cells |
| `z_core_skin_depths` | 1.0 | core depth in data skin depths (`Inf` = full column) |
| `z_core_cells` | 0 | core depth as the top N layers instead, when N > 0 |
| `keep_models`, `keep_dpred` | | keep every trial model and predicted data file |
| `model_save_every` | | > 0 keeps the winning trial model every N iterations |
| `sigma_scale`, `sigma_scale_deep` | | RBF 1σ in cells at the top and bottom of the core |
| `trunc_sigmas` | | RBF kernels are cut beyond this many σ |
| `ctrl_depth_power` | | control placement weight (depth + z₁)^(−p) |
| `water_log10` | | start-model cells below this log10 ρ are frozen water (`NaN` = off) |
| `bathymetry_file` | | frozen water from a bathymetry file instead (`""` = off) |
| `distortion_mode` | `:off` | `:on` fits a per-site 2 × 2 galvanic distortion matrix at every misfit evaluation |
| `distortion_damping` | 0.01 | pull of that matrix towards the identity; `Inf` = correction off |

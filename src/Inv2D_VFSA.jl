# 2D MT VFSA inversion
# Author: @pankajkmishra
# Very fast simulated annealing over Gaussian-RBF control points, parameterised as VFSA3DMT does in 3D: controls
# only in the core (the uniform lateral block, optionally grown by whole cells, down to a skin-depth or layer
# limit), the lateral padding blended from the core edge back to the start model and the cells below the core
# carried down a third per layer. Five files (start, data, fwd.ctrl, the VFSA control and mask.ctrl): no
# covariance and no prior, the mask says which cells are frozen. The mesh, topography, water, data rows and run
# directory are those of GN and NLCG. Independent chains run on threads; their best models form the ensemble
# whose mean, median, spread and 5-95% range are the uncertainty estimate

using Dates
using LinearAlgebra
using Printf
using Random
using SparseArrays
using Statistics

"""
    VFSA2DConfig(; kwargs...)

VFSA controls, named as in `VFSA3DMTConfig`.
- `n_chains`, `max_iter`, `n_trials`: independent chains, iterations per chain, proposals per iteration
- `n_ctrl`: RBF control points per chain, drawn among the free core cells
- `log_bounds`, `step_scale`, `frac_update_controls`: log10 ρ box, proposal width as a share of
  the box, share of the controls moved per proposal
- `T0`, `cool_ratio`: temperature at the start and its ratio at `max_iter`; one schedule
  drives proposal width and acceptance
- `target_rms`: a chain stops once its best rms reaches it
- `pad_tol`, `core_expand_cells`: lateral core, the uniform-width block (as in 3D), grown
  by `core_expand_cells` per side
- `z_core_skin_depths`, `z_core_cells`: core depth, in skin depths of the data (median
  off-diagonal apparent resistivity, longest period; Inf = whole model) or as the top
  `z_core_cells` layers when that is positive
- `sigma_scale`, `sigma_scale_deep`, `trunc_sigmas`: RBF widths in cells at the top and
  the bottom of the core (linear in depth between them) and their truncation
- `ctrl_depth_power`: control placement weight (depth + z₁)^(-p), 0 = uniform
- `padding_decay_length`: e-fold, in core cells, of the blend from the core edge to the start model
- `mode`, `seed`, `snapshot_interval` (0 = off), `verbose`
"""
Base.@kwdef struct VFSA2DConfig
    n_chains::Int = 2
    n_ctrl::Int = 400
    max_iter::Int = 3000
    n_trials::Int = 1
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    frac_update_controls::Float64 = 0.2
    step_scale::Float64 = 0.2
    T0::Float64 = 1.0
    cool_ratio::Float64 = 1e-3
    target_rms::Float64 = 1.0
    seed::Int = 20260308
    pad_tol::Float64 = 0.20
    core_expand_cells::Int = 0
    z_core_skin_depths::Float64 = 1.0
    z_core_cells::Int = 0
    sigma_scale::Float64 = 2.0
    sigma_scale_deep::Float64 = sigma_scale
    trunc_sigmas::Float64 = 3.0
    ctrl_depth_power::Float64 = 0.0
    padding_decay_length::Float64 = 8.0
    mode::Symbol = :TETM
    snapshot_interval::Int = 0
    verbose::Bool = true
end

function _vfsa2d_validate(c::VFSA2DConfig)
    c.n_chains >= 1 && c.n_ctrl >= 1 && c.max_iter >= 0 && c.n_trials >= 1 ||
        throw(ArgumentError("VFSA needs n_chains, n_ctrl, n_trials ≥ 1 and max_iter ≥ 0"))
    c.log_bounds[1] < c.log_bounds[2] || throw(ArgumentError("log_bounds must be increasing"))
    0 < c.frac_update_controls <= 1 || throw(ArgumentError("frac_update_controls must lie in (0, 1]"))
    c.T0 > 0 && 0 < c.cool_ratio <= 1 || throw(ArgumentError("need T0 > 0 and 0 < cool_ratio ≤ 1"))
    c.core_expand_cells >= 0 && c.z_core_cells >= 0 && c.z_core_skin_depths > 0 ||
        throw(ArgumentError("need core_expand_cells ≥ 0, z_core_cells ≥ 0 and z_core_skin_depths > 0"))
    c.sigma_scale > 0 && c.sigma_scale_deep > 0 && c.trunc_sigmas > 0 && c.ctrl_depth_power >= 0 ||
        throw(ArgumentError("need positive RBF widths and truncation, and ctrl_depth_power ≥ 0"))
    c.padding_decay_length > 0 || throw(ArgumentError("padding_decay_length must be positive"))
    c.mode in (:TE, :TM, :TETM) || throw(ArgumentError("mode must be :TE, :TM or :TETM"))
    c
end

"""
    VFSA2DConfig(ctrl::VFSACtrl2D; mode=ctrl.mode, n_ctrl=ctrl.control_points)

VFSA controls from a VFSA control file.
"""
VFSA2DConfig(c::VFSACtrl2D; mode::Symbol = c.mode, n_ctrl::Int = c.control_points) =
    VFSA2DConfig(n_chains = c.chains, n_ctrl = n_ctrl, max_iter = c.max_iter, n_trials = c.trials,
                 log_bounds = c.log_bounds, frac_update_controls = c.share_moved, step_scale = c.step_scale,
                 T0 = c.temperature, cool_ratio = c.cooling, target_rms = c.target_rms, seed = c.seed,
                 core_expand_cells = c.core_expansion, z_core_skin_depths = c.core_skin_depths,
                 z_core_cells = c.core_layers, sigma_scale = c.rbf_top, sigma_scale_deep = c.rbf_bottom,
                 ctrl_depth_power = c.depth_power, padding_decay_length = c.padding_decay, mode = mode,
                 snapshot_interval = c.snapshots)

#---------- core ----------

# lateral core columns and depth core rows (full-mesh indices, earth only), as core_ranges and z_core_range in 3D
function _vfsa2d_core(mesh::MT2DMesh, observed, config::VFSA2DConfig)
    ny, na = length(mesh.y_cell_sizes), mesh.n_air_cells
    c = _core_range(mesh.y_cell_sizes; tol = config.pad_tol)
    e = config.core_expand_cells
    iy = max(1, first(c) - e):min(ny, last(c) + e)
    depth = mt2d_z_centers(mesh)[na+1:end] .- mesh.z_nodes[na+1]
    nz = length(depth)
    last_row = if config.z_core_cells > 0
        min(config.z_core_cells, nz)
    elseif isfinite(config.z_core_skin_depths)
        ρa = filter(x -> isfinite(x) && x > 0, vcat(vec(observed.rho_xy), vec(observed.rho_yx)))
        δ = isempty(ρa) ? Inf : config.z_core_skin_depths * 503.0 * sqrt(median(ρa) * maximum(1 ./ observed.frequencies))
        max(1, searchsortedlast(depth, δ))
    else
        nz
    end
    iy, na+1:na+last_row
end

#---------- rbf parameterisation ----------

# core change as W p, as build_rbf_map in 3D: controls drawn among the free core cells with weight
# (depth + z₁)^(-p) (Efraimidis-Spirakis), gaussian kernels in cell-index space whose width grows linearly
# in depth from sigma_scale to sigma_scale_deep over the core, truncated at trunc_sigmas, normalised per
# cell; a cell out of every kernel's reach takes its nearest control
struct VFSA2DMap
    cells::Vector{CartesianIndex{2}}
    frozen::BitVector
    controls::Vector{CartesianIndex{2}}
    W::SparseMatrixCSC{Float64, Int}
    iy::UnitRange{Int}
    kz::UnitRange{Int}
end

function _vfsa2d_map(mesh::MT2DMesh, protected::AbstractMatrix{Bool}, iy, kz, config::VFSA2DConfig, rng::AbstractRNG)
    cells = vec([CartesianIndex(k, j) for k in kz, j in iy])
    frozen = BitVector([protected[c] for c in cells])
    free = cells[.!frozen]
    isempty(free) && throw(ArgumentError("every core cell is frozen"))
    ztop = mesh.z_nodes[mesh.n_air_cells+1]
    depth = mt2d_z_centers(mesh) .- ztop
    z1 = depth[first(kz)] + eps()
    p = config.ctrl_depth_power
    keys = [rand(rng)^(1 / (p == 0 ? 1.0 : (depth[c[1]] + z1)^(-p))) for c in free]
    controls = free[partialsortperm(keys, 1:min(config.n_ctrl, length(free)); rev = true)]
    config.n_ctrl > length(controls) && @warn "VFSA: only $(length(controls)) of $(config.n_ctrl) controls placed"

    span = depth[last(kz)] - depth[first(kz)]
    σ = [config.sigma_scale + (span > 0 ? clamp((depth[q[1]] - depth[first(kz)]) / span, 0, 1) : 0.0) *
         (config.sigma_scale_deep - config.sigma_scale) for q in controls]
    I, J, V = Int[], Int[], Float64[]
    for (r, c) in enumerate(cells)
        r2 = [((c[1] - q[1])^2 + (c[2] - q[2])^2) / σ[n]^2 for (n, q) in enumerate(controls)]
        near = findall(<=(config.trunc_sigmas^2), r2)
        isempty(near) && (near = [argmin(r2)])
        w = exp.(-0.5 .* r2[near])
        append!(I, fill(r, length(near))); append!(J, near); append!(V, w ./ sum(w))
    end
    VFSA2DMap(cells, frozen, controls, sparse(I, J, V, length(cells), length(controls)), iy, kz)
end

#---------- padding ----------

# below the core each core column continues its last value, a third per layer, towards the start model
function _vfsa2d_decay_z!(m, iy, kz, protected, start)
    k0 = last(kz)
    @inbounds for k in k0+1:size(m, 1), j in iy
        protected[k, j] && continue
        w = (1 / 3)^(k - k0)
        m[k, j] = protected[k0, j] ? start[k, j] : m[k0, j] * w + start[k, j] * (1 - w)
    end
    m
end

# lateral padding: the median of the edge_window core columns at each edge, row by row, blended into the
# start model with e-fold L (m) of the true distance from the edge column
function _vfsa2d_decay_y!(m, iy, yc, L, protected, start; edge_window::Int = 2)
    y1, y2 = first(iy), last(iy)
    buf = Float64[]
    @inbounds for k in axes(m, 1)
        edge = map(((y1, y1:min(y2, y1 + edge_window)), (y2, max(y1, y2 - edge_window):y2))) do (_, cols)
            empty!(buf)
            foreach(j -> protected[k, j] || push!(buf, m[k, j]), cols)
            isempty(buf) ? NaN : median!(buf)
        end
        for j in axes(m, 2)
            (y1 <= j <= y2 || protected[k, j]) && continue
            jc, src = j < y1 ? (y1, edge[1]) : (y2, edge[2])
            w = exp(-abs(yc[j] - yc[jc]) / L)
            m[k, j] = isfinite(src) ? src * w + start[k, j] * (1 - w) : start[k, j]
        end
    end
    m
end

#---------- one chain ----------

# log10 model of control changes p and its data misfit (chi2 over real data, as GN): the core from the RBF
# map within the bounds, frozen cells at the start, then the blends below and beside the core
function _vfsa2d_evaluate(problem, map::VFSA2DMap, p)
    lo, hi = problem.config.log_bounds
    m = copy(problem.m0)
    v = clamp.(problem.m0[map.cells] .+ map.W * p, lo, hi)
    v[map.frozen] .= problem.m0[map.cells[map.frozen]]
    m[map.cells] .= v
    _vfsa2d_decay_z!(m, map.iy, map.kz, problem.protected, problem.m0)
    _vfsa2d_decay_y!(m, map.iy, problem.yc, problem.L, problem.protected, problem.m0)
    ρ = 10.0 .^ m
    response, _ = _mt2d_forward_cache(problem.mesh, ρ; mode = problem.config.mode, cache_fields = false)
    r = _inv2d_residual(response, problem.rows)
    chi2 = sum(abs2, r)
    (; m, ρ, response, chi2, rms = sqrt(chi2 / length(r)))
end

function _vfsa2d_chain(k::Int, problem, dir::Union{Nothing, String})
    config = problem.config
    rng = MersenneTwister(config.seed + 1000 * (k - 1))
    map = _vfsa2d_map(problem.mesh, problem.protected, problem.iy, problem.kz, config, rng)
    lo, hi = config.log_bounds
    base = problem.m0[map.controls]
    p = zeros(length(map.controls))
    n_move = max(1, round(Int, config.frac_update_controls * length(p)))
    current = _vfsa2d_evaluate(problem, map, p)
    best = current
    history = [(chain = k, iteration = 0, temperature = config.T0, trial_rms = current.rms,
                rms = current.rms, best_rms = best.rms, accepted = true)]
    ak = _resolve_ak(config.cool_ratio, config.max_iter)
    every = max(1, config.max_iter ÷ 20)
    chaindir = dir === nothing ? nothing : mkpath(joinpath(dir, @sprintf("chain_%02d", k)))
    for it in 1:config.max_iter
        best.rms <= config.target_rms && break
        T = _T_schedule(it; T0 = config.T0, ak)
        # all trials branch from the current state, the best one takes one metropolis test
        trial, trial_p = nothing, p
        for _ in 1:config.n_trials
            q = copy(p)
            propose_controls!(q, T, lo, hi, base, n_move, rng; step_scale = config.step_scale)
            t = _vfsa2d_evaluate(problem, map, q)
            (trial === nothing || !(t.rms >= trial.rms)) && ((trial, trial_p) = (t, q))
        end
        dE = (trial.rms^2 - current.rms^2) / max(current.rms^2, eps())
        accepted = isfinite(trial.rms) && rand(rng) < (dE <= 0 ? 1.0 : exp(-dE / max(T, 1e-12)))
        if accepted
            p, current = trial_p, trial
            current.chi2 < best.chi2 && (best = current)
        end
        push!(history, (chain = k, iteration = it, temperature = T, trial_rms = trial.rms,
                        rms = current.rms, best_rms = best.rms, accepted))
        if chaindir !== nothing && config.snapshot_interval > 0 && it % config.snapshot_interval == 0
            WriteModel2D(joinpath(chaindir, @sprintf("best_iter_%05d.rho", it)), problem.mesh, best.ρ)
        end
        config.verbose && (it % every == 0 || it == config.max_iter) &&
            @printf("VFSA chain %2d  iter %6d/%d  T %.3g  RMS %.4f  best %.4f\n", k, it, config.max_iter, T, current.rms, best.rms)
    end
    if chaindir !== nothing
        WriteModel2D(joinpath(chaindir, "best.rho"), problem.mesh, best.ρ)
        _vfsa2d_write_history(joinpath(chaindir, "History.csv"), history)
    end
    (; chain = k, best, current, history, n_controls = length(map.controls), controls = map.controls,
       acceptance = length(history) > 1 ? mean(h.accepted for h in history[2:end]) : NaN)
end

function _vfsa2d_write_history(path, history)
    open(path, "w") do io
        println(io, "chain,iteration,temperature,trial_rms,rms,best_rms,accepted")
        for h in history
            @printf(io, "%d,%d,%.6g,%.6f,%.6f,%.6f,%d\n", h.chain, h.iteration, h.temperature, h.trial_rms, h.rms,
                    h.best_rms, h.accepted)
        end
    end
    path
end

#---------- ensemble ----------

"""
    mt2d_ensemble(models) -> (; mean, median, std, p05, p95, count)

Cell-wise statistics of log10 resistivity over an ensemble of equally shaped
resistivity models (ohm m): mean and median, standard deviation and the 5 and 95%
quantiles, all in log10 ρ.
"""
function mt2d_ensemble(models::AbstractVector{<:AbstractMatrix{<:Real}})
    isempty(models) && throw(ArgumentError("empty ensemble"))
    L = cat((log10.(m) for m in models)...; dims = 3)
    q(x, p) = quantile(vec(x), p)
    (mean = dropdims(mean(L; dims = 3); dims = 3), median = dropdims(median(L; dims = 3); dims = 3),
     std = length(models) > 1 ? dropdims(std(L; dims = 3); dims = 3) : zeros(size(L)[1:2]),
     p05 = [q(L[i, j, :], 0.05) for i in axes(L, 1), j in axes(L, 2)],
     p95 = [q(L[i, j, :], 0.95) for i in axes(L, 1), j in axes(L, 2)], count = length(models))
end

# ensemble models, CSV table and chain table into dir
function _vfsa2d_write_ensemble(dir, mesh, ens, active, chains)
    mkpath(dir)
    for (k, name) in ((:mean, "model.mean.rho"), (:median, "model.median.rho"), (:p05, "model.p05.rho"), (:p95, "model.p95.rho"))
        WriteModel2D(joinpath(dir, name), mesh, 10.0 .^ getproperty(ens, k))
    end
    yc, zc = mt2d_y_centers(mesh), mt2d_z_centers(mesh) .- mesh.z_nodes[mesh.n_air_cells+1]
    air = mt2d_air_mask(mesh)
    open(joinpath(dir, "Uncertainty.csv"), "w") do io
        println(io, "iz,iy,y_m,depth_m,active,log10_mean,log10_median,log10_std,log10_p05,log10_p95")
        for iy in axes(active, 2), iz in mesh.n_air_cells+1:size(active, 1)
            air[iz, iy] && continue
            @printf(io, "%d,%d,%.2f,%.2f,%d,%.5f,%.5f,%.5f,%.5f,%.5f\n", iz - mesh.n_air_cells, iy, yc[iy], zc[iz],
                    active[iz, iy], ens.mean[iz, iy], ens.median[iz, iy], ens.std[iz, iy], ens.p05[iz, iy], ens.p95[iz, iy])
        end
    end
    open(joinpath(dir, "Chains.csv"), "w") do io
        println(io, "chain,controls,iterations,acceptance,best_rms,final_rms")
        for c in chains
            @printf(io, "%d,%d,%d,%.4f,%.6f,%.6f\n", c.chain, c.n_controls, c.history[end].iteration, c.acceptance,
                    c.best.rms, c.current.rms)
        end
    end
    dir
end

#---------- driver ----------

"""
    VFSA2D(mesh, initial_resistivity, observed::DataFile2D; config=VFSA2DConfig(),
           active_cells=nothing, water_cells=nothing, run_dir=nothing)

Very fast simulated annealing over RBF control points, parameterised as `VFSA3DMT`: the
controls sit in the core only (the uniform lateral block grown by `core_expand_cells`,
down to `z_core_skin_depths` skin depths or `z_core_cells` layers); the lateral padding is
blended from the core edge back to the start model over `padding_decay_length` core
cells, and the cells below the core carry its bottom down a third per layer. Air, water
and inactive (mask 0) cells never change. The misfit is the error-weighted impedance
chi2 of GN and NLCG. Chains run on the Julia threads, each with its own random controls
and seed. Their best models form the ensemble: mean and median (log10), standard
deviation and 5-95% range.

With `run_dir` the chains write `vfsa/chain_XX/{best.rho, History.csv}` and the ensemble
`vfsa/{model.mean.rho, model.median.rho, model.p05.rho, model.p95.rho, model.best.rho,
Uncertainty.csv, Chains.csv}`. Returns the ensemble mean model (`resistivity`), its
response and fit, the best chain's model, response and fit, the ensemble and the chains.
"""
function VFSA2D(mesh::MT2DMesh, initial_resistivity::AbstractMatrix{<:Real}, observed::DataFile2D;
                config::VFSA2DConfig = VFSA2DConfig(), active_cells = nothing,
                water_cells::Union{Nothing, AbstractMatrix{Bool}} = nothing,
                run_dir::Union{Nothing, AbstractString} = nothing)
    _vfsa2d_validate(config)
    ρ0 = Matrix{Float64}(initial_resistivity)
    air = mt2d_air_mask(mesh)
    ρ0[air] .= mesh.air_resistivity
    cells = _inv2d_cells(mesh, ρ0, active_cells)
    water_cells === nothing || filter!(i -> !water_cells[i], cells)
    active = falses(size(ρ0))
    active[cells] .= true
    rows = _inv2d_data(mesh, observed, config.mode)
    protected = .!active                        # air, water and mask-0 cells: no controls, no change

    # the core carries the controls; the padding and the cells below it follow by the 3D blends
    iy, kz = _vfsa2d_core(mesh, observed, config)
    lo, hi = config.log_bounds
    all(i -> lo <= log10(ρ0[i]) <= hi, (CartesianIndex(k, j) for k in kz, j in iy if !protected[k, j])) ||
        throw(ArgumentError("free starting resistivities of the core are outside log_bounds"))
    L = config.padding_decay_length * median(mesh.y_cell_sizes[_core_range(mesh.y_cell_sizes; tol = config.pad_tol)])
    problem = (; mesh, rows, active, protected, iy, kz, yc = mt2d_y_centers(mesh), L, m0 = log10.(ρ0), config)
    depth = mesh.z_nodes[last(kz)+1] - mesh.z_nodes[mesh.n_air_cells+1]
    config.verbose && @printf("VFSA core: columns %d-%d of %d, layers 1-%d of %d (to %.0f m), %d free core cells\n",
                              first(iy), last(iy), size(ρ0, 2), length(kz), size(ρ0, 1) - mesh.n_air_cells, depth,
                              count(!, protected[kz, iy]))
    dir = run_dir === nothing ? nothing : mkpath(joinpath(String(run_dir), "vfsa"))

    # one chain per task; the sparse solves stay single threaded underneath
    nblas = BLAS.get_num_threads()
    Threads.nthreads() > 1 && BLAS.set_num_threads(1)
    chains = Vector{Any}(undef, config.n_chains)
    try
        Threads.@threads :dynamic for k in 1:config.n_chains
            chains[k] = _vfsa2d_chain(k, problem, dir)
        end
    finally
        BLAS.set_num_threads(nblas)
    end
    chains = identity.(chains)

    ens = mt2d_ensemble([c.best.ρ for c in chains])
    mean_ρ = 10.0 .^ ens.mean
    mean_ρ[air] .= mesh.air_resistivity
    meanfit = let r = _mt2d_forward_cache(mesh, mean_ρ; mode = config.mode, cache_fields = false)[1]
        res = _inv2d_residual(r, rows)
        (; response = r, chi2 = sum(abs2, res), rms = sqrt(sum(abs2, res) / length(res)))
    end
    bestchain = chains[argmin([c.best.chi2 for c in chains])]
    if dir !== nothing
        _vfsa2d_write_ensemble(dir, mesh, ens, active, chains)
        WriteModel2D(joinpath(dir, "model.best.rho"), mesh, bestchain.best.ρ)
    end
    config.verbose && @printf("VFSA ensemble of %d chains: mean-model RMS %.4f, best chain %d RMS %.4f\n",
                              length(chains), meanfit.rms, bestchain.chain, bestchain.best.rms)
    (; resistivity = mean_ρ, response = meanfit.response, rms = meanfit.rms, chi2 = meanfit.chi2,
       best_resistivity = bestchain.best.ρ, best_response = bestchain.best.response, best_rms = bestchain.best.rms,
       best_chain = bestchain.chain, ensemble = ens, chains, active_cells = cells, config)
end

#---------- file workflow ----------

"""
    VFSA2D(start_path, data_path, fwd_path, vfsa_path, mask_path; run_dir=nothing)

VFSA inversion from five files: the start model, the observed data, `fwd.ctrl`, the VFSA
control (`InvCtrl.VFSA`) and `mask.ctrl` (0 = air or fixed, 9 = water, others free).
There is no covariance and no prior: the start model is the centre of the search.
Topographic air must carry mask 0, and no station may stand over water.

Writes into `run_dir` (default `run_YYYYmmdd_HHMMSS/` next to the data) the ensemble
mean as `model.rho` and its prediction as `data.pred`, the chains and ensemble in
`vfsa/` (see the in-memory method) with the best chain's prediction as
`vfsa/data.best.pred`, `Summary.txt` and the inputs. Returns a named tuple like the
six-file `Invert2D`, with the VFSA result in `vfsa`.
"""
function VFSA2D(start_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString,
                vfsa_path::AbstractString, mask_path::AbstractString; run_dir::Union{Nothing, AbstractString} = nothing)
    observed = load_data2d(data_path)
    fwd, ctrl, mask = ReadFwdCtrl2D(fwd_path), ReadVFSACtrl2D(vfsa_path), ReadMask2D(mask_path)
    s = _inv2d_file_setup(ReadModel2D(start_path), observed, fwd, mask, "mask")
    dir = _inv2d_open_run(run_dir, data_path, (start_path, data_path, fwd_path, vfsa_path, mask_path))

    vfsa = VFSA2D(s.mesh, s.ρ0, observed; config = VFSA2DConfig(ctrl), active_cells = s.active,
                  water_cells = any(s.water) ? s.water : nothing, run_dir = dir)
    predicted = _inv2d_predicted(vfsa.response, observed, ctrl.mode)
    fit = chi2_rms2d(observed, predicted; components = _inv2d_components(ctrl.mode))
    converged = fit.rms <= ctrl.target_rms
    reason = converged ? :target_rms : :max_iter

    WriteModel2D(joinpath(dir, "model.rho"), s.mesh, vfsa.resistivity)
    write_data2d(joinpath(dir, "data.pred"), predicted)
    write_data2d(joinpath(dir, "vfsa", "data.best.pred"), _inv2d_predicted(vfsa.best_response, observed, ctrl.mode))
    open(joinpath(dir, "Summary.txt"), "w") do io
        println(io, "Algorithm: VFSA")
        println(io, "Termination: ", reason)
        println(io, "Converged: ", converged)
        @printf(io, "RMS: %.6f\n", fit.rms)
        println(io, "Real data count: ", fit.count)
        println(io, "Model: VFSA ensemble mean (log10) of the chains' best models")
        @printf(io, "Best chain: %d, RMS %.6f\n", vfsa.best_chain, vfsa.best_rms)
        for c in vfsa.chains
            @printf(io, "  chain %2d: %d controls, %d iterations, acceptance %.3f, best RMS %.6f\n", c.chain,
                    c.n_controls, c.history[end].iteration, c.acceptance, c.best.rms)
        end
        _inv2d_summary_cells(io, s)
    end
    (; run_dir = dir, algorithm = :vfsa, ctrl, fwd, mesh = s.mesh, observed, predicted, start = s.ρ0, prior = nothing,
       final = vfsa.resistivity, active = s.active, water = s.water, history = nothing, vfsa, rms = fit.rms, reason,
       converged)
end

"""
    AnalyseEnsemble2D(run_dir) -> (; ensemble, paths)

Recompute the ensemble statistics of a VFSA run from the chains' `best.rho` files in
`run_dir/vfsa`, e.g. after stopping it early; writes the ensemble models next to them.
"""
function AnalyseEnsemble2D(run_dir::AbstractString)
    vdir = isdir(joinpath(run_dir, "vfsa")) ? joinpath(run_dir, "vfsa") : run_dir
    paths = sort(filter(isfile, [joinpath(vdir, d, "best.rho") for d in readdir(vdir) if startswith(d, "chain_")]))
    isempty(paths) && error("no chain_XX/best.rho under $vdir")
    models = ReadModel2D.(paths)
    ens = mt2d_ensemble([m.resistivity for m in models])
    first_model = models[1]
    out = map((:mean, :median, :p05, :p95)) do k
        ρ = 10.0 .^ getproperty(ens, k)
        ρ[first_model.resistivity .> MT2D_AIR_THRESHOLD] .= MT2D_AIR_TAG
        WriteModel2D(joinpath(vdir, "model.$k.rho"), first_model.y_cell_sizes, first_model.z_cell_sizes, ρ)
    end
    (; ensemble = ens, paths = out, chains = paths)
end

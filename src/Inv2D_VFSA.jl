# 2D MT VFSA inversion
# Author: @pankajkmishra
# Very fast simulated annealing over RBF control points, from five files (start, data, fwd.ctrl, the VFSA
# control and mask.ctrl): no covariance and no prior, the mask alone says which cells move. The mesh,
# topography, water, data rows and run directory are those of GN and NLCG. Independent chains run on
# threads; their best models form the ensemble whose mean, median, spread and 5-95% range are the
# uncertainty estimate

using Dates
using LinearAlgebra
using Printf
using Random
using SparseArrays
using Statistics

"""
    VFSA2DConfig(; kwargs...)

VFSA controls.
- `n_chains`, `max_iter`, `n_trials`: independent chains, iterations per chain, proposals per iteration
- `n_ctrl`: RBF control points per chain, drawn at random among the active core cells
- `rbf_sigma_scale_y`, `rbf_sigma_scale_z`: RBF widths in cells; `trunc_sigmas`: truncation in widths
- `log_bounds`, `step_scale`, `frac_update_controls`: log10 ρ box, proposal width as a share of
  the box, share of the controls moved per proposal
- `temp_kappa`, `cool_ratio`: temperature at the start and its ratio at `max_iter`; one
  schedule drives proposal width and acceptance
- `target_rms`: a chain stops once its best rms reaches it
- `padding_decay_length`: decay, in core cells, of the model change into the lateral padding
- `perturb_depth_m`: deepest cell bottom the controls reach (Inf = whole model)
- `pad_tolerance`: width tolerance of the lateral core, as in 3D
- `mode`, `seed`, `snapshot_interval` (0 = off), `verbose`
"""
Base.@kwdef struct VFSA2DConfig
    n_chains::Int = 2
    n_ctrl::Int = 400
    max_iter::Int = 3000
    n_trials::Int = 1
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    frac_update_controls::Float64 = 1.0
    step_scale::Float64 = 0.11
    temp_kappa::Float64 = 1.0
    cool_ratio::Float64 = 1e-3
    target_rms::Float64 = 1.0
    seed::Int = 20260308
    pad_tolerance::Float64 = 0.20
    padding_decay_length::Float64 = 8.0
    rbf_sigma_scale_y::Float64 = 2.0
    rbf_sigma_scale_z::Float64 = 2.5
    trunc_sigmas::Float64 = 3.0
    perturb_depth_m::Float64 = Inf
    mode::Symbol = :TETM
    snapshot_interval::Int = 0
    verbose::Bool = true
end

function _vfsa2d_validate(c::VFSA2DConfig)
    c.n_chains >= 1 && c.n_ctrl >= 1 && c.max_iter >= 0 && c.n_trials >= 1 ||
        throw(ArgumentError("VFSA needs n_chains, n_ctrl, n_trials ≥ 1 and max_iter ≥ 0"))
    c.log_bounds[1] < c.log_bounds[2] || throw(ArgumentError("log_bounds must be increasing"))
    0 < c.frac_update_controls <= 1 || throw(ArgumentError("frac_update_controls must lie in (0, 1]"))
    c.temp_kappa > 0 && 0 < c.cool_ratio <= 1 || throw(ArgumentError("need temp_kappa > 0 and 0 < cool_ratio ≤ 1"))
    c.mode in (:TE, :TM, :TETM) || throw(ArgumentError("mode must be :TE, :TM or :TETM"))
    c
end

"""
    VFSA2DConfig(ctrl::VFSACtrl2D; mode=ctrl.mode, n_ctrl=ctrl.control_points)

VFSA controls from a VFSA control file.
"""
VFSA2DConfig(c::VFSACtrl2D; mode::Symbol = c.mode, n_ctrl::Int = c.control_points) =
    VFSA2DConfig(n_chains = c.chains, n_ctrl = n_ctrl, max_iter = c.max_iter, n_trials = c.trials,
                 log_bounds = c.log_bounds, step_scale = c.step_scale, temp_kappa = c.temperature,
                 cool_ratio = c.cooling, rbf_sigma_scale_y = c.rbf_y, rbf_sigma_scale_z = c.rbf_z, seed = c.seed,
                 target_rms = c.target_rms, mode = mode, snapshot_interval = c.snapshots)

#---------- parameterisation ----------

# model change on the parameterised cells as W p: core cells take normalised gaussian RBF
# weights of the controls within trunc_sigmas (distances in cells, so the widths follow the
# geometric layers), lateral padding cells a decayed copy of the nearest core cell in their row
struct VFSA2DMap
    cells::Vector{CartesianIndex{2}}
    controls::Vector{CartesianIndex{2}}
    W::SparseMatrixCSC{Float64, Int}
end

function _vfsa2d_map(mesh::MT2DMesh, active::AbstractMatrix{Bool}, config::VFSA2DConfig, rng::AbstractRNG)
    core = _core_range(mesh.y_cell_sizes; tol = config.pad_tolerance)
    zbottom = mesh.z_nodes[2:end] .- mesh.z_nodes[mesh.n_air_cells+1]
    deep = isfinite(config.perturb_depth_m) ? zbottom .> config.perturb_depth_m + 1e-9 : falses(length(zbottom))
    usable = copy(active)
    usable[deep, :] .= false
    corecells = [i for i in findall(usable) if i[2] in core]
    isempty(corecells) && throw(ArgumentError("no active cells in the lateral core"))
    controls = corecells[randperm(rng, length(corecells))[1:min(config.n_ctrl, length(corecells))]]

    sy, sz, cut2 = config.rbf_sigma_scale_y, config.rbf_sigma_scale_z, config.trunc_sigmas^2
    weights = Dict{CartesianIndex{2}, Tuple{Vector{Int}, Vector{Float64}}}()
    for c in corecells
        r2 = [((c[2] - k[2]) / sy)^2 + ((c[1] - k[1]) / sz)^2 for k in controls]
        near = findall(<=(cut2), r2)
        isempty(near) && (near = [argmin(r2)])
        w = exp.(-0.5 .* r2[near])
        weights[c] = (near, w ./ sum(w))
    end
    cells = copy(corecells)

    # padding: decayed copy of the nearest core cell of the same row
    yc = mt2d_y_centers(mesh)
    scale = config.padding_decay_length * median(mesh.y_cell_sizes[core])
    for i in findall(usable)
        i[2] in core && continue
        inward = i[2] < first(core) ? (first(core):last(core)) : (last(core):-1:first(core))
        k = findfirst(iy -> haskey(weights, CartesianIndex(i[1], iy)), inward)
        k === nothing && continue
        src = CartesianIndex(i[1], inward[k])
        near, w = weights[src]
        weights[i] = (near, exp(-abs(yc[i[2]] - yc[src[2]]) / max(scale, eps())) .* w)
        push!(cells, i)
    end
    I = reduce(vcat, [fill(r, length(weights[c][1])) for (r, c) in enumerate(cells)])
    J = reduce(vcat, [weights[c][1] for c in cells])
    V = reduce(vcat, [weights[c][2] for c in cells])
    VFSA2DMap(cells, controls, sparse(I, J, V, length(cells), length(controls)))
end

@inline function _vfsa_y(u::Float64, temperature::Float64)
    sign = u >= 0.5 ? 1.0 : -1.0
    sign * temperature * ((1 + 1 / temperature)^abs(2u - 1.0) - 1.0)
end

# Ingber's proposal on a random share of the controls, clipped to the log10 box
function _vfsa2d_propose!(p, base, T, lo, hi, n_move, rng, step_scale)
    span = (hi - lo) * step_scale
    for j in randperm(rng, length(p))[1:n_move]
        p[j] = clamp(p[j] + _vfsa_y(rand(rng), T) * span, lo - base[j], hi - base[j])
    end
    p
end

#---------- one chain ----------

# log10 model of control values p, and its data misfit (chi2 over real data, as GN)
function _vfsa2d_evaluate(problem, map::VFSA2DMap, p)
    m = copy(problem.m0)
    lo, hi = problem.config.log_bounds
    m[map.cells] .= clamp.(problem.m0[map.cells] .+ map.W * p, lo, hi)
    ρ = 10.0 .^ m
    response, _ = _mt2d_forward_cache(problem.mesh, ρ; mode = problem.config.mode, cache_fields = false)
    r = _inv2d_residual(response, problem.rows)
    chi2 = sum(abs2, r)
    (; m, ρ, response, chi2, rms = sqrt(chi2 / length(r)))
end

function _vfsa2d_chain(k::Int, problem, dir::Union{Nothing, String})
    config = problem.config
    rng = MersenneTwister(config.seed + 1000 * (k - 1))
    map = _vfsa2d_map(problem.mesh, problem.active, config, rng)
    lo, hi = config.log_bounds
    base = problem.m0[map.controls]
    p = zeros(length(map.controls))
    n_move = max(1, round(Int, config.frac_update_controls * length(p)))
    current = _vfsa2d_evaluate(problem, map, p)
    best = current
    history = [(chain = k, iteration = 0, temperature = config.temp_kappa, trial_rms = current.rms,
                rms = current.rms, best_rms = best.rms, accepted = true)]
    ak = _resolve_ak(config.cool_ratio, config.max_iter)
    every = max(1, config.max_iter ÷ 20)
    chaindir = dir === nothing ? nothing : mkpath(joinpath(dir, @sprintf("chain_%02d", k)))
    for it in 1:config.max_iter
        best.rms <= config.target_rms && break
        T = _T_schedule(it; T0 = config.temp_kappa, ak)
        # all trials branch from the current state, the best one takes one metropolis test
        trial, trial_p = nothing, p
        for _ in 1:config.n_trials
            q = _vfsa2d_propose!(copy(p), base, T, lo, hi, n_move, rng, config.step_scale)
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

Very fast simulated annealing over RBF control points, on the same problem as
`Invert2D`: only the active earth cells change (default all earth cells; air, water and
mask-0 cells stay fixed), the misfit is the error-weighted impedance chi2 of GN and NLCG.
Chains run on the Julia threads, each with its own random controls and seed. Their best
models form the ensemble: mean and median (log10), standard deviation and 5-95% range.

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
    lo, hi = config.log_bounds
    all(i -> lo <= log10(ρ0[i]) <= hi, cells) || throw(ArgumentError("active starting resistivities are outside log_bounds"))
    problem = (; mesh, rows, active, m0 = log10.(ρ0), config)
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

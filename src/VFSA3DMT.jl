# 3D VFSA MT inversion engine
# Author: @pankajkmishra
# Very Fast Simulated Annealing for MT, Gaussian-RBF parameterization, ModEM forward solver
# Includes ensemble statistics for multi-chain analysis

using Random
using Dates

#---------- configuration ----------

"""
    VFSA3DMTConfig

Configuration for the 3D VFSA MT inversion. All VFSA hyper-parameters,
file-management switches, and external-solver settings live here. Every field is
a required keyword with no default — values are supplied by the run script (see
`examples/run_vfsa3dmt.jl`), so the engine carries no hidden defaults.

`out_root::String` is the base name/path of the run directory holding the
per-chain ModEM working folders. A timestamp is appended, so the actual run
directory is `<out_root>_<yyyymmdd_HHMMSS>` (e.g. `out_root="run"` gives
`run_20260626_143000`). A relative `out_root` is created alongside the starting
model file; an absolute path is used as-is. The iteration/trial logs and the
per-chain best models are written directly inside this run directory.
"""
Base.@kwdef mutable struct VFSA3DMTConfig
    nchains::Int
    nprocs::Int
    mpirun_cmd::String
    modem_exe::String
    out_root::String
    n_ctrl::Int
    log_bounds::Tuple{Float64,Float64}
    frac_update_controls::Float64
    step_scale::Float64
    max_iter::Int
    n_trials::Int
    temp_kappa::Float64
    explore_frac::Float64
    acc_freeze_frac::Float64
    cool_ratio::Float64
    ak::Union{Nothing,Float64}
    ak_acc::Union{Nothing,Float64}
    target_rms::Float64
    seed::Int
    pad_tol::Float64
    padding_decay_length::Float64
    keep_models::Bool
    keep_dpred::Bool
    # model_save_every: 0 = behave per keep_models; >0 = retain trial models only
    # on iterations where iter % model_save_every == 0 (best model always kept)
    model_save_every::Int
    sigma_scale::Float64
    trunc_sigmas::Float64
end

#---------- internal helpers ----------

function _cleanup_modem_artifacts!(dir::AbstractString)
    isdir(dir) || return
    for f in readdir(dir; join=true)
        bn = basename(f)
        if startswith(bn, "Nodes") || startswith(bn, "BICG")
            rm(f; force=true, recursive=true)
        end
    end
end

function _find_model_files(dir::AbstractString, pattern::Regex)
    isdir(dir) || error("Directory not found: $dir")
    sort([joinpath(dir, f) for f in readdir(dir) if occursin(pattern, f)])
end

#---------- core region helpers ----------

function core_ranges(m::WS3DModel; tol::Real=0.2)
    ix = core_indices(m.cx; tol=tol)
    iy = core_indices(m.cy; tol=tol)
    return ix, iy
end

extract_core_array(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int}) = @view m.A[ix, iy, :]

function embed_core!(m::WS3DModel, A_core::AbstractArray{<:Real,3},
                     ix::UnitRange{Int}, iy::UnitRange{Int})
    @assert size(A_core,1)==length(ix) && size(A_core,2)==length(iy) && size(A_core,3)==m.nz
    @views m.A[ix, iy, :] .= A_core
    return m
end

#---------- horizontal padding decay (x,y only, z untouched) ----------

function smooth_padding_decay_xy!(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int},
                                  background_res_log10::Float64, decay_length::Float64)
    nx, ny, nz = size(m.A)
    xi1, xi2 = first(ix), last(ix)
    yi1, yi2 = first(iy), last(iy)

    dx_median = median(m.dx)
    dy_median = median(m.dy)
    L = decay_length * min(dx_median, dy_median)

    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        inside_x = (xi1 <= i <= xi2)
        inside_y = (yi1 <= j <= yi2)
        if inside_x && inside_y
            continue
        end
        ic = clamp(i, xi1, xi2)
        jc = clamp(j, yi1, yi2)
        dist = hypot((i - ic) * dx_median, (j - jc) * dy_median)
        weight = exp(-dist / max(L, eps()))
        boundary_val = m.A[ic, jc, k]
        m.A[i, j, k] = boundary_val * weight + background_res_log10 * (1 - weight)
    end
    return m
end

#---------- gaussian RBF interpolation (compact support) ----------

struct RBFMap
    nx::Int
    ny::Int
    nz::Int
    ci::Vector{Int}
    cj::Vector{Int}
    ck::Vector{Int}
    ctrl_at::Array{Int,3}
    ptr::Vector{Int}
    nbrs::Vector{Int}
    wts::Vector{Float64}
end

"""
    build_rbf_map(m, ix, iy, n_ctrl, rng; sigma_scale=2.0, trunc_sigmas=3.0)

Randomly pick up to `n_ctrl` control voxels in the core and precompute
Gaussian-RBF weights for every core voxel.
"""
function build_rbf_map(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int}, n_ctrl::Int, rng::AbstractRNG;
                       sigma_scale::Float64 = 2.0, trunc_sigmas::Float64 = 3.0)

    nx, ny, nz = length(ix), length(iy), m.nz
    N = nx * ny * nz
    n_sel = min(n_ctrl, N)

    idxs = randperm(rng, N)[1:n_sel]
    ci = Vector{Int}(undef, n_sel)
    cj = similar(ci)
    ck = similar(ci)
    @inbounds for (q, id) in enumerate(idxs)
        k = Int(ceil(id / (nx*ny)))
        rem1 = id - (k-1)*nx*ny
        j = Int(ceil(rem1 / nx))
        i = rem1 - (j-1)*nx
        ci[q] = i; cj[q] = j; ck[q] = k
    end

    ctrl_at = zeros(Int, nx, ny, nz)
    @inbounds for q in 1:n_sel
        ctrl_at[ci[q], cj[q], ck[q]] = q
    end

    dx_med = median(m.dx)
    dy_med = median(m.dy)
    dz_med = median(m.dz)

    σx = sigma_scale * dx_med
    σy = sigma_scale * dy_med
    σz = sigma_scale * dz_med

    rx = ceil(Int, trunc_sigmas * σx / dx_med)
    ry = ceil(Int, trunc_sigmas * σy / dy_med)
    rz = ceil(Int, trunc_sigmas * σz / dz_med)

    ptr  = Vector{Int}(undef, N + 1)
    nbrs = Int[];       sizehint!(nbrs, 8 * N)
    wts  = Float64[];   sizehint!(wts,  8 * N)

    row = 1
    ptr[1] = 1
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        i1 = max(1, i - rx);  i2 = min(nx, i + rx)
        j1 = max(1, j - ry);  j2 = min(ny, j + ry)
        k1 = max(1, k - rz);  k2 = min(nz, k + rz)

        local_ids = Int[]
        local_wts = Float64[]

        for kk in k1:k2, jj in j1:j2, ii in i1:i2
            q = ctrl_at[ii, jj, kk]
            if q == 0; continue; end
            Δx = (i - ii) * dx_med
            Δy = (j - jj) * dy_med
            Δz = (k - kk) * dz_med
            r2 = (Δx/σx)^2 + (Δy/σy)^2 + (Δz/σz)^2
            if r2 <= trunc_sigmas^2
                push!(local_ids, q)
                push!(local_wts, exp(-0.5 * r2))
            end
        end

        if isempty(local_ids)
            best_q = 1
            best_r2 = typemax(Float64)
            for q in 1:n_sel
                Δx = (i - ci[q]) * dx_med
                Δy = (j - cj[q]) * dy_med
                Δz = (k - ck[q]) * dz_med
                r2 = (Δx/σx)^2 + (Δy/σy)^2 + (Δz/σz)^2
                if r2 < best_r2
                    best_r2 = r2; best_q = q
                end
            end
            push!(local_ids, best_q)
            push!(local_wts, 1.0)
        end

        s = sum(local_wts)
        invs = (s > 0) ? (1.0/s) : 1.0
        for t in eachindex(local_wts)
            local_wts[t] *= invs
        end

        append!(nbrs, local_ids)
        append!(wts,  local_wts)
        ptr[row+1] = length(nbrs) + 1
        row += 1
    end

    return RBFMap(nx, ny, nz, ci, cj, ck, ctrl_at, ptr, nbrs, wts)
end

"""
    apply_rbf_map!(delta_values, rbfmap, delta_params)

Compute the core 3D field from control deltas via precomputed Gaussian-RBF weights.
"""
function apply_rbf_map!(delta_values::Array{Float64,3}, rbfmap::RBFMap, delta_params::AbstractVector{<:Real})
    @assert size(delta_values,1) == rbfmap.nx &&
            size(delta_values,2) == rbfmap.ny &&
            size(delta_values,3) == rbfmap.nz

    row = 1
    @inbounds for k in 1:rbfmap.nz, j in 1:rbfmap.ny, i in 1:rbfmap.nx
        s = 0.0
        p1 = rbfmap.ptr[row]; p2 = rbfmap.ptr[row+1] - 1
        for p in p1:p2
            q = rbfmap.nbrs[p]
            s += rbfmap.wts[p] * delta_params[q]
        end
        delta_values[i,j,k] = s
        row += 1
    end
    return delta_values
end

#---------- vfsa machinery ----------

@inline function vfsa_y(u::Float64, T::Float64)
    s = ifelse(u >= 0.5, 1.0, -1.0)
    return s * T * ((1 + 1/T)^(abs(2u - 1.0)) - 1.0)
end

# perturb nsel controls by a full-width cauchy step, resample up to 100 tries to stay in [lo,hi]
function propose_controls!(delta_params::Vector{Float64}, T::Float64,
                           lo::Float64, hi::Float64, v0_at_ctrl::Vector{Float64},
                           nsel::Int, rng::AbstractRNG;
                           step_scale::Float64=1.0)
    M = length(delta_params)
    idxs = randperm(rng, M)[1:nsel]
    dm = (hi - lo) * step_scale
    mid = 0.5 * (lo + hi)
    @inbounds for id in idxs
        base = v0_at_ctrl[id] + delta_params[id]
        cand = mid
        for _ in 1:100
            c = base + vfsa_y(rand(rng), T) * dm
            if lo <= c <= hi
                cand = c
                break
            end
        end
        delta_params[id] = cand - v0_at_ctrl[id]
    end
    return idxs
end

function forward_and_misfit!(m::WS3DModel;
                             run_dir::AbstractString,
                             model_filename::String,
                             dobs_filename::String,
                             dpred_filename::String,
                             origin::Vector{Float64},
                             rotation::Float64,
                             cfg::VFSA3DMTConfig)
    model_abs = joinpath(run_dir, model_filename)
    dpred_abs = joinpath(run_dir, dpred_filename)
    dobs_abs  = joinpath(run_dir, dobs_filename)
    write_ws3d_model(model_abs, m.dx, m.dy, m.dz, m.A, origin; rotation=rotation, type_str="LOGE")
    ok = true
    cd(run_dir) do
        cmd = `$(cfg.mpirun_cmd) -n $(cfg.nprocs) $(cfg.modem_exe) -F $(model_filename) $(dobs_filename) $(dpred_filename)`
        try
            run(pipeline(cmd, stdout="0runlog.dat"))
        catch e
            @warn("[forward] ModEM failed", error=e)
            ok = false
        end
    end
    if !ok || !isfile(dpred_abs)
        _cleanup_modem_artifacts!(run_dir)
        return Inf
    end
    s = chi2_and_rms(dobs_abs, dpred_abs; use_impedance=true, use_tipper=true, components=String[])
    _cleanup_modem_artifacts!(run_dir)
    return s.rms
end

#---------- logging ----------

function _write_trials_header_3d(path::AbstractString; timestamp::String, cfg::VFSA3DMTConfig, chain_seed::Int, T0_chain::Float64, ak_prop::Float64, ak_acc::Float64)
    open(path, "w") do io
        println(io, "# VFSA 3D MT detailed trials — ", timestamp)
        println(io, "# chains=", cfg.nchains,
                    "  seed=", chain_seed,
                    "  n_ctrl=", cfg.n_ctrl,
                    "  frac_update_controls=", cfg.frac_update_controls,
                    "  temp_kappa=", cfg.temp_kappa,
                    "  T0=", T0_chain,
                    "  explore_frac=", cfg.explore_frac,
                    "  acc_freeze_frac=", cfg.acc_freeze_frac,
                    "  cool_ratio=", cfg.cool_ratio,
                    "  ak_prop=", ak_prop,
                    "  ak_acc=", ak_acc,
                    "  max_iter=", cfg.max_iter,
                    "  target_rms=", cfg.target_rms,
                    "  n_trials=", cfg.n_trials)
        println(io, repeat("-", 114))
        @printf(io, "%8s %8s %12s %12s %11s %12s %10s %10s %5s %11s %s\n",
                "Iter","Trial","Tprop","Tacc",
                "RMS","dE_rms2","Pacc","Uacc","Acc","RMSBest","Model")
        println(io, repeat("-", 114))
    end
end

function _append_trial_row_3d(path::AbstractString; iter::Int, trial::Int,
                              Tprop::Float64, Tacc::Float64,
                           rms_prop::Float64,
                           dE::Float64, p_acc::Float64, u_acc::Float64,
                           rms_best::Float64,
                           acc::Int, model_rel::String)
    # one row per trial: the proposal's own misfit, energy dE, metropolis
    # prob/draw, accept flag, running best, and the proposal model filename
    open(path, "a") do io
        @printf(io, "%8d %8d %12.4g %12.4g %11.5f %12.5f %10.5f %10.5f %5d %11.5f %s\n",
                iter, trial, Tprop, Tacc, rms_prop, dE, p_acc, u_acc, acc, rms_best, model_rel)
    end
end

function _write_iter_header_3d(path::AbstractString; timestamp::String, cfg::VFSA3DMTConfig, chain_seed::Int, T0_chain::Float64, ak_prop::Float64, ak_acc::Float64)
    open(path, "w") do io
        println(io, "# VFSA 3D MT iteration best — ", timestamp)
        println(io, "# chains=", cfg.nchains,
                    "  seed=", chain_seed,
                    "  n_ctrl=", cfg.n_ctrl,
                    "  frac_update_controls=", cfg.frac_update_controls,
                    "  temp_kappa=", cfg.temp_kappa,
                    "  T0=", T0_chain,
                    "  explore_frac=", cfg.explore_frac,
                    "  acc_freeze_frac=", cfg.acc_freeze_frac,
                    "  cool_ratio=", cfg.cool_ratio,
                    "  ak_prop=", ak_prop,
                    "  ak_acc=", ak_acc,
                    "  max_iter=", cfg.max_iter,
                    "  target_rms=", cfg.target_rms,
                    "  n_trials=", cfg.n_trials)
        println(io, repeat("-", 100))
        # current accepted state after the iteration, Nacc = accepted trials this iter
        @printf(io, "%8s %12s %12s %11s %11s %6s %s\n",
                "Iter","Tprop","Tacc","RMS","RMSBest","Nacc","Model")
        println(io, repeat("-", 100))
    end
end

function _append_iter_row_3d(path::AbstractString; iter::Int,
                          Tprop::Float64, Tacc::Float64,
                          rms_curr::Float64,
                          rms_best::Float64,
                          nacc::Int, model_rel::String)
    open(path, "a") do io
        @printf(io, "%8d %12.4g %12.4g %11.5f %11.5f %6d %s\n",
                iter, Tprop, Tacc, rms_curr, rms_best, nacc, model_rel)
    end
end

_T_schedule(k::Int; T0::Float64, ak::Float64) = T0 * exp(-sqrt(ak) * (k - 1))

function _resolve_ak(cfg::VFSA3DMTConfig, freeze_frac::Float64, override::Union{Nothing,Float64})
    override === nothing || return override
    k_freeze = freeze_frac * cfg.max_iter
    return (log(1.0 / cfg.cool_ratio) / max(k_freeze - 1.0, 1.0))^2
end

#---------- single-chain vfsa loop ----------

function _run_chain_3d(chain_id::Int, start_model_path::String,
                       dobs_filename::String, timestamp::String,
                       run_root::String, results_root::String; cfg::VFSA3DMTConfig)
    chain_seed = cfg.seed + 1000*(chain_id-1)
    rng = MersenneTwister(chain_seed)
    chain_dir = joinpath(run_root, @sprintf("chain_%02d", chain_id))
    mkpath(chain_dir)

    # each chain keeps its own logs so chains stay independent and parallel-safe
    trials_log = joinpath(chain_dir, "0vfsa3DMT_detailed.log")
    iter_log   = joinpath(chain_dir, "0vfsa3DMT.log")

    dobs_chain_path = joinpath(chain_dir, dobs_filename)
    if !isfile(dobs_chain_path)
        cp(joinpath(run_root, dobs_filename), dobs_chain_path; force=true)
    end

    m = load_ws3d_model(start_model_path)
    _, _, _, _, _, _, origin, rotation = read_ws3d_model(start_model_path, true)

    ix, iy = core_ranges(m; tol=cfg.pad_tol)
    Acore = extract_core_array(m, ix, iy)
    v0_core_log10 = Array(Acore)

    #---------- background resistivity from the padding ring ----------
    nx_full, ny_full, _ = size(m.A)
    background_resistivities = Float64[]
    for i in 1:nx_full
        if i < first(ix) || i > last(ix)
            append!(background_resistivities, vec(m.A[i, :, :]))
        end
    end
    for j in 1:ny_full
        if j < first(iy) || j > last(iy)
            append!(background_resistivities, vec(m.A[:, j, :]))
        end
    end
    finite_bg = filter(isfinite, background_resistivities)
    background_log10 = median(finite_bg)

    # build the RBF map
    rbfmap = build_rbf_map(m, ix, iy, cfg.n_ctrl, rng;
                           sigma_scale=cfg.sigma_scale, trunc_sigmas=cfg.trunc_sigmas)
    M = length(rbfmap.ci)

    v0_ctrl = Vector{Float64}(undef, M)
    @inbounds for t in 1:M
        v0_ctrl[t] = v0_core_log10[rbfmap.ci[t], rbfmap.cj[t], rbfmap.ck[t]]
    end

    delta_params_current = zeros(Float64, M)
    delta_values_buffer = zeros(Float64, rbfmap.nx, rbfmap.ny, rbfmap.nz)
    nsel_default = max(1, round(Int, cfg.frac_update_controls * M))

    # initial forward on the start model
    model0_filename = @sprintf("model_%02d_%03d_%02d.rho", chain_id, 0, 0)
    dpred0_filename = @sprintf("dpred_%02d_%03d_%02d.dat", chain_id, 0, 0)

    embed_core!(m, v0_core_log10, ix, iy)
    smooth_padding_decay_xy!(m, ix, iy, background_log10, cfg.padding_decay_length)

    dp0 = forward_and_misfit!(m; run_dir=chain_dir, model_filename=model0_filename,
                              dobs_filename=dobs_filename, dpred_filename=dpred0_filename,
                              origin=origin, rotation=rotation, cfg=cfg)
    rms_current = dp0
    best_rms    = rms_current
    T0_chain = cfg.temp_kappa
    ak_prop_chain = _resolve_ak(cfg, cfg.explore_frac, cfg.ak)
    ak_acc_chain = _resolve_ak(cfg, cfg.acc_freeze_frac, cfg.ak_acc)
    _write_trials_header_3d(trials_log; timestamp=timestamp, cfg=cfg, chain_seed=chain_seed, T0_chain=T0_chain, ak_prop=ak_prop_chain, ak_acc=ak_acc_chain)
    _write_iter_header_3d(iter_log; timestamp=timestamp, cfg=cfg, chain_seed=chain_seed, T0_chain=T0_chain, ak_prop=ak_prop_chain, ak_acc=ak_acc_chain)
    best_model_abs = joinpath(results_root, @sprintf("best_model_chain%02d.rho", chain_id))
    cp(joinpath(chain_dir, model0_filename), best_model_abs; force=true)
    # filename of the current accepted model, carried forward so the iteration
    # log always names a model even on iterations with no accepted trial
    current_model_rel = model0_filename

    lo, hi = cfg.log_bounds

    for k in 1:cfg.max_iter
        T_prop = _T_schedule(k; T0=T0_chain, ak=ak_prop_chain)
        T_acc = _T_schedule(k; T0=T0_chain, ak=ak_acc_chain)

        #---------- per-trial metropolis against the current state ----------
        n_accepted_iter = 0
        last_accepted_trial = 0
        trial_cache = NamedTuple[]

        for t in 1:cfg.n_trials
            # propose from the current accepted state
            delta_params_trial = copy(delta_params_current)
            propose_controls!(delta_params_trial, T_prop, lo, hi, v0_ctrl, nsel_default, rng;
                              step_scale=cfg.step_scale)

            apply_rbf_map!(delta_values_buffer, rbfmap, delta_params_trial)
            v_trial = clamp.(v0_core_log10 .+ delta_values_buffer, lo, hi)

            embed_core!(m, v_trial, ix, iy)
            smooth_padding_decay_xy!(m, ix, iy, background_log10, cfg.padding_decay_length)

            model_filename = @sprintf("model_%02d_%03d_%02d.rho", chain_id, k, t)
            dpred_filename = @sprintf("dpred_%02d_%03d_%02d.dat", chain_id, k, t)
            dp = forward_and_misfit!(m; run_dir=chain_dir, model_filename=model_filename,
                                     dobs_filename=dobs_filename, dpred_filename=dpred_filename,
                                     origin=origin, rotation=rotation, cfg=cfg)

            dE = (dp^2 - rms_current^2) / max(rms_current^2, eps())
            u_acc = rand(rng)
            # Pacc = 1 downhill, exp(-dE/T_acc) uphill, kept finite for the log
            p_acc = dE <= 0 ? 1.0 : exp(-dE / max(T_acc, 1e-12))
            accept_trial = isfinite(dp) && (u_acc < p_acc)

            push!(trial_cache, (
                iter=k, trial=t, Tprop=T_prop, Tacc=T_acc,
                rms_prop=dp,
                dE=dE, p_acc=p_acc, u_acc=u_acc, accepted=accept_trial,
                rms_best=best_rms,
                model_rel=model_filename
            ))

            if accept_trial
                # advance the chain, next trial starts here
                delta_params_current .= delta_params_trial
                rms_current  = dp
                n_accepted_iter += 1
                last_accepted_trial = t
                current_model_rel = model_filename
                # record the global best on improvement, never feed it back
                if rms_current < best_rms
                    best_rms  = rms_current
                    cp(joinpath(chain_dir, model_filename), best_model_abs; force=true)
                end
            end

            # legacy per-trial cleanup; with model_save_every>0 the whole
            # iteration is pruned below instead
            if cfg.model_save_every <= 0
                if !cfg.keep_models && t > 1
                    prev_model = joinpath(chain_dir, @sprintf("model_%02d_%03d_%02d.rho", chain_id, k, t-1))
                    isfile(prev_model) && rm(prev_model; force=true)
                end
                if !cfg.keep_dpred && t > 1
                    prev_dpred = joinpath(chain_dir, @sprintf("dpred_%02d_%03d_%02d.dat", chain_id, k, t-1))
                    isfile(prev_dpred) && rm(prev_dpred; force=true)
                end
            end
        end

        #---------- checkpointed pruning (model_save_every > 0) ----------
        # keep all trial files only on checkpoint iters; on a checkpoint keep
        # just the last accepted trial; otherwise drop the whole iteration.
        # the best model lives in best_model_abs and is never touched here.
        if cfg.model_save_every > 0
            is_checkpoint = (k % cfg.model_save_every == 0)
            keep_trial = last_accepted_trial > 0 ? last_accepted_trial : cfg.n_trials
            for t in 1:cfg.n_trials
                if is_checkpoint && t == keep_trial
                    continue
                end
                model_t = joinpath(chain_dir, @sprintf("model_%02d_%03d_%02d.rho", chain_id, k, t))
                isfile(model_t) && rm(model_t; force=true)
                dpred_t = joinpath(chain_dir, @sprintf("dpred_%02d_%03d_%02d.dat", chain_id, k, t))
                isfile(dpred_t) && rm(dpred_t; force=true)
            end
        end

        # write one row per trial
        for tr in trial_cache
            _append_trial_row_3d(trials_log; iter=tr.iter, trial=tr.trial,
                              Tprop=tr.Tprop, Tacc=tr.Tacc,
                              rms_prop=tr.rms_prop,
                              dE=tr.dE, p_acc=tr.p_acc, u_acc=tr.u_acc,
                              rms_best=tr.rms_best,
                              acc=(tr.accepted ? 1 : 0), model_rel=tr.model_rel)
        end

        # iteration summary row: current accepted state, named by the model that
        # represents it (carried forward when no trial was accepted this iter)
        _append_iter_row_3d(iter_log; iter=k,
                         Tprop=T_prop, Tacc=T_acc,
                         rms_curr=rms_current,
                         rms_best=best_rms,
                         nacc=n_accepted_iter, model_rel=current_model_rel)

        # stop once best rms hits the target
        if best_rms <= cfg.target_rms
            @info @sprintf("[chain %02d] target RMS %.3f reached at iter %d (best RMS %.5f); stopping.",
                           chain_id, cfg.target_rms, k, best_rms)
            break
        end
    end
    return nothing
end

#---------- main entry point ----------

"""
    VFSA3DMT(start_model_path; dobs_path, cfg)

Run a 3D VFSA MT inversion. Returns `(best_model_path, iter_log_path)`.

# Arguments
- `start_model_path::AbstractString`: Path to a WS3D-format starting model (.rho).
- `dobs_path::AbstractString`: Path to the ModEM observed data file.
- `cfg::VFSA3DMTConfig`: Configuration struct (required; see `VFSA3DMTConfig`).
"""
function VFSA3DMT(start_model_path::AbstractString;
                  dobs_path::AbstractString,
                  cfg::VFSA3DMTConfig)
    start_model_abs = abspath(start_model_path)
    model_dir = dirname(start_model_abs)

    run_root_base = isabspath(cfg.out_root) ? cfg.out_root : joinpath(model_dir, cfg.out_root)
    t_now = now()
    run_root = string(run_root_base, "_", Dates.format(t_now, "yyyymmdd_HHMMSS"))
    isdir(run_root) || mkpath(run_root)

    dobs_filename = basename(dobs_path)
    dobs_abs_target = abspath(joinpath(run_root, dobs_filename))
    cp(abspath(dobs_path), dobs_abs_target; force=true)

    # per-chain logs live in chain_NN/, best models flat in run_root for the ensemble readers
    timestamp = Dates.format(t_now, "yyyy-mm-dd HH:MM:SS")

    for c in 1:cfg.nchains
        _run_chain_3d(c, start_model_abs, dobs_filename, timestamp,
                      run_root, run_root; cfg=cfg)
    end

    best_model = joinpath(run_root, @sprintf("best_model_chain%02d.rho", 1))
    iter_log   = joinpath(run_root, @sprintf("chain_%02d", 1), "0vfsa3DMT.log")
    return best_model, iter_log
end

#---------- ensemble statistics ----------

function _intersect_ranges(r1::UnitRange{Int}, r2::UnitRange{Int})
    a = max(first(r1), first(r2)); b = min(last(r1), last(r2))
    b < a && error("Empty intersection across chains.")
    a:b
end

function _maxabsdiff(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    length(a) == length(b) || return Inf
    d = 0.0
    @inbounds for i in eachindex(a, b)
        v = abs(a[i] - b[i])
        v > d && (d = v)
    end
    d
end

function _assert_grid_match!(mref, m, ix, iy, iz; grid_atol::Float64, grid_rtol::Float64)
    dxr = mref.dx[ix]; dyr = mref.dy[iy]; dzr = mref.dz[iz]
    dxm = m.dx[ix];    dym = m.dy[iy];    dzm = m.dz[iz]
    sdx = _maxabsdiff(dxr, dxm)
    sdy = _maxabsdiff(dyr, dym)
    sdz = _maxabsdiff(dzr, dzm)
    if !(sdx <= grid_atol + grid_rtol * maximum(abs.(dxr)) &&
         sdy <= grid_atol + grid_rtol * maximum(abs.(dyr)) &&
         sdz <= grid_atol + grid_rtol * maximum(abs.(dzr)))
        error("Grid mismatch across chains within the overlapping index ranges.")
    end
    if length(mref.origin) >= 3 && length(m.origin) >= 3
        dor = maximum(abs.(mref.origin[1:3] .- m.origin[1:3]))
        dor <= 1e-6 || error("Origin mismatch across chains.")
    end
    nothing
end

"""
    core_statistics(cores::Vector{Array{Float64,3}})

Compute element-wise mean, median, std over a vector of 3D arrays.
"""
function core_statistics(cores::Vector{Array{Float64,3}})
    nx, ny, nz = size(cores[1])
    meanA   = Array{Float64}(undef, nx, ny, nz)
    medianA = Array{Float64}(undef, nx, ny, nz)
    stdA    = Array{Float64}(undef, nx, ny, nz)
    vals = Vector{Float64}(undef, length(cores))

    for k in 1:nz, j in 1:ny, i in 1:nx
        cnt = 0
        for c in 1:length(cores)
            v = cores[c][i,j,k]
            if isfinite(v)
                cnt += 1
                vals[cnt] = v
            end
        end
        if cnt == 0
            meanA[i,j,k] = NaN
            medianA[i,j,k] = NaN
            stdA[i,j,k] = NaN
        else
            μ = sum(vals[1:cnt]) / cnt
            meanA[i,j,k] = μ
            sort!(vals[1:cnt])
            if isodd(cnt)
                medianA[i,j,k] = vals[(cnt + 1) >>> 1]
            else
                medianA[i,j,k] = 0.5 * (vals[cnt >>> 1] + vals[(cnt >>> 1) + 1])
            end
            stdA[i,j,k] = cnt == 1 ? 0.0 : sqrt(sum((vals[1:cnt] .- μ).^2) / (cnt - 1))
        end
    end
    meanA, medianA, stdA
end

"""
    AnalyseEnsemble3D(input_dir; pattern=r"^best_model.*\\.rho\$",
                      grid_atol=1e-8, grid_rtol=1e-8)

Load all matching model files from `input_dir`, validate grid consistency,
compute ensemble mean/median/std, and write them as WS3D model files. The default
pattern matches the per-chain best models (`best_model_chainNN.rho`) that VFSA3DMT
writes into the starting model's directory.

Returns `(mean_path, median_path, std_path)`.
"""
function AnalyseEnsemble3D(input_dir::AbstractString;
                           pattern::Regex = r"^best_model.*\.rho$",
                           grid_atol::Float64 = 1e-8,
                           grid_rtol::Float64 = 1e-8)
    files = _find_model_files(input_dir, pattern)
    isempty(files) && error("No model files matching pattern in $input_dir")

    models = map(load_ws3d_model, files)
    mref = models[1]

    ix = 1:size(mref.A, 1)
    iy = 1:size(mref.A, 2)
    iz = 1:size(mref.A, 3)
    for m in models[2:end]
        ix = _intersect_ranges(ix, 1:size(m.A, 1))
        iy = _intersect_ranges(iy, 1:size(m.A, 2))
        iz = _intersect_ranges(iz, 1:size(m.A, 3))
    end

    for m in models[2:end]
        _assert_grid_match!(mref, m, ix, iy, iz; grid_atol=grid_atol, grid_rtol=grid_rtol)
    end

    cubes = [m.A[ix, iy, iz] for m in models]
    meanA, medianA, stdA = core_statistics(cubes)

    dx = mref.dx[ix]
    dy = mref.dy[iy]
    dz = mref.dz[iz]

    origin = copy(mref.origin)
    if length(origin) >= 1 && first(ix) > 1
        origin[1] += sum(mref.dx[1:first(ix)-1])
    end
    if length(origin) >= 2 && first(iy) > 1
        origin[2] += sum(mref.dy[1:first(iy)-1])
    end
    if length(origin) >= 3 && first(iz) > 1
        origin[3] += sum(mref.dz[1:first(iz)-1])
    end

    mean_path   = joinpath(input_dir, "model.mean")
    median_path = joinpath(input_dir, "model.median")
    std_path    = joinpath(input_dir, "model.std")

    write_ws3d_model(mean_path,   dx, dy, dz, meanA,   origin; rotation=0.0, type_str="LOGE")
    write_ws3d_model(median_path, dx, dy, dz, medianA, origin; rotation=0.0, type_str="LOGE")
    write_ws3d_model(std_path,    dx, dy, dz, stdA,    origin; rotation=0.0, type_str="LOGE")

    println("Wrote: $mean_path")
    println("Wrote: $median_path")
    println("Wrote: $std_path")

    return mean_path, median_path, std_path
end

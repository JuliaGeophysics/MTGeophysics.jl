# 3D VFSA MT Inversion Engine.
# Author: @pankajkmishra
# This file provides the full 3D Very Fast Simulated Annealing (VFSA) inversion
# workflow for magnetotelluric data, using Gaussian-RBF parameterization and
# ModEM as the external forward solver. It also includes ensemble statistics
# for multi-chain analysis.

using Random
using Dates

# =========================================================================== #
#  Configuration
# =========================================================================== #

"""
    VFSA3DMTConfig

Configuration for the 3D VFSA MT inversion. All VFSA hyper-parameters,
file-management switches, and external-solver settings live here.
"""
Base.@kwdef mutable struct VFSA3DMTConfig
    nchains::Int                  = 1
    nprocs::Int                   = 21
    mpirun_cmd::String            = "mpirun"
    modem_exe::String             = "Mod3DMT_2025"
    out_root::String              = ""
    n_ctrl::Int                   = 900
    log_bounds::Tuple{Float64,Float64} = (0.0, 5.0)
    frac_update_controls::Float64 = 1.0
    step_scale::Float64           = 0.05
    max_iter::Int                 = 3000
    n_trials::Int                 = 4
    T0_prop::Float64              = 1.0
    Tf_prop::Float64              = 1e-3
    T0_acc::Float64               = 1.0
    Tf_acc::Float64               = 1e-3
    seed::Int                     = 1911
    pad_tol::Float64              = 0.2
    padding_decay_length::Float64 = 10.0
    keep_models::Bool             = true
    keep_dpred::Bool              = false
    sigma_scale::Float64          = 2.0
    trunc_sigmas::Float64         = 3.0
end

# =========================================================================== #
#  Internal helpers
# =========================================================================== #

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

# =========================================================================== #
#  Core region helpers (WS3DModel-specific)
# =========================================================================== #

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

# =========================================================================== #
#  Horizontal-only padding decay (x,y only; z untouched)
# =========================================================================== #

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

# =========================================================================== #
#  Gaussian RBF interpolation (compact support)
# =========================================================================== #

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

# =========================================================================== #
#  VFSA machinery
# =========================================================================== #

@inline function vfsa_y(u::Float64, T::Float64)
    s = ifelse(u >= 0.5, 1.0, -1.0)
    return s * T * ((1 + 1/T)^(abs(2u - 1.0)) - 1.0)
end

function propose_controls!(delta_params::Vector{Float64}, T::Float64,
                           lo::Float64, hi::Float64, v0_at_ctrl::Vector{Float64},
                           nsel::Int, rng::AbstractRNG;
                           step_scale::Float64=0.05)
    M = length(delta_params)
    idxs = randperm(rng, M)[1:nsel]
    step_size = (hi - lo) * step_scale
    @inbounds for id in idxs
        u = rand(rng)
        y = vfsa_y(u, T)
        trial = delta_params[id] + y * step_size
        v_abs = v0_at_ctrl[id] + trial
        if v_abs < lo
            trial = lo - v0_at_ctrl[id]
        elseif v_abs > hi
            trial = hi - v0_at_ctrl[id]
        end
        delta_params[id] = trial
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
        return (chi2=Inf, rms=Inf)
    end
    s = chi2_and_rms(dobs_abs, dpred_abs; use_impedance=true, use_tipper=true, components=String[])
    _cleanup_modem_artifacts!(run_dir)
    return (chi2=s.χ², rms=s.rms)
end

# =========================================================================== #
#  Logging
# =========================================================================== #

function _write_trials_header_3d(path::AbstractString; timestamp::String, cfg::VFSA3DMTConfig)
    open(path, "w") do io
        println(io, "# VFSA 3D MT detailed trials — ", timestamp)
        println(io, "# chains=", cfg.nchains,
                    "  n_ctrl=", cfg.n_ctrl,
                    "  frac_update_controls=", cfg.frac_update_controls,
                    "  Tproposal=", cfg.T0_prop, "→", cfg.Tf_prop,
                    "  Tacceptance=", cfg.T0_acc, "→", cfg.Tf_acc,
                    "  n_trials=", cfg.n_trials)
        println(io, repeat("-", 230))
        @printf(io, "%6s %8s %8s %12s %14s %14s %11s %14s %12s %12s %10s %12s %11s %8s %s\n",
                "Chain","Iter","Trial","Tproposal","Tacceptance","Chi2_before","RMS_before",
                "Chi2Proposal","dChi2","RMSProposal","dRMS","Chi2Best","RMSBest","AccIdx","Model")
        println(io, repeat("-", 230))
    end
end

function _append_trial_row_3d(path::AbstractString; chain::Int, iter::Int, trial::Int,
                              Tprop::Float64, Tacc::Float64,
                           chi2_before::Float64, rms_before::Float64,
                           chi2_prop::Float64, dchi2::Float64,
                           rms_prop::Float64, drms::Float64,
                           chi2_best::Float64, rms_best::Float64,
                           acc::Int, model_rel::String)
    open(path, "a") do io
        @printf(io, "%6d %8d %8d %12.4g %14.4g %14.4f %11.5f %14.4f %12.4f %12.5f %10.5f %12.4f %11.5f %8d %s\n",
                chain, iter, trial, Tprop, Tacc, chi2_before, rms_before, chi2_prop, dchi2, rms_prop, drms, chi2_best, rms_best, acc, model_rel)
    end
end

function _write_iter_header_3d(path::AbstractString; timestamp::String, cfg::VFSA3DMTConfig)
    open(path, "w") do io
        println(io, "# VFSA 3D MT iteration best — ", timestamp)
        println(io, "# chains=", cfg.nchains,
                    "  n_ctrl=", cfg.n_ctrl,
                    "  frac_update_controls=", cfg.frac_update_controls,
                    "  Tproposal=", cfg.T0_prop, "→", cfg.Tf_prop,
                    "  Tacceptance=", cfg.T0_acc, "→", cfg.Tf_acc,
                    "  n_trials=", cfg.n_trials)
        println(io, repeat("-", 200))
        @printf(io, "%6s %8s %12s %14s %12s %12s %11s %12s %11s %8s %s\n",
                "Chain","Iter","Tproposal","Chi2Proposal","dChi2","RMSProposal","dRMS","Chi2Best","RMSBest","AccIdx","Model")
        println(io, repeat("-", 200))
    end
end

function _append_iter_row_3d(path::AbstractString; chain::Int, iter::Int,
                          Tprop::Float64,
                          chi2_prop::Float64, dchi2::Float64,
                          rms_prop::Float64, drms::Float64,
                          chi2_best::Float64, rms_best::Float64,
                          acc::Int, model_rel::String)
    open(path, "a") do io
        @printf(io, "%6d %8d %12.4g %14.4f %12.4f %12.5f %11.5f %14.4f %11.5f %8d %s\n",
                chain, iter, Tprop, chi2_prop, dchi2, rms_prop, drms, chi2_best, rms_best, acc, model_rel)
    end
end

_T_schedule(k::Int; T0::Float64, Tend::Float64, Nref::Int) =
    T0 * (Tend/T0)^((k-1)/max(Nref-1,1))

# =========================================================================== #
#  Single-chain VFSA loop
# =========================================================================== #

function _run_chain_3d(chain_id::Int, start_model_path::String,
                       dobs_filename::String, timestamp::String,
                       trials_log::String, iter_log::String,
                       run_root::String; cfg::VFSA3DMTConfig)
    rng = MersenneTwister(cfg.seed + 1000*(chain_id-1))
    chain_dir = joinpath(run_root, @sprintf("chain_%02d", chain_id))
    mkpath(chain_dir)

    dobs_chain_path = joinpath(chain_dir, dobs_filename)
    if !isfile(dobs_chain_path)
        cp(joinpath(run_root, dobs_filename), dobs_chain_path; force=true)
    end

    m = load_ws3d_model(start_model_path)
    _, _, _, _, _, _, origin, rotation = read_ws3d_model(start_model_path, true)

    ix, iy = core_ranges(m; tol=cfg.pad_tol)
    Acore = extract_core_array(m, ix, iy)
    v0_core_log10 = Array(Acore)

    # background resistivity from padding
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

    # build RBF
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

    # initial forward
    model0_filename = @sprintf("model_%02d_%03d_%02d.rho", chain_id, 0, 0)
    dpred0_filename = @sprintf("dpred_%02d_%03d_%02d.dat", chain_id, 0, 0)

    embed_core!(m, v0_core_log10, ix, iy)
    smooth_padding_decay_xy!(m, ix, iy, background_log10, cfg.padding_decay_length)

    dp0 = forward_and_misfit!(m; run_dir=chain_dir, model_filename=model0_filename,
                              dobs_filename=dobs_filename, dpred_filename=dpred0_filename,
                              origin=origin, rotation=rotation, cfg=cfg)
    chi2_current = dp0.chi2
    rms_current  = dp0.rms
    best_chi2 = chi2_current
    best_rms  = rms_current
    best_model_abs = joinpath(run_root, @sprintf("best_model_chain%02d.rho", chain_id))
    cp(joinpath(chain_dir, model0_filename), best_model_abs; force=true)

    lo, hi = cfg.log_bounds

    for k in 1:cfg.max_iter
        Tprop = _T_schedule(k; T0=cfg.T0_prop, Tend=cfg.Tf_prop, Nref=cfg.max_iter)
        Tacc  = _T_schedule(k; T0=cfg.T0_acc,  Tend=cfg.Tf_acc,  Nref=cfg.max_iter)

        chi2_before_iter = chi2_current
        rms_before_iter  = rms_current

        best_trial_chi2 = Inf
        best_trial_rms  = Inf
        best_trial_delta_params = nothing
        best_trial_model_filename = ""
        best_trial_idx = 0
        trial_cache = NamedTuple[]

        for t in 1:cfg.n_trials
            delta_params_trial = copy(delta_params_current)
            propose_controls!(delta_params_trial, Tprop, lo, hi, v0_ctrl, nsel_default, rng;
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

            delta_chi2 = dp.chi2 - chi2_before_iter
            delta_rms  = dp.rms - rms_before_iter

            push!(trial_cache, (
                chain=chain_id, iter=k, trial=t, Tprop=Tprop, Tacc=Tacc,
                chi2_before=chi2_before_iter, rms_before=rms_before_iter,
                chi2_prop=dp.chi2, dchi2=delta_chi2,
                rms_prop=dp.rms, drms=delta_rms,
                chi2_best=best_chi2, rms_best=best_rms,
                model_rel=model_filename
            ))

            if dp.chi2 < best_trial_chi2
                best_trial_chi2 = dp.chi2
                best_trial_rms  = dp.rms
                best_trial_delta_params = copy(delta_params_trial)
                best_trial_model_filename = model_filename
                best_trial_idx = t
            end

            if !cfg.keep_models && t > 1
                prev_model = joinpath(chain_dir, @sprintf("model_%02d_%03d_%02d.rho", chain_id, k, t-1))
                isfile(prev_model) && rm(prev_model; force=true)
            end
            if !cfg.keep_dpred && t > 1
                prev_dpred = joinpath(chain_dir, @sprintf("dpred_%02d_%03d_%02d.dat", chain_id, k, t-1))
                isfile(prev_dpred) && rm(prev_dpred; force=true)
            end
        end

        # Metropolis-Hastings acceptance
        delta_chi2_iter = best_trial_chi2 - chi2_before_iter
        delta_rel = delta_chi2_iter / max(chi2_before_iter, eps(Float64))
        accept = isfinite(best_trial_chi2) &&
                 ((delta_rel <= 0) || (rand(rng) < exp(-delta_rel / max(Tacc, 1e-12))))
        acc_idx = accept ? 1 : 0

        if accept
            delta_params_current .= best_trial_delta_params
            chi2_current = best_trial_chi2
            rms_current  = best_trial_rms
            if chi2_current < best_chi2
                best_chi2 = chi2_current
                best_rms  = best_trial_rms
                cp(joinpath(chain_dir, best_trial_model_filename), best_model_abs; force=true)
            end
        end

        # write logs
        for tr in trial_cache
            acc_flag = (tr.trial == best_trial_idx) ? acc_idx : 0
            _append_trial_row_3d(trials_log; chain=tr.chain, iter=tr.iter, trial=tr.trial,
                              Tprop=tr.Tprop, Tacc=tr.Tacc,
                              chi2_before=tr.chi2_before, rms_before=tr.rms_before,
                              chi2_prop=tr.chi2_prop, dchi2=tr.dchi2,
                              rms_prop=tr.rms_prop, drms=tr.drms,
                              chi2_best=tr.chi2_best, rms_best=tr.rms_best,
                              acc=acc_flag, model_rel=tr.model_rel)
        end

        _append_iter_row_3d(iter_log; chain=chain_id, iter=k,
                         Tprop=Tprop,
                         chi2_prop=best_trial_chi2, dchi2=delta_chi2_iter,
                         rms_prop=best_trial_rms, drms=(best_trial_rms - rms_before_iter),
                         chi2_best=best_chi2, rms_best=best_rms,
                         acc=acc_idx, model_rel=best_trial_model_filename)
    end
    return nothing
end

# =========================================================================== #
#  Main entry point
# =========================================================================== #

"""
    VFSA3DMT(start_model_path; dobs_path, cfg=VFSA3DMTConfig())

Run a 3D VFSA MT inversion. Returns `(best_model_path, iter_log_path)`.

# Arguments
- `start_model_path::AbstractString`: Path to a WS3D-format starting model (.rho).
- `dobs_path::AbstractString`: Path to the ModEM observed data file.
- `cfg::VFSA3DMTConfig`: Configuration struct (see `VFSA3DMTConfig`).
"""
function VFSA3DMT(start_model_path::AbstractString;
                  dobs_path::AbstractString,
                  cfg::VFSA3DMTConfig = VFSA3DMTConfig())
    if isempty(cfg.out_root)
        cfg.out_root = joinpath("runs", Dates.format(now(), "yyyymmdd_HHMMSS"))
    end
    isdir(cfg.out_root) || mkpath(cfg.out_root)

    dobs_filename = basename(dobs_path)
    dobs_abs_target = abspath(joinpath(cfg.out_root, dobs_filename))
    cp(abspath(dobs_path), dobs_abs_target; force=true)

    trials_log = joinpath(cfg.out_root, "0vfsa3DMT_detailed.log")
    iter_log   = joinpath(cfg.out_root, "0vfsa3DMT.log")
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    _write_trials_header_3d(trials_log; timestamp=timestamp, cfg=cfg)
    _write_iter_header_3d(iter_log; timestamp=timestamp, cfg=cfg)

    for c in 1:cfg.nchains
        _run_chain_3d(c, abspath(start_model_path), dobs_filename, timestamp,
                      trials_log, iter_log, cfg.out_root; cfg=cfg)
    end

    best_model = joinpath(cfg.out_root, @sprintf("best_model_chain%02d.rho", 1))
    return best_model, iter_log
end

# =========================================================================== #
#  Ensemble statistics for multi-chain analysis
# =========================================================================== #

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
    AnalyseEnsemble3D(input_dir; pattern=r"^chain.*\\.rho\$",
                      grid_atol=1e-8, grid_rtol=1e-8)

Load all matching model files from `input_dir`, validate grid consistency,
compute ensemble mean/median/std, and write them as WS3D model files.

Returns `(mean_path, median_path, std_path)`.
"""
function AnalyseEnsemble3D(input_dir::AbstractString;
                           pattern::Regex = r"^chain.*\.rho$",
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

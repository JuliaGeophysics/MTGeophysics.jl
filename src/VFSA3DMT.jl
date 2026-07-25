# 3D VFSA MT inversion engine
# Author: @pankajkmishra
# Very Fast Simulated Annealing for MT, Gaussian-RBF parameterization, ModEM forward solver
# Includes ensemble statistics for multi-run analysis

using Random
using Dates

#---------- configuration ----------

"""
    VFSA3DMTConfig

Configuration for the 3D VFSA MT inversion. All VFSA hyper-parameters,
file-management switches, and external-solver settings live here. Every field is
a required keyword with no default — values are supplied by the run script (see
`examples/run_vfsa3D.jl`), so the engine carries no hidden defaults. The only
exceptions are fields added after that convention was set (`z_core_skin_depths`,
`z_core_cells`, `core_expand_cells`, the galvanic-distortion fields,
`fwd_ctrl`), which default to
sensible behavior and keep existing run scripts working unchanged.

The perturbable core is truncated in depth: deep layers exist for the forward
solver, not for the inversion. The core extends to `z_core_skin_depths` skin
depths `δ = 503·sqrt(ρ_bg·T_max)` derived from the observed data (median
off-diagonal apparent resistivity, longest period), or to the top `z_core_cells`
z layers when set. Controls are placed and perturbation applied only inside
that depth; `Inf` skin depths restores the full-column behavior.

With `distortion_mode = :on` each misfit evaluation solves, per site, a real
frequency-independent 2x2 matrix `C` minimizing the error-weighted misfit
`Zobs ~ C*Zpred` (variable projection, damped toward the identity by
`distortion_damping`), and anneals on the corrected rms. The per-site matrices
of the current best model are written to `distortion_best.txt`, sorted
by distortion severity, identifying which sites are galvanically distorted.

`out_root::String` is the base name/path of the run directory, which is also
the ModEM working folder. A timestamp is appended, so the actual run
directory is `<out_root>_<yyyymmdd_HHMMSS>` (e.g. `out_root="run"` gives
`run_20260626_143000`). A relative `out_root` is created alongside the starting
model file; an absolute path is used as-is. The iteration/trial logs and the
best model are written directly inside this run directory. One process runs one
chain; run several jobs with different seeds and out_roots for an ensemble.

Each iteration proposes `n_trials` models branching from the iteration-start
state, selects the lowest-RMS trial, and runs a single Metropolis test on that
winner (`n_trials = 1` is classic VFSA). One temperature schedule drives both
proposal width and acceptance, cooling from `T0` to
`cool_ratio*T0` at `max_iter`. Pick `T0` on the scale of the
typical uphill `dE_rms2` so the Metropolis test is selective from the start.

Control-point density can decrease with depth: sampling weight per cell is
`(depth + z1)^(-ctrl_depth_power)`, depth measured from the model top (0 =
uniform per cell). The weight is per cell, not per volume, so on meshes whose
dz grows with depth the mesh's own coarsening compounds with the exponent —
tune `ctrl_depth_power` for the grid at hand. Kernel widths grade from
`sigma_scale` (cells) at the top of the core to `sigma_scale_deep` at the
bottom, interpolated in log-depth over the core's own depth range; equal
values give uniform widths.

Frozen cells never receive control points, are never perturbed, and are never
bound-clamped. Topographic air (tagged as NaN in the loaded model) is always
frozen. Water is frozen additionally when a mask is active: from
`bathymetry_file` (written by [`write_bathymetry`](@ref)) if set, otherwise from
start-model cells with `log10 ρ < water_log10` (NaN disables); a model-derived
water mask is exported to `bathymetry.dat` in the run directory for reuse.
"""
Base.@kwdef mutable struct VFSA3DMTConfig
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
    T0::Float64
    cool_ratio::Float64
    target_rms::Float64
    seed::Int
    pad_tol::Float64
    # core_expand_cells: grow the perturbable lateral core by N cells per side
    # beyond the uniform-cell block, for sites that sit over the transition ring
    core_expand_cells::Int = 0
    padding_decay_length::Float64
    # z_core_skin_depths: perturbable core depth in data skin depths; Inf = full column
    z_core_skin_depths::Float64 = 1.0
    # z_core_cells: explicit core depth as the top N z layers; 0 = derive from data
    z_core_cells::Int = 0
    keep_models::Bool
    keep_dpred::Bool
    # model_save_every: 0 = behave per keep_models; >0 = retain trial models only
    # on iterations where iter % model_save_every == 0 (best model always kept)
    model_save_every::Int
    sigma_scale::Float64
    sigma_scale_deep::Float64
    trunc_sigmas::Float64
    # ctrl_depth_power: control-point sampling weight ∝ depth^(-p); 0 = uniform
    ctrl_depth_power::Float64
    # water_log10: freeze start-model cells below this log10 ρ; NaN = no mask
    water_log10::Float64
    # bathymetry_file: build the water mask from this file instead ("" = off)
    bathymetry_file::String
    # distortion_mode: :off = no distortion handling; :on = per-site 2x2 galvanic C
    distortion_mode::Symbol = :off
    # distortion_damping: pull of the per-site distortion matrix C toward identity,
    # lambda = a*tr(N)/2; 0 = free fit, larger = weaker correction, Inf = C fixed at
    # identity (correction off)
    distortion_damping::Float64 = 0.01
    # fwd_ctrl: ModEM forward-solver control file (F.dat); "" = omit the argument
    # and let the binary use its compiled-in defaults. Worth setting explicitly:
    # those defaults differ between ModEM builds (BICG in the maintained tree,
    # QMR in ModEM-GPU), and the CUDA build only accelerates BICG.
    fwd_ctrl::String = ""
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
    files = String[]
    for (root, _, fs) in walkdir(dir), f in fs
        occursin(pattern, f) && push!(files, joinpath(root, f))
    end
    sort(files)
end

#---------- core region helpers ----------

# expand grows the detected core by whole cells per side, clamped to the grid, so
# sites sitting over the transition ring just outside the uniform block still get
# perturbable cells under them
function core_ranges(m::WS3DModel; tol::Real=0.2, expand::Integer=0)
    ix = core_indices(m.cx; tol=tol)
    iy = core_indices(m.cy; tol=tol)
    e = max(0, Int(expand))
    e == 0 && return ix, iy
    return (max(1, first(ix) - e):min(m.nx, last(ix) + e),
            max(1, first(iy) - e):min(m.ny, last(iy) + e))
end

# bad z-core / decay settings fall back to defaults instead of aborting a queued job
function _sanitize_cfg!(cfg::VFSA3DMTConfig)
    if !(cfg.z_core_skin_depths > 0)
        @warn "z_core_skin_depths invalid ($(cfg.z_core_skin_depths)); using 1.0"
        cfg.z_core_skin_depths = 1.0
    end
    if cfg.z_core_cells < 0
        @warn "z_core_cells invalid ($(cfg.z_core_cells)); using 0 (derive from data)"
        cfg.z_core_cells = 0
    end
    if cfg.core_expand_cells < 0
        @warn "core_expand_cells invalid ($(cfg.core_expand_cells)); using 0"
        cfg.core_expand_cells = 0
    end
    if !(isfinite(cfg.padding_decay_length) && cfg.padding_decay_length > 0)
        @warn "padding_decay_length invalid ($(cfg.padding_decay_length)); using 10.0"
        cfg.padding_decay_length = 10.0
    end
    return cfg
end

# δ = 503·sqrt(ρ_bg·T_max), ρ_bg = median finite off-diagonal apparent resistivity
function _z_core_max_depth(dobs_abs::AbstractString, cfg::VFSA3DMTConfig)
    isfinite(cfg.z_core_skin_depths) || return Inf
    d = try
        load_data_modem(dobs_abs)
    catch e
        @warn "z-core: could not load observed data, using full model depth" error=e
        return Inf
    end
    rho_off = filter(x -> isfinite(x) && x > 0, vec(d.ρ[:, 2:3, :]))
    if isempty(rho_off) || isempty(d.T)
        @warn "z-core: no finite off-diagonal apparent resistivities, using full model depth"
        return Inf
    end
    return cfg.z_core_skin_depths * 503.0 * sqrt(median(rho_off) * maximum(d.T))
end

# cz carries origin[3], so cz[k] <= max_depth counts depth below the datum
function z_core_range(m::WS3DModel, max_depth::Real)
    kz_last = searchsortedlast(m.cz, max_depth)
    kz_last >= 1 || error("z-core depth $(max_depth) m is above the first layer center")
    return 1:kz_last
end

extract_core_array(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int}, kz::UnitRange{Int}) =
    @view m.A[ix, iy, kz]

function embed_core!(m::WS3DModel, A_core::AbstractArray{<:Real,3},
                     ix::UnitRange{Int}, iy::UnitRange{Int}, kz::UnitRange{Int})
    @assert size(A_core,1)==length(ix) && size(A_core,2)==length(iy) && size(A_core,3)==length(kz)
    @views m.A[ix, iy, kz] .= A_core
    return m
end

#---------- horizontal padding decay (x,y only, z untouched) ----------

function smooth_padding_decay_xy!(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int},
                                  background_res_log10::Float64, decay_length::Float64,
                                  protected::Union{Nothing,AbstractArray{Bool,3}}=nothing)
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
        protected !== nothing && protected[i, j, k] && continue
        ic = clamp(i, xi1, xi2)
        jc = clamp(j, yi1, yi2)
        dist = hypot((i - ic) * dx_median, (j - jc) * dy_median)
        weight = exp(-dist / max(L, eps()))
        boundary_val = m.A[ic, jc, k]
        # a frozen source must not bleed outward: seawater is finite, so without
        # this the ocean edge paints a conductive halo across dry padding
        usable = isfinite(boundary_val) &&
            !(protected !== nothing && protected[ic, jc, k])
        m.A[i, j, k] = usable ?
            boundary_val * weight + background_res_log10 * (1 - weight) :
            background_res_log10
    end
    return m
end

#---------- vertical padding decay (below the z core) ----------

# each core column continues its last perturbed value downward, blending toward
# background; runs before the xy decay so deep pad corners see the filled columns
#
# the carry-down halves every layer (50%, 25%, 12.5%, ...) all the way to the
# bottom of the model, counted in layers rather than metres: dz below the core
# grades so steeply that a physical e-fold collapses the whole blend into the
# first layer and leaves the rest of the column structureless
function smooth_padding_decay_z!(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int},
                                 kz::UnitRange{Int}, background_res_log10::Float64,
                                 protected::Union{Nothing,AbstractArray{Bool,3}}=nothing)
    klast = last(kz)
    klast == m.nz && return m
    @inbounds for k in (klast+1):m.nz
        weight = 0.5 ^ (k - klast)
        for j in iy, i in ix
            protected !== nothing && protected[i, j, k] && continue
            boundary_val = m.A[i, j, klast]
            usable = isfinite(boundary_val) &&
                !(protected !== nothing && protected[i, j, klast])
            m.A[i, j, k] = usable ?
                boundary_val * weight + background_res_log10 * (1 - weight) :
                background_res_log10
        end
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
    build_rbf_map(m, ix, iy, n_ctrl, rng; kz=1:m.nz, sigma_scale=2.0, trunc_sigmas=3.0,
                  sigma_scale_deep=sigma_scale, depth_power=0.0, exclude=nothing)

Randomly pick up to `n_ctrl` control voxels in the core and precompute
Gaussian-RBF weights for every core voxel. `kz` restricts the core in depth.
Kernel widths (in cells) grade from `sigma_scale` at the top of the core to
`sigma_scale_deep` at the bottom, interpolated per control in log-depth over the
core's own depth range, so the profile adapts to any survey without
absolute-depth settings; equal values give uniform widths. `depth_power` biases
control placement shallow and `exclude` marks cells that may not receive
controls.
"""
function build_rbf_map(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int}, n_ctrl::Int, rng::AbstractRNG;
                       kz::UnitRange{Int} = 1:m.nz,
                       sigma_scale::Float64 = 2.0, trunc_sigmas::Float64 = 3.0,
                       sigma_scale_deep::Float64 = sigma_scale,
                       depth_power::Float64 = 0.0,
                       exclude::Union{Nothing,AbstractArray{Bool,3}} = nothing)

    nx, ny, nz = length(ix), length(iy), length(kz)
    N = nx * ny * nz
    n_avail = exclude === nothing ? N : N - count(exclude)
    n_avail > 0 || error("build_rbf_map: every core cell is excluded")
    n_sel = min(n_ctrl, n_avail)
    n_sel < n_ctrl && @warn "build_rbf_map: only $n_sel of $n_ctrl controls placed (exclusion)"

    # depth-weighted, mask-aware placement: sampling weight ∝ (depth + z1)^(-depth_power)
    # with depth measured from the model-top edge so it stays monotonic across the
    # datum on air/topo grids; the weight is per cell, so a dz grading that coarsens
    # with depth compounds with depth_power. zero weight in excluded cells;
    # depth_power = 0 with no mask is uniform placement
    # (weighted sampling without replacement via exponential keys, Efraimidis–Spirakis)
    zc = view(m.cz, kz)
    z_top_edge = m.z[1]
    z1 = (zc[1] - z_top_edge) + eps()
    keys = Vector{Float64}(undef, N)
    @inbounds for id in 1:N
        k = Int(ceil(id / (nx*ny)))
        rem1 = id - (k-1)*nx*ny
        j = Int(ceil(rem1 / nx))
        i = rem1 - (j-1)*nx
        w = (exclude !== nothing && exclude[i, j, k]) ? 0.0 :
            (depth_power == 0.0 ? 1.0 : ((zc[k] - z_top_edge) + z1)^(-depth_power))
        keys[id] = w > 0 ? rand(rng)^(1.0 / w) : -Inf
    end
    idxs = collect(partialsortperm(keys, 1:n_sel; rev=true))
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

    z_top = zc[1] - z_top_edge
    z_bot = zc[nz] - z_top_edge
    denom = log(z_bot + z1) - log(z_top + z1)
    sig = Vector{Float64}(undef, n_sel)
    @inbounds for q in 1:n_sel
        t = denom > 0 ? clamp((log((zc[ck[q]] - z_top_edge) + z1) - log(z_top + z1)) / denom, 0.0, 1.0) : 0.0
        sig[q] = sigma_scale + t * (sigma_scale_deep - sigma_scale)
    end

    cell_ids = [Int[] for _ in 1:N]
    cell_wts = [Float64[] for _ in 1:N]
    @inbounds for q in 1:n_sel
        sq = sig[q]
        rq = ceil(Int, trunc_sigmas * sq)
        for kk in max(1, ck[q]-rq):min(nz, ck[q]+rq),
            jj in max(1, cj[q]-rq):min(ny, cj[q]+rq),
            ii in max(1, ci[q]-rq):min(nx, ci[q]+rq)
            # index-space metric: isotropic in cells, anisotropic in meters on graded dz
            r2 = ((ii - ci[q])^2 + (jj - cj[q])^2 + (kk - ck[q])^2) / sq^2
            if r2 <= trunc_sigmas^2
                id = ii + (jj-1)*nx + (kk-1)*nx*ny
                push!(cell_ids[id], q)
                push!(cell_wts[id], exp(-0.5 * r2))
            end
        end
    end

    ptr  = Vector{Int}(undef, N + 1)
    nbrs = Int[];       sizehint!(nbrs, 8 * N)
    wts  = Float64[];   sizehint!(wts,  8 * N)

    ptr[1] = 1
    row = 1
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        local_ids = cell_ids[row]
        local_wts = cell_wts[row]

        if isempty(local_ids)
            best_q = 1
            best_r2 = typemax(Float64)
            for q in 1:n_sel
                r2 = ((i - ci[q])^2 + (j - cj[q])^2 + (k - ck[q])^2) / sig[q]^2
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
    # absolute: the command runs with run_dir as cwd
    fwd_ctrl_abs = isempty(cfg.fwd_ctrl) ? "" : abspath(cfg.fwd_ctrl)
    cd(run_dir) do
        # -F takes the control file in slot 5, behind wFile_EMsoln; "n" is a
        # placeholder that leaves it unwritten (ModEM gates the write on
        # len_trim(wFile_EMsoln) > 1), so no EM solution is dumped per eval
        cmd = isempty(fwd_ctrl_abs) ?
            `$(cfg.mpirun_cmd) -n $(cfg.nprocs) $(cfg.modem_exe) -F $(model_filename) $(dobs_filename) $(dpred_filename)` :
            `$(cfg.mpirun_cmd) -n $(cfg.nprocs) $(cfg.modem_exe) -F $(model_filename) $(dobs_filename) $(dpred_filename) n $(fwd_ctrl_abs)`
        try
            run(pipeline(cmd, stdout="0runlog.dat"))
        catch e
            @warn("[forward] ModEM failed", error=e)
            ok = false
        end
    end
    if !ok || !isfile(dpred_abs)
        _cleanup_modem_artifacts!(run_dir)
        return Inf, nothing
    end
    if cfg.distortion_mode == :on
        fit = chi2_and_rms_distorted(dobs_abs, dpred_abs; damping=cfg.distortion_damping)
        _cleanup_modem_artifacts!(run_dir)
        return fit.rms, fit
    end
    s = chi2_and_rms(dobs_abs, dpred_abs; use_impedance=true, use_tipper=true, components=String[])
    _cleanup_modem_artifacts!(run_dir)
    return s.rms, nothing
end

#---------- logging ----------

function _write_cfg_header_3d(io::IO; cfg::VFSA3DMTConfig, ak::Float64,
                              start_model::String, dobs_filename::String)
    names = fieldnames(VFSA3DMTConfig)
    w = maximum(length(String(f)) for f in names)
    println(io, "# ", rpad("start_model", w), " = ", start_model)
    println(io, "# ", rpad("dobs", w), " = ", dobs_filename)
    println(io, "# ", rpad("ak", w), " = ", ak, "   (derived: cool_ratio, max_iter)")
    for f in names
        println(io, "# ", rpad(String(f), w), " = ", getfield(cfg, f))
    end
end

function _write_trials_header_3d(path::AbstractString; timestamp::String, cfg::VFSA3DMTConfig,
                                 ak::Float64, start_model::String, dobs_filename::String)
    open(path, "w") do io
        println(io, "# VFSA 3D MT detailed trials — ", timestamp)
        _write_cfg_header_3d(io; cfg=cfg, ak=ak,
                             start_model=start_model, dobs_filename=dobs_filename)
        println(io, repeat("-", 102))
        @printf(io, "%8s %8s %12s %11s %12s %10s %10s %5s %11s %s\n",
                "Iter","Trial","Temp",
                "RMS","dE_rms2","Pacc","Uacc","Acc","RMSBest","Model")
        println(io, repeat("-", 102))
    end
end

function _append_trial_row_3d(path::AbstractString; iter::Int, trial::Int,
                           T::Float64,
                           rms_prop::Float64,
                           dE::Float64, p_acc::Float64, u_acc::Float64,
                           rms_best::Float64,
                           acc::Int, model_rel::String)
    # one row per trial: the proposal's own misfit, energy dE, metropolis
    # prob/draw, accept flag, running best, and the proposal model filename
    # non-winning greedy trials carry no metropolis draw, print an em dash
    pacc_str = isnan(p_acc) ? lpad("—", 10) : @sprintf("%10.5f", p_acc)
    uacc_str = isnan(u_acc) ? lpad("—", 10) : @sprintf("%10.5f", u_acc)
    open(path, "a") do io
        @printf(io, "%8d %8d %12.4g %11.5f %12.5f %s %s %5d %11.5f %s\n",
                iter, trial, T, rms_prop, dE, pacc_str, uacc_str, acc, rms_best, model_rel)
    end
end

function _write_iter_header_3d(path::AbstractString; timestamp::String, cfg::VFSA3DMTConfig,
                               ak::Float64, start_model::String, dobs_filename::String)
    open(path, "w") do io
        println(io, "# VFSA 3D MT iteration best — ", timestamp)
        _write_cfg_header_3d(io; cfg=cfg, ak=ak,
                             start_model=start_model, dobs_filename=dobs_filename)
        println(io, repeat("-", 88))
        # current accepted state after the iteration, Nacc = accepted trials this iter
        @printf(io, "%8s %12s %11s %11s %6s %s\n",
                "Iter","Temp","RMS","RMSBest","Nacc","Model")
        println(io, repeat("-", 88))
    end
end

function _append_iter_row_3d(path::AbstractString; iter::Int,
                          T::Float64,
                          rms_curr::Float64,
                          rms_best::Float64,
                          nacc::Int, model_rel::String)
    open(path, "a") do io
        @printf(io, "%8d %12.4g %11.5f %11.5f %6d %s\n",
                iter, T, rms_curr, rms_best, nacc, model_rel)
    end
end

_T_schedule(k::Int; T0::Float64, ak::Float64) = T0 * exp(-sqrt(ak) * (k - 1))

# decay rate so T reaches cool_ratio*T0 at max_iter (shared with the 2D driver)
_resolve_ak(cool_ratio::Float64, max_iter::Int) = (log(1.0 / cool_ratio) / max(max_iter - 1.0, 1.0))^2

#---------- vfsa loop ----------

function _run_vfsa3d(start_model_path::String,
                     dobs_filename::String, timestamp::String,
                     run_root::String; cfg::VFSA3DMTConfig)
    rng = MersenneTwister(cfg.seed)

    trials_log = joinpath(run_root, "0vfsa3DMT_detailed.log")
    iter_log   = joinpath(run_root, "0vfsa3DMT.log")

    m = load_ws3d_model(start_model_path)
    _, _, _, _, _, _, origin, rotation = read_ws3d_model(start_model_path, true)

    ix, iy = core_ranges(m; tol=cfg.pad_tol, expand=cfg.core_expand_cells)
    kz = cfg.z_core_cells > 0 ? (1:min(cfg.z_core_cells, m.nz)) :
         z_core_range(m, _z_core_max_depth(joinpath(run_root, dobs_filename), cfg))
    last(kz) < m.nz && @info @sprintf("z-core: perturbing %d of %d layers (bottom %.0f m)",
                                      length(kz), m.nz, m.cz[last(kz)])
    Acore = extract_core_array(m, ix, iy, kz)
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

    #---------- frozen cells: no controls, no perturbation, no clamp ----------
    # topography and bathymetry are tracked as separate full-grid masks so the
    # padding ring is protected too, and both are exported for reuse/plotting
    air_full = air_mask_from_model(m)
    water_full = falses(size(m.A))
    bathy = nothing
    if !isempty(cfg.bathymetry_file)
        bathy = read_bathymetry(cfg.bathymetry_file)
        water_full .= water_mask_from_bathymetry(m, bathy)
    elseif !isnan(cfg.water_log10)
        water_full .= water_mask_from_model(m; water_log10=cfg.water_log10)
        any(water_full) && (bathy = extract_bathymetry(m; water_log10=cfg.water_log10))
    end
    protected_full = air_full .| water_full

    if any(air_full)
        write_topography(joinpath(run_root, "topography.dat"), extract_topography(m);
                         origin=m.origin)
    end
    if bathy !== nothing && !isempty(bathy.depth)
        write_bathymetry(joinpath(run_root, "bathymetry.dat"), bathy;
                         water_log10=cfg.water_log10, origin=m.origin)
    end

    mask = BitArray(protected_full[ix, iy, kz])
    n_air = count(air_full[ix, iy, kz])
    n_water = count(water_full[ix, iy, kz])
    n_frozen = count(mask)

    # build the RBF map
    rbfmap = build_rbf_map(m, ix, iy, cfg.n_ctrl, rng;
                           kz=kz,
                           sigma_scale=cfg.sigma_scale, trunc_sigmas=cfg.trunc_sigmas,
                           sigma_scale_deep=cfg.sigma_scale_deep,
                           depth_power=cfg.ctrl_depth_power,
                           exclude=(n_frozen > 0 ? mask : nothing))
    M = length(rbfmap.ci)
    if n_frozen > 0 || cfg.ctrl_depth_power != 0.0
        z_ctrl = [m.cz[k] for k in rbfmap.ck]
        @info @sprintf("mask: %d air + %d water cells frozen; %d controls, median depth %.0f m, %d above 6 km",
                       n_air, n_water, M, median(z_ctrl), count(<(6000.0), z_ctrl))
    end

    v0_ctrl = Vector{Float64}(undef, M)
    @inbounds for t in 1:M
        v0_ctrl[t] = v0_core_log10[rbfmap.ci[t], rbfmap.cj[t], rbfmap.ck[t]]
    end

    delta_params_current = zeros(Float64, M)
    delta_values_buffer = zeros(Float64, rbfmap.nx, rbfmap.ny, rbfmap.nz)
    nsel_default = max(1, round(Int, cfg.frac_update_controls * M))

    # initial forward on the start model
    model0_filename = @sprintf("model_%03d_%02d.rho", 0, 0)
    dpred0_filename = @sprintf("dpred_%03d_%02d.dat", 0, 0)

    embed_core!(m, v0_core_log10, ix, iy, kz)
    smooth_padding_decay_z!(m, ix, iy, kz, background_log10, protected_full)
    smooth_padding_decay_xy!(m, ix, iy, background_log10, cfg.padding_decay_length,
                             protected_full)

    dp0, fit0 = forward_and_misfit!(m; run_dir=run_root, model_filename=model0_filename,
                                    dobs_filename=dobs_filename, dpred_filename=dpred0_filename,
                                    origin=origin, rotation=rotation, cfg=cfg)
    rms_current = dp0
    best_rms    = rms_current
    T0 = cfg.T0
    ak = _resolve_ak(cfg.cool_ratio, cfg.max_iter)
    _write_trials_header_3d(trials_log; timestamp=timestamp, cfg=cfg, ak=ak,
                            start_model=start_model_path, dobs_filename=dobs_filename)
    _write_iter_header_3d(iter_log; timestamp=timestamp, cfg=cfg, ak=ak,
                          start_model=start_model_path, dobs_filename=dobs_filename)
    best_model_abs = joinpath(run_root, "best_model.rho")
    cp(joinpath(run_root, model0_filename), best_model_abs; force=true)
    # per-site distortion of the current best model, refreshed on every new best
    distortion_abs = joinpath(run_root, "distortion_best.txt")
    fit0 !== nothing && write_distortion_file(distortion_abs, fit0;
                                              iter=0, rms=dp0, damping=cfg.distortion_damping)
    # filename of the current accepted model, carried forward so the iteration
    # log always names a model even on iterations with no accepted trial
    current_model_rel = model0_filename

    lo, hi = cfg.log_bounds

    for k in 1:cfg.max_iter
        # one schedule drives both proposal width and acceptance
        T = _T_schedule(k; T0=T0, ak=ak)

        #---------- generate and test the trials ----------
        # all trials branch from the iteration-start state, the lowest-rms trial
        # gets one deferred metropolis test
        n_accepted_iter = 0
        trial_cache = NamedTuple[]
        rms_iter_start = rms_current
        best_trial_rms = Inf
        best_trial_delta_params = similar(delta_params_current)
        best_trial_model_rel = ""
        best_trial_idx = 0
        best_trial_fit = nothing

        for t in 1:cfg.n_trials
            # propose from the iteration-start state, the chain only advances
            # after the trial loop
            delta_params_trial = copy(delta_params_current)
            propose_controls!(delta_params_trial, T, lo, hi, v0_ctrl, nsel_default, rng;
                              step_scale=cfg.step_scale)

            apply_rbf_map!(delta_values_buffer, rbfmap, delta_params_trial)
            v_trial = clamp.(v0_core_log10 .+ delta_values_buffer, lo, hi)
            # frozen water cells: undo kernel bleed-in and bound clamping
            n_frozen > 0 && (v_trial[mask] .= v0_core_log10[mask])

            embed_core!(m, v_trial, ix, iy, kz)
            smooth_padding_decay_z!(m, ix, iy, kz, background_log10, protected_full)
            smooth_padding_decay_xy!(m, ix, iy, background_log10, cfg.padding_decay_length,
                                     protected_full)

            model_filename = @sprintf("model_%03d_%02d.rho", k, t)
            dpred_filename = @sprintf("dpred_%03d_%02d.dat", k, t)
            dp, fit_t = forward_and_misfit!(m; run_dir=run_root, model_filename=model_filename,
                                            dobs_filename=dobs_filename, dpred_filename=dpred_filename,
                                            origin=origin, rotation=rotation, cfg=cfg)

            dE = (dp^2 - rms_iter_start^2) / max(rms_iter_start^2, eps())
            push!(trial_cache, (
                iter=k, trial=t, T=T,
                rms_prop=dp,
                dE=dE, p_acc=NaN, u_acc=NaN, accepted=false,
                rms_best=best_rms,
                model_rel=model_filename
            ))
            if dp < best_trial_rms
                best_trial_rms = dp
                best_trial_delta_params .= delta_params_trial
                best_trial_model_rel = model_filename
                best_trial_idx = t
                best_trial_fit = fit_t
            end
        end

        #---------- single metropolis test on the best trial ----------
        dE_best = (best_trial_rms^2 - rms_iter_start^2) / max(rms_iter_start^2, eps())
        u_acc_best = rand(rng)
        # Pacc = 1 downhill, exp(-dE/T) uphill, kept finite for the log
        p_acc_best = dE_best <= 0 ? 1.0 : exp(-dE_best / max(T, 1e-12))
        accept_best = isfinite(best_trial_rms) && (u_acc_best < p_acc_best)
        if accept_best
            delta_params_current .= best_trial_delta_params
            rms_current = best_trial_rms
            n_accepted_iter = 1
            current_model_rel = best_trial_model_rel
            # record the global best on improvement, never feed it back
            if rms_current < best_rms
                best_rms = rms_current
                cp(joinpath(run_root, best_trial_model_rel), best_model_abs; force=true)
                best_trial_fit !== nothing && write_distortion_file(distortion_abs, best_trial_fit;
                                                                    iter=k, rms=best_rms,
                                                                    damping=cfg.distortion_damping)
            end
        end
        # backfill the winning trial's metropolis draw into its log row
        if best_trial_idx > 0
            tr = trial_cache[best_trial_idx]
            trial_cache[best_trial_idx] = merge(tr, (p_acc=p_acc_best, u_acc=u_acc_best, accepted=accept_best))
        end
        # trial files can only be removed after the best-model copy above;
        # in-loop deletion would lose the winning trial's model file
        if cfg.model_save_every <= 0
            for t in 1:cfg.n_trials
                if !cfg.keep_models
                    model_t = joinpath(run_root, @sprintf("model_%03d_%02d.rho", k, t))
                    isfile(model_t) && rm(model_t; force=true)
                end
                if !cfg.keep_dpred
                    dpred_t = joinpath(run_root, @sprintf("dpred_%03d_%02d.dat", k, t))
                    isfile(dpred_t) && rm(dpred_t; force=true)
                end
            end
        end

        #---------- checkpointed pruning (model_save_every > 0) ----------
        # keep all trial files only on checkpoint iters; on a checkpoint keep
        # just the winning trial (accepted or not); otherwise drop the whole
        # iteration. the best model lives in best_model_abs and is never touched here.
        if cfg.model_save_every > 0
            is_checkpoint = (k % cfg.model_save_every == 0)
            keep_trial = best_trial_idx > 0 ? best_trial_idx : cfg.n_trials
            for t in 1:cfg.n_trials
                if is_checkpoint && t == keep_trial
                    continue
                end
                model_t = joinpath(run_root, @sprintf("model_%03d_%02d.rho", k, t))
                isfile(model_t) && rm(model_t; force=true)
                dpred_t = joinpath(run_root, @sprintf("dpred_%03d_%02d.dat", k, t))
                isfile(dpred_t) && rm(dpred_t; force=true)
            end
        end

        # write one row per trial
        for tr in trial_cache
            _append_trial_row_3d(trials_log; iter=tr.iter, trial=tr.trial,
                              T=tr.T,
                              rms_prop=tr.rms_prop,
                              dE=tr.dE, p_acc=tr.p_acc, u_acc=tr.u_acc,
                              rms_best=tr.rms_best,
                              acc=(tr.accepted ? 1 : 0), model_rel=tr.model_rel)
        end

        # iteration summary row: current accepted state, named by the model that
        # represents it (carried forward when no trial was accepted this iter)
        _append_iter_row_3d(iter_log; iter=k,
                         T=T,
                         rms_curr=rms_current,
                         rms_best=best_rms,
                         nacc=n_accepted_iter,
                         model_rel=current_model_rel)

        # stop once best rms hits the target
        if best_rms <= cfg.target_rms
            @info @sprintf("target RMS %.3f reached at iter %d (best RMS %.5f); stopping.",
                           cfg.target_rms, k, best_rms)
            break
        end
    end
    return best_model_abs, iter_log
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
    cfg.ctrl_depth_power >= 0 ||
        error("cfg.ctrl_depth_power must be >= 0, got $(cfg.ctrl_depth_power)")
    cfg.sigma_scale_deep > 0 ||
        error("cfg.sigma_scale_deep must be > 0, got $(cfg.sigma_scale_deep)")
    _sanitize_cfg!(cfg)
    isempty(cfg.bathymetry_file) || isfile(cfg.bathymetry_file) ||
        error("cfg.bathymetry_file not found: $(cfg.bathymetry_file)")
    isempty(cfg.fwd_ctrl) || isfile(cfg.fwd_ctrl) ||
        error("cfg.fwd_ctrl not found: $(cfg.fwd_ctrl)")
    cfg.distortion_mode in (:off, :on) ||
        error("cfg.distortion_mode must be :off or :on, got $(cfg.distortion_mode)")
    cfg.distortion_mode == :off || cfg.distortion_damping >= 0 ||
        error("cfg.distortion_damping must be >= 0, got $(cfg.distortion_damping)")
    start_model_abs = abspath(start_model_path)
    model_dir = dirname(start_model_abs)

    run_root_base = isabspath(cfg.out_root) ? cfg.out_root : joinpath(model_dir, cfg.out_root)
    t_now = now()
    run_root = string(run_root_base, "_", Dates.format(t_now, "yyyymmdd_HHMMSS"))
    isdir(run_root) || mkpath(run_root)

    dobs_filename = basename(dobs_path)
    dobs_abs_target = abspath(joinpath(run_root, dobs_filename))
    cp(abspath(dobs_path), dobs_abs_target; force=true)

    isfile(PROGRAM_FILE) &&
        cp(abspath(PROGRAM_FILE), joinpath(run_root, basename(PROGRAM_FILE)); force=true)

    timestamp = Dates.format(t_now, "yyyy-mm-dd HH:MM:SS")

    return _run_vfsa3d(start_model_abs, dobs_filename, timestamp, run_root; cfg=cfg)
end

#---------- ensemble statistics ----------

function _intersect_ranges(r1::UnitRange{Int}, r2::UnitRange{Int})
    a = max(first(r1), first(r2)); b = min(last(r1), last(r2))
    b < a && error("Empty intersection across runs.")
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
        error("Grid mismatch across runs within the overlapping index ranges.")
    end
    if length(mref.origin) >= 3 && length(m.origin) >= 3
        dor = maximum(abs.(mref.origin[1:3] .- m.origin[1:3]))
        dor <= 1e-6 || error("Origin mismatch across runs.")
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

Load all matching model files from `input_dir` (searched recursively), validate
grid consistency, compute ensemble mean/median/std, and write them as WS3D model
files. Each VFSA3DMT job writes `best_model.rho` into its own run directory, so
point `input_dir` at the parent directory holding the run directories.

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

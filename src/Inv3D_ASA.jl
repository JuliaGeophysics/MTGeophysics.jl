# 3D ASA MT inversion engine
# Author: @pankajkmishra
# Adaptive simulated annealing (ASA) for 3D MT on the Gaussian-RBF control points of VFSA3DMT (same core, padding
# blend, masks and external forward solver), written from the ASA equations rather than from any implementation.
#
# Generating function. A control q in [A, B] = log_bounds moves by y*(B - A), y drawn from the ASA density
# g(y) = 1 / (2 (|y| + T_q) ln(1 + 1/T_q)) on [-1, 1], i.e. y = sgn(u - 1/2) T_q ((1 + 1/T_q)^|2u-1| - 1), redrawn
# until the control stays inside the bounds. Generating all controls at once with this fat-tailed density puts
# hundreds of large jumps in every state, so states are generated sequentially in blocks, the alternative ASA offers:
# D = 16 controls per state, taken in turn from a shuffled sweep so every control is generated equally often.
#
# Schedule. Each control carries its own annealing index k_q, the number of times it has been generated, and
# temperature T_q(k_q) = T0 exp(-c k_q^(1/D)), c = m exp(-n/D), which gives T_f = T0 exp(-m) at k_f = exp(n). The
# sufficiency argument sum_k prod_i g(y_i) ~ sum_k 1/k = infinity takes the product over the D parameters generated
# together, so D is the block size and the schedule is ASA proper, no quenching. As in the paper T0 = 1 spans the
# whole range, m = ln(1e5) (temperature ratio scale 1e-5) and n = ln(k_f) with k_f the number of generations each
# control receives in max_iter states, so every control ends its run at T0 exp(-m). The acceptance test is
# Boltzmann, exp(-dL/T_acc) on the cost L = rms^2, with the acceptance temperature on the same schedule but indexed by
# the number of accepted states, starting at the cost of the initial trial, T_acc0 = L(start model).
#
# Reannealing, every 100 accepted states. ASA rescales each parameter temperature by s_max/s_q, the sensitivity of
# the most sensitive parameter over its own, so that insensitive parameters are searched more widely. The sensitivity
# used here is the diagonal of the Gauss-Newton Hessian in control space, H = W' J' Cd^-1 J W (W the RBF map), the
# curvature whose inverse square root ASA reports as the parameter standard deviation. It is estimated with two
# Rademacher probes per reanneal averaged into a running estimate: diag(H) = E[(W' J' Cd^-1/2 z)^2] for random signs z
# on the data, each probe one J' product of the forward solver, about 1.5 forward runs. Applied repeatedly, the
# paper's rescaling compounds and drives every control but the stiffest back to T0 when D is large; here it is applied
# to the scheduled temperature instead, T_q = min(T0, w_q T(k_q)), w_q = H_ref / H_q clipped to [1, 100], H_ref the
# 99th percentile. The cost temperature is reannealed as in the paper: T_acc0 is lowered to the scale of the last and
# best costs, max(L_last, L_best, |L_last - L_best|), T_acc is kept but no higher than the new T_acc0, and its index
# is reset so the temperature is continuous.
#
# References:
# Ingber, L. (1989). Very fast simulated re-annealing. Mathl. Comput. Modelling, 12(8), 967-973.
# Ingber, L. (1993). Simulated annealing: practice versus theory. Mathl. Comput. Modelling, 18(11), 29-57.
# Ingber, L. (1996). Adaptive simulated annealing (ASA): lessons learned. Control and Cybernetics, 25(1), 33-54.
# Hutchinson, M. F. (1990). A stochastic estimator of the trace of the influence matrix for Laplacian smoothing
# splines. Commun. Stat. Simul. Comput., 19(2), 433-450.

using Random
using Dates
using Printf
using Statistics

const _ASA_BLOCK = 16
const _ASA_T0 = 1.0
const _ASA_TEMP_RATIO_SCALE = 1e-5
const _ASA_REANNEAL_EVERY = 100
const _ASA_PROBES = 2
const _ASA_REHEAT_MAX = 100.0
const _ASA_REF_QUANTILE = 0.99

"""
    ASA3DMTConfig

Configuration for the 3D adaptive simulated annealing (ASA) MT inversion. The
model parameterization, masks, padding and solver fields mean exactly what they
mean in [`VFSA3DMTConfig`](@ref), so a VFSA run script converts by deleting the
annealing fields. There are no temperatures, cooling rates, step widths, update
fractions or trial counts to set: the schedules follow the ASA paper (T0 = 1,
final ratio 1e-5 at the end of the budget, acceptance temperature starting at
the initial cost), every control's generating temperature is rescaled from the
data sensitivities, and the remaining constants are fixed in the source (see the
header of `src/Inv3D_ASA.jl`).

- `nprocs`, `mpirun_cmd`, `modem_exe`, `fwd_ctrl`: forward solver launch, as in VFSA3DMT
- `out_root`: run directory base, a timestamp is appended (relative = next to the start model)
- `max_iter`: generated states, one forward run each; it also sets the annealing timescale
- `target_rms`: stop once the best rms reaches it
- `seed`: random seed (control placement, proposals, probes)
- `n_ctrl`, `log_bounds`: RBF control points and the log10 ρ box
- `pad_tol`, `core_expand_cells`, `padding_decay_length`: lateral core and padding blend
- `z_core_skin_depths`, `z_core_cells`: perturbable core depth
- `sigma_scale`, `sigma_scale_deep`, `trunc_sigmas`, `ctrl_depth_power`: RBF widths and placement
- `water_log10`, `bathymetry_file`: frozen water cells
- `keep_models`, `keep_dpred`, `model_save_every`: file retention, as in VFSA3DMT
- `distortion_mode`, `distortion_damping`: per-site galvanic distortion, as in VFSA3DMT
"""
Base.@kwdef mutable struct ASA3DMTConfig
    nprocs::Int
    mpirun_cmd::String
    modem_exe::String
    fwd_ctrl::String = ""
    out_root::String
    max_iter::Int
    target_rms::Float64
    seed::Int
    n_ctrl::Int
    log_bounds::Tuple{Float64,Float64}
    pad_tol::Float64
    core_expand_cells::Int = 0
    padding_decay_length::Float64
    z_core_skin_depths::Float64 = 1.0
    z_core_cells::Int = 0
    sigma_scale::Float64
    sigma_scale_deep::Float64
    trunc_sigmas::Float64
    ctrl_depth_power::Float64
    water_log10::Float64
    bathymetry_file::String
    keep_models::Bool
    keep_dpred::Bool
    model_save_every::Int
    distortion_mode::Symbol = :off
    distortion_damping::Float64 = 0.01
end

function _asa_check_cfg(cfg::ASA3DMTConfig)
    cfg.max_iter >= 1 || error("cfg.max_iter must be >= 1, got $(cfg.max_iter)")
    cfg.n_ctrl >= 1 || error("cfg.n_ctrl must be >= 1, got $(cfg.n_ctrl)")
    lo, hi = cfg.log_bounds
    hi > lo || error("cfg.log_bounds must be increasing, got $(cfg.log_bounds)")
    cfg.ctrl_depth_power >= 0 || error("cfg.ctrl_depth_power must be >= 0, got $(cfg.ctrl_depth_power)")
    cfg.sigma_scale > 0 && cfg.sigma_scale_deep > 0 ||
        error("cfg.sigma_scale and cfg.sigma_scale_deep must be > 0")
    cfg.z_core_skin_depths > 0 || error("cfg.z_core_skin_depths must be > 0, got $(cfg.z_core_skin_depths)")
    cfg.z_core_cells >= 0 || error("cfg.z_core_cells must be >= 0, got $(cfg.z_core_cells)")
    cfg.core_expand_cells >= 0 || error("cfg.core_expand_cells must be >= 0, got $(cfg.core_expand_cells)")
    isfinite(cfg.padding_decay_length) && cfg.padding_decay_length > 0 ||
        error("cfg.padding_decay_length must be > 0, got $(cfg.padding_decay_length)")
    isempty(cfg.bathymetry_file) || isfile(cfg.bathymetry_file) ||
        error("cfg.bathymetry_file not found: $(cfg.bathymetry_file)")
    isempty(cfg.fwd_ctrl) || isfile(cfg.fwd_ctrl) || error("cfg.fwd_ctrl not found: $(cfg.fwd_ctrl)")
    cfg.distortion_mode in (:off, :on) ||
        error("cfg.distortion_mode must be :off or :on, got $(cfg.distortion_mode)")
    cfg.distortion_mode == :off || cfg.distortion_damping >= 0 ||
        error("cfg.distortion_damping must be >= 0, got $(cfg.distortion_damping)")
    return cfg
end

_asa_temperature(T0::Float64, c::Float64, k::Real, D::Int) = T0 * exp(-c * float(k)^(1 / D))

_asa_index(T0::Float64, T::Float64, c::Float64, D::Int) = T >= T0 ? 0.0 : (log(T0 / T) / c)^D

function _asa_generate!(delta::Vector{Float64}, idxs::AbstractVector{Int}, Tq::Vector{Float64},
                        v0_ctrl::Vector{Float64}, lo::Float64, hi::Float64, rng::AbstractRNG)
    width = hi - lo
    @inbounds for q in idxs
        base = v0_ctrl[q] + delta[q]
        cand = base
        for _ in 1:100
            c = base + vfsa_y(rand(rng), Tq[q]) * width
            if lo <= c <= hi
                cand = c
                break
            end
        end
        delta[q] = cand - v0_ctrl[q]
    end
    return delta
end

function _asa_rbf_transpose(rbfmap::RBFMap, g::AbstractArray{<:Real,3})
    out = zeros(Float64, length(rbfmap.ci))
    row = 1
    @inbounds for k in 1:rbfmap.nz, j in 1:rbfmap.ny, i in 1:rbfmap.nx
        gv = g[i, j, k]
        for p in rbfmap.ptr[row]:(rbfmap.ptr[row+1] - 1)
            out[rbfmap.nbrs[p]] += rbfmap.wts[p] * gv
        end
        row += 1
    end
    return out
end

function _asa_write_probe(dobs_abs::AbstractString, probe_abs::AbstractString, rng::AbstractRNG)
    open(probe_abs, "w") do io
        for line in eachline(dobs_abs)
            t = split(line)
            if length(t) >= 11 && !startswith(t[1], "#") && !startswith(t[1], ">") &&
               tryparse(Float64, t[1]) !== nothing
                err = something(tryparse(Float64, t[11]), NaN)
                w = (isfinite(err) && err > 0) ? 1.0 / err : 0.0
                re = rand(rng, Bool) ? w : -w
                im = rand(rng, Bool) ? w : -w
                println(io, join(t[1:8], " "), @sprintf(" %15.6E %15.6E %15.6E", re, im, err))
            else
                println(io, line)
            end
        end
    end
    return probe_abs
end

function _asa_read_record(io::IO)
    n = read(io, Int32)
    b = read(io, n)
    read(io, Int32) == n || error("corrupt record in J' output")
    return b
end

function _asa_read_jt(path::AbstractString, nx::Int, ny::Int, nz::Int)
    tot = zeros(Float64, nx, ny, nz)
    open(path) do io
        r = _asa_read_record(io)
        while length(r) != 4
            r = _asa_read_record(io)
        end
        ntx = reinterpret(Int32, r)[1]
        for _ in 1:ntx
            _asa_read_record(io)
            dims = Int.(reinterpret(Int32, _asa_read_record(io)))
            dims == [nx, ny, nz] || error("J' output grid $(dims) does not match the model ($nx, $ny, $nz)")
            for _ in 1:4
                _asa_read_record(io)
            end
            tot .+= reshape(reinterpret(Float64, _asa_read_record(io)), nx, ny, nz)
        end
    end
    return tot
end

function _asa_gn_diagonal(m::WS3DModel, rbfmap::RBFMap, ix::UnitRange{Int}, iy::UnitRange{Int},
                          kz::UnitRange{Int}, mask::BitArray{3}; run_dir::AbstractString,
                          dobs_filename::String, origin::Vector{Float64}, rotation::Float64,
                          cfg::ASA3DMTConfig, rng::AbstractRNG)
    model_filename = "sens_model.rho"
    write_ws3d_model(joinpath(run_dir, model_filename), m.dx, m.dy, m.dz, m.A, origin;
                     rotation=rotation, type_str="LOGE")
    fwd_ctrl_abs = isempty(cfg.fwd_ctrl) ? "" : abspath(cfg.fwd_ctrl)
    h = zeros(Float64, length(rbfmap.ci))
    n_ok = 0
    for p in 1:_ASA_PROBES
        probe_filename = "sens_probe.dat"
        out_filename = "sens_jt.bin"
        _asa_write_probe(joinpath(run_dir, dobs_filename), joinpath(run_dir, probe_filename), rng)
        ok = true
        cd(run_dir) do
            cmd = isempty(fwd_ctrl_abs) ?
                `$(cfg.mpirun_cmd) -n $(cfg.nprocs) $(cfg.modem_exe) -x $(model_filename) $(probe_filename) $(out_filename)` :
                `$(cfg.mpirun_cmd) -n $(cfg.nprocs) $(cfg.modem_exe) -x $(model_filename) $(probe_filename) $(out_filename) $(fwd_ctrl_abs)`
            try
                run(pipeline(cmd, stdout="0runlog_sens.dat"))
            catch e
                @warn("[sensitivity] J' product failed", error=e)
                ok = false
            end
        end
        out_abs = joinpath(run_dir, out_filename)
        if ok && isfile(out_abs)
            g = try
                _asa_read_jt(out_abs, m.nx, m.ny, m.nz)
            catch e
                @warn("[sensitivity] could not read J' output", error=e)
                nothing
            end
            if g !== nothing
                gc = g[ix, iy, kz]
                gc[mask] .= 0.0
                h .+= _asa_rbf_transpose(rbfmap, gc) .^ 2
                n_ok += 1
            end
        end
        rm(joinpath(run_dir, probe_filename); force=true)
        rm(out_abs; force=true)
        _cleanup_modem_artifacts!(run_dir)
    end
    rm(joinpath(run_dir, model_filename); force=true)
    return n_ok > 0 ? h ./ n_ok : nothing, n_ok
end

function _asa_reheat(h::Vector{Float64})
    pos = filter(x -> isfinite(x) && x > 0, h)
    isempty(pos) && return ones(length(h)), NaN
    href = quantile(pos, _ASA_REF_QUANTILE)
    w = [(isfinite(x) && x > 0) ? clamp(href / x, 1.0, _ASA_REHEAT_MAX) : _ASA_REHEAT_MAX for x in h]
    return w, href
end

function _asa_write_header(path::AbstractString, title::AbstractString; timestamp::String,
                           cfg::ASA3DMTConfig, start_model::String, dobs_filename::String,
                           derived::Vector{Pair{String,Any}}, columns::String)
    names = fieldnames(ASA3DMTConfig)
    hard = ["block_D" => _ASA_BLOCK, "T0" => _ASA_T0, "temp_ratio_scale" => _ASA_TEMP_RATIO_SCALE,
            "reanneal_every" => _ASA_REANNEAL_EVERY, "probes" => _ASA_PROBES,
            "reheat_max" => _ASA_REHEAT_MAX, "ref_quantile" => _ASA_REF_QUANTILE]
    w = maximum(length, vcat([String(f) for f in names], first.(hard), first.(derived), ["start_model"]))
    open(path, "w") do io
        println(io, "# ", title, " — ", timestamp)
        println(io, "# ", rpad("start_model", w), " = ", start_model)
        println(io, "# ", rpad("dobs", w), " = ", dobs_filename)
        for f in names
            println(io, "# ", rpad(String(f), w), " = ", getfield(cfg, f))
        end
        for (k, v) in hard
            println(io, "# ", rpad(k, w), " = ", v, "   (fixed)")
        end
        for (k, v) in derived
            println(io, "# ", rpad(k, w), " = ", v, "   (derived)")
        end
        println(io, repeat("-", length(columns)))
        println(io, columns)
        println(io, repeat("-", length(columns)))
    end
end

function _run_asa3d(start_model_path::String, dobs_filename::String, timestamp::String,
                    run_root::String; cfg::ASA3DMTConfig)
    rng = MersenneTwister(cfg.seed)
    detailed_log = joinpath(run_root, "0asa3DMT_detailed.log")
    iter_log     = joinpath(run_root, "0asa3DMT.log")
    reanneal_log = joinpath(run_root, "0asa3DMT_reanneal.log")

    m = load_ws3d_model(start_model_path)
    _, _, _, _, _, _, origin, rotation = read_ws3d_model(start_model_path, true)

    ix, iy = core_ranges(m; tol=cfg.pad_tol, expand=cfg.core_expand_cells)
    kz = cfg.z_core_cells > 0 ? (1:min(cfg.z_core_cells, m.nz)) :
         z_core_range(m, _z_core_max_depth(joinpath(run_root, dobs_filename), cfg))
    last(kz) < m.nz && @info @sprintf("z-core: perturbing %d of %d layers (bottom %.0f m)",
                                      length(kz), m.nz, m.cz[last(kz)])
    v0_core = Array(extract_core_array(m, ix, iy, kz))

    ring_values = vcat([vec(m.A[i, :, :]) for i in 1:m.nx if !(i in ix)]...,
                       [vec(m.A[:, j, :]) for j in 1:m.ny if !(j in iy)]...)
    background_log10 = median(filter(isfinite, ring_values))
    A_start = copy(m.A)

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
    any(air_full) && write_topography(joinpath(run_root, "topography.dat"), extract_topography(m);
                                      origin=m.origin)
    bathy !== nothing && !isempty(bathy.depth) &&
        write_bathymetry(joinpath(run_root, "bathymetry.dat"), bathy;
                         water_log10=cfg.water_log10, origin=m.origin)
    mask = BitArray(protected_full[ix, iy, kz])
    n_frozen = count(mask)

    rbfmap = build_rbf_map(m, ix, iy, cfg.n_ctrl, rng; kz=kz,
                           sigma_scale=cfg.sigma_scale, trunc_sigmas=cfg.trunc_sigmas,
                           sigma_scale_deep=cfg.sigma_scale_deep,
                           depth_power=cfg.ctrl_depth_power,
                           exclude=(n_frozen > 0 ? mask : nothing))
    M = length(rbfmap.ci)
    v0_ctrl = [v0_core[rbfmap.ci[t], rbfmap.cj[t], rbfmap.ck[t]] for t in 1:M]
    @info @sprintf("ASA: %d controls, %d frozen core cells", M, n_frozen)

    lo, hi = cfg.log_bounds
    B = min(_ASA_BLOCK, M)
    T0 = _ASA_T0
    m_scale = log(1 / _ASA_TEMP_RATIO_SCALE)
    k_f = max(1.0, cfg.max_iter * B / M)
    n_scale = log(k_f)
    c = m_scale * exp(-n_scale / B)

    delta_current = zeros(Float64, M)
    buffer = zeros(Float64, rbfmap.nx, rbfmap.ny, rbfmap.nz)
    k_gen = zeros(Int, M)
    reheat = ones(Float64, M)
    Tq = fill(T0, M)
    h_est = Float64[]

    function embed!(delta)
        apply_rbf_map!(buffer, rbfmap, delta)
        v = clamp.(v0_core .+ buffer, lo, hi)
        n_frozen > 0 && (v[mask] .= v0_core[mask])
        embed_core!(m, v, ix, iy, kz)
        smooth_padding_decay_z!(m, ix, iy, kz, background_log10, protected_full; target=A_start)
        smooth_padding_decay_xy!(m, ix, iy, background_log10, cfg.padding_decay_length,
                                 protected_full; target=A_start)
        return m
    end
    update_Tq!() = (@inbounds for q in 1:M
        Tq[q] = min(T0, reheat[q] * _asa_temperature(T0, c, k_gen[q], B))
    end)

    derived = Pair{String,Any}["n_controls" => M, "m" => m_scale, "k_f" => k_f, "n" => n_scale, "c" => c]
    hdr = (; timestamp=timestamp, cfg=cfg, start_model=start_model_path,
           dobs_filename=dobs_filename, derived=derived)
    _asa_write_header(iter_log, "ASA 3D MT iterations"; hdr...,
        columns=@sprintf("%8s %11s %11s %11s %11s %6s %7s %s",
                         "Iter", "TgenMed", "Tacc", "RMS", "RMSBest", "Nacc", "Nsolve", "Model"))
    _asa_write_header(detailed_log, "ASA 3D MT generated states"; hdr...,
        columns=@sprintf("%8s %11s %11s %11s %12s %10s %10s %5s %11s %s",
                         "Iter", "TgenBlk", "Tacc", "RMS", "dL", "Pacc", "Uacc", "Acc", "RMSBest", "Model"))
    _asa_write_header(reanneal_log, "ASA 3D MT reannealing"; hdr...,
        columns=@sprintf("%8s %7s %7s %6s %11s %10s %8s %11s %11s",
                         "Iter", "Nacc", "Nsolve", "Probes", "Href", "ReheatMed", "AtCap", "Tacc0", "Tacc"))

    embed!(delta_current)
    model0 = @sprintf("model_%03d_%02d.rho", 0, 0)
    rms0, fit0 = forward_and_misfit!(m; run_dir=run_root, model_filename=model0,
                                     dobs_filename=dobs_filename,
                                     dpred_filename=@sprintf("dpred_%03d_%02d.dat", 0, 0),
                                     origin=origin, rotation=rotation, cfg=cfg)
    n_solve = 1
    rms_current = rms0
    best_rms = rms0
    best_model_abs = joinpath(run_root, "best_model.rho")
    cp(joinpath(run_root, model0), best_model_abs; force=true)
    distortion_abs = joinpath(run_root, "distortion_best.txt")
    fit0 !== nothing && write_distortion_file(distortion_abs, fit0; iter=0, rms=rms0,
                                              damping=cfg.distortion_damping)
    current_model_rel = model0

    Tacc0 = rms0^2
    Tacc = Tacc0
    k_acc = 0.0
    n_acc = 0
    n_acc_last_reanneal = 0

    function reanneal!(k)
        embed!(delta_current)
        h, n_ok = _asa_gn_diagonal(m, rbfmap, ix, iy, kz, mask; run_dir=run_root,
                                   dobs_filename=dobs_filename, origin=origin, rotation=rotation,
                                   cfg=cfg, rng=rng)
        n_solve += n_ok
        href = NaN
        if h !== nothing
            h_est = isempty(h_est) ? h : 0.5 .* (h_est .+ h)
            w, href = _asa_reheat(h_est)
            reheat .= w
            update_Tq!()
        end
        if k > 0 && isfinite(rms_current)
            L_last, L_best = rms_current^2, best_rms^2
            Tacc0 = min(Tacc0, max(L_last, L_best, abs(L_last - L_best), eps()))
            Tacc = min(Tacc, Tacc0)
            k_acc = _asa_index(Tacc0, Tacc, c, B)
        end
        open(reanneal_log, "a") do io
            @printf(io, "%8d %7d %7d %6d %11.4e %10.3f %8.3f %11.4e %11.4e\n",
                    k, n_acc, n_solve, n_ok, href, median(reheat),
                    count(>=(_ASA_REHEAT_MAX), reheat) / M, Tacc0, Tacc)
        end
    end

    reanneal!(0)

    order = randperm(rng, M)
    pos = 1
    for k in 1:cfg.max_iter
        if pos + B - 1 > M
            order = randperm(rng, M)
            pos = 1
        end
        block = view(order, pos:pos+B-1)
        pos += B
        Tblk = median(Tq[block])

        delta_trial = copy(delta_current)
        _asa_generate!(delta_trial, block, Tq, v0_ctrl, lo, hi, rng)
        @inbounds for q in block
            k_gen[q] += 1
            Tq[q] = min(T0, reheat[q] * _asa_temperature(T0, c, k_gen[q], B))
        end
        embed!(delta_trial)
        model_filename = @sprintf("model_%03d_%02d.rho", k, 1)
        dpred_filename = @sprintf("dpred_%03d_%02d.dat", k, 1)
        rms_trial, fit_trial = forward_and_misfit!(m; run_dir=run_root, model_filename=model_filename,
                                                   dobs_filename=dobs_filename,
                                                   dpred_filename=dpred_filename,
                                                   origin=origin, rotation=rotation, cfg=cfg)
        n_solve += 1

        dL = rms_trial^2 - rms_current^2
        u = rand(rng)
        p_acc = !isfinite(dL) ? 0.0 : dL <= 0 ? 1.0 : exp(-dL / max(Tacc, 1e-300))
        accept = u < p_acc
        if accept
            delta_current .= delta_trial
            rms_current = rms_trial
            current_model_rel = model_filename
            n_acc += 1
            k_acc += 1
            Tacc = _asa_temperature(Tacc0, c, k_acc, B)
            if rms_current < best_rms
                best_rms = rms_current
                cp(joinpath(run_root, model_filename), best_model_abs; force=true)
                fit_trial !== nothing && write_distortion_file(distortion_abs, fit_trial; iter=k,
                                                               rms=best_rms,
                                                               damping=cfg.distortion_damping)
            end
        end

        keep = cfg.model_save_every > 0 ? (k % cfg.model_save_every == 0) : cfg.keep_models
        keep_d = cfg.model_save_every > 0 ? (k % cfg.model_save_every == 0) : cfg.keep_dpred
        keep || rm(joinpath(run_root, model_filename); force=true)
        keep_d || rm(joinpath(run_root, dpred_filename); force=true)

        open(detailed_log, "a") do io
            @printf(io, "%8d %11.4e %11.4e %11.5f %12.5f %10.5f %10.5f %5d %11.5f %s\n",
                    k, Tblk, Tacc, rms_trial, dL, p_acc, u, accept ? 1 : 0, best_rms, model_filename)
        end
        open(iter_log, "a") do io
            @printf(io, "%8d %11.4e %11.4e %11.5f %11.5f %6d %7d %s\n",
                    k, median(Tq), Tacc, rms_current, best_rms, n_acc, n_solve, current_model_rel)
        end

        if best_rms <= cfg.target_rms
            @info @sprintf("target RMS %.3f reached at iter %d (best RMS %.5f); stopping.",
                           cfg.target_rms, k, best_rms)
            break
        end
        if n_acc - n_acc_last_reanneal >= _ASA_REANNEAL_EVERY && k < cfg.max_iter
            n_acc_last_reanneal = n_acc
            reanneal!(k)
        end
    end
    return best_model_abs, iter_log
end

"""
    ASA3DMT(start_model_path; dobs_path, cfg)

Run a 3D adaptive simulated annealing MT inversion. Returns
`(best_model_path, iter_log_path)`.

The run directory `<out_root>_<yyyymmdd_HHMMSS>` holds `0asa3DMT.log` (state
after every generated model), `0asa3DMT_detailed.log` (every generated model
and its acceptance test), `0asa3DMT_reanneal.log` (every reannealing),
`best_model.rho`, and the retained trial models. One process runs one chain;
run several jobs with different seeds for an ensemble and combine them with
[`AnalyseEnsemble3D`](@ref).

# Arguments
- `start_model_path::AbstractString`: WS3D-format starting model (.rho/.ws)
- `dobs_path::AbstractString`: observed data file
- `cfg::ASA3DMTConfig`: configuration (see `ASA3DMTConfig`)
"""
function ASA3DMT(start_model_path::AbstractString; dobs_path::AbstractString, cfg::ASA3DMTConfig)
    _asa_check_cfg(cfg)
    start_model_abs = abspath(start_model_path)
    run_root_base = isabspath(cfg.out_root) ? cfg.out_root :
                    joinpath(dirname(start_model_abs), cfg.out_root)
    t_now = now()
    run_root = string(run_root_base, "_", Dates.format(t_now, "yyyymmdd_HHMMSS"))
    isdir(run_root) || mkpath(run_root)
    dobs_filename = basename(dobs_path)
    cp(abspath(dobs_path), joinpath(run_root, dobs_filename); force=true)
    isfile(PROGRAM_FILE) &&
        cp(abspath(PROGRAM_FILE), joinpath(run_root, basename(PROGRAM_FILE)); force=true)
    timestamp = Dates.format(t_now, "yyyy-mm-dd HH:MM:SS")
    return _run_asa3d(start_model_abs, dobs_filename, timestamp, run_root; cfg=cfg)
end

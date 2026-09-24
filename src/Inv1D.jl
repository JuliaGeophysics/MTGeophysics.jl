# 1D MT inversion
# Author: @pankajkmishra
# Every site of a data file inverted on its own, on the layering MakeMesh1D gives it, from its background
# halfspace: damped Gauss-Newton with vertical smoothing, or VFSA with every layer a parameter (no sparse
# parameterisation). Self-contained: only the data and model file formats and the plot style come from 2D.
# One small control file, the rest fixed

using CairoMakie
using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

#---------- control ----------

"""
    InvCtrl1D

1D inversion control, read from `InvCtrl.GN` or `InvCtrl.VFSA` (`examples/ctrl/1D`).
- `algorithm`: `:gn` or `:vfsa`
- `mode`: impedance fitted, `:XY`, `:YX`, `:XYYX` or `:DET`
- `target_rms`, `max_iter`: stopping controls
- `lambda`: GN regularization weight (vertical smoothing, small smallness), fixed
- `chains`, `log_bounds`: VFSA chains and the log10 ρ search box; every layer is a VFSA
  parameter. GN is unbounded

Fixed inside: GN damping 0.01, smallness 0.01, max log10 step 0.5; VFSA temperature 0.03,
cooling ratio 0.001, a fifth of the layers moved per proposal, step scale 0.11.
"""
Base.@kwdef struct InvCtrl1D
    algorithm::Symbol
    target_rms::Float64
    max_iter::Int
    mode::Symbol = :XYYX
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    lambda::Float64 = 1.0
    chains::Int = 10
end

_ctrl_algorithm1d(s) = (a = Symbol(lowercase(strip(s))); a in (:gn, :vfsa) ? a :
                        throw(ArgumentError("1D algorithm must be GN or VFSA")))

const _INV1D_CTRL_SPEC = (
    "Algorithm"                          => (:algorithm, _ctrl_algorithm1d, true),
    "Mode"                               => (:mode, _ctrl_mode1d, false),
    "Exit search when rms is less than"  => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"       => (:max_iter, _ctrl_int, true),
    "Log10 resistivity bounds"           => (:log_bounds, _ctrl_parse_pair, false),
    "Initial damping factor lambda"      => (:lambda, _ctrl_float, false),
    "Number of chains"                   => (:chains, _ctrl_int, false),
)

function _validate_ctrl(c::InvCtrl1D)
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c.log_bounds[1] < c.log_bounds[2] || throw(ArgumentError("log10 bounds must be increasing"))
    c.lambda >= 0 || throw(ArgumentError("lambda must be nonnegative"))
    c.chains >= 1 || throw(ArgumentError("the number of chains must be at least 1"))
    c
end

"""
    ReadInvCtrl1D(path) -> InvCtrl1D

Read a 1D inversion control file. GN files may not carry the VFSA keys (chains, bounds),
nor VFSA files lambda.
"""
function ReadInvCtrl1D(path::AbstractString)
    v = _read_ctrl(path, _INV1D_CTRL_SPEC, "1D inversion control")
    other = v[:algorithm] == :gn ? (:chains, :log_bounds) : (:lambda,)
    foreign = filter(k -> haskey(v, k), other)
    isempty(foreign) || error("$path: $(join(foreign, ", ")) do not apply to $(uppercase(string(v[:algorithm])))")
    _validate_ctrl(InvCtrl1D(; v...))
end

"""
    WriteInvCtrl1D(path, ctrl::InvCtrl1D) -> path

Write a 1D inversion control file with its algorithm's keys.
"""
function WriteInvCtrl1D(path::AbstractString, c::InvCtrl1D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g(x) = @sprintf("%.6g", x)
    rows = [("Algorithm", uppercase(string(c.algorithm))), ("Mode", string(c.mode)),
            ("Exit search when rms is less than", g(c.target_rms)), ("Maximum number of iterations", string(c.max_iter))]
    c.algorithm == :gn ? push!(rows, ("Initial damping factor lambda", g(c.lambda))) :
        append!(rows, [("Log10 resistivity bounds", "$(g(c.log_bounds[1])) $(g(c.log_bounds[2]))"),
                       ("Number of chains", string(c.chains))])
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

#---------- misfit ----------

# the fitted data of one site: impedance, sign of Z, error; Re and Im count separately
function _inv1d_rows(site::DataFile2D, mode::Symbol)
    rows = [(k = k, s = s, d = getproperty(site, key)[k, 1], σ = getproperty(site, ekey)[k, 1])
            for (key, ekey, s) in _mt1d_components(mode) for k in eachindex(site.frequencies)]
    filter!(r -> isfinite(r.d) && isfinite(r.σ) && r.σ > 0, rows)
    isempty(rows) && throw(ArgumentError("site $(site.site_names[1]) has no valid impedances for mode $mode"))
    rows
end

# weighted residual (Re, Im per datum) of log10 model m
function _inv1d_residual(m, h, f, rows)
    Z = mt1d_impedance(f, 10 .^ m, h)
    r = [(row.s * Z[row.k] - row.d) / row.σ for row in rows]
    vcat(real.(r), imag.(r))
end

_inv1d_rms(r) = sqrt(sum(abs2, r) / length(r))

#---------- gauss-newton ----------

# vertical first differences weighted by 1/√(mean thickness) and a small smallness weighted by √h,
# both relative to the median thickness, so λ has the meaning it has in 2D
function _inv1d_regularizer(h, smallness)
    n, h0 = length(h), median(h)
    R = zeros(2n - 1, n)
    for i in 1:n
        R[i, i] = sqrt(smallness * h[i] / h0)
    end
    for i in 1:n-1
        w = sqrt(h0 / ((h[i] + h[i+1]) / 2))
        R[n+i, i], R[n+i, i+1] = -w, w
    end
    R
end

function _inv1d_gn(h, f, rows, m0, ctrl; name = "")
    λ = ctrl.lambda
    R = _inv1d_regularizer(h, 0.01)
    φ(m, r) = (sum(abs2, r) + λ * sum(abs2, R * (m - m0))) / 2
    m, r = copy(m0), _inv1d_residual(m0, h, f, rows)
    obj, μ = φ(m, r), 1e-2
    history = [(iteration = 0, rms = _inv1d_rms(r), objective = obj, step = 0.0, damping = μ)]
    reason = :max_iter
    for it in 1:ctrl.max_iter
        _inv1d_rms(r) <= ctrl.target_rms && (reason = :target_rms; break)
        J = ForwardDiff.jacobian(x -> _inv1d_residual(x, h, f, rows), m)
        H = J' * J + λ * (R' * R)
        g = J' * r + λ * (R' * (R * (m - m0)))
        accepted = false
        for _ in 1:6
            δ = -((H + μ * Diagonal(max.(diag(H), 1e-12))) \ g)
            δ .*= min(1.0, 0.5 / max(maximum(abs, δ), eps()))
            mt = m .+ δ
            rt = _inv1d_residual(mt, h, f, rows)
            if φ(mt, rt) < obj
                step, m, r, obj_old = maximum(abs, mt - m), mt, rt, obj
                obj, μ, accepted = φ(m, r), max(μ / 3, 1e-12), true
                push!(history, (iteration = it, rms = _inv1d_rms(r), objective = obj, step, damping = μ))
                obj_old - obj <= 1e-6 * obj_old && (reason = :stalled)
                break
            end
            μ *= 10
        end
        accepted || (reason = :line_search_failed; break)
        reason == :stalled && break
        _inv1d_rms(r) <= ctrl.target_rms && (reason = :target_rms; break)
    end
    _inv1d_rms(r) <= ctrl.target_rms && (reason = :target_rms)
    @printf("GN 1D %-10s %3d iterations  RMS %.4f  %s\n", name, history[end].iteration, history[end].rms, reason)
    (; m, history, reason)
end

#---------- vfsa ----------

# Ingber's proposal on a random share of the layers, each redrawn into the log10 box
function _inv1d_propose(m, T, lo, hi, span, share, rng)
    t = copy(m)
    for j in randperm(rng, length(m))[1:max(1, round(Int, share * length(m)))]
        for _ in 1:100
            u = rand(rng)
            c = m[j] + sign(u - 0.5) * T * ((1 + 1 / T)^abs(2u - 1) - 1) * span
            lo <= c <= hi && (t[j] = c; break)
        end
    end
    t
end

function _inv1d_vfsa_chain(k, h, f, rows, m0, ctrl, rng)
    # tuned on 1D-I (10 chains × 2000): best chain 0.74, median chain 0.83; T0 1 or all layers at once fit worse
    T0, cool, step, share = 0.03, 1e-3, 0.11, 0.2
    lo, hi = ctrl.log_bounds
    rate = log(1 / cool) / max(ctrl.max_iter - 1, 1)
    m, rms = copy(m0), _inv1d_rms(_inv1d_residual(m0, h, f, rows))
    best, best_rms = copy(m), rms
    history = [(chain = k, iteration = 0, temperature = T0, trial_rms = rms, rms, best_rms, accepted = true)]
    for it in 1:ctrl.max_iter
        best_rms <= ctrl.target_rms && break
        T = T0 * exp(-rate * (it - 1))
        trial = _inv1d_propose(m, T, lo, hi, (hi - lo) * step, share, rng)
        trms = _inv1d_rms(_inv1d_residual(trial, h, f, rows))
        dE = (trms^2 - rms^2) / max(rms^2, eps())
        accepted = isfinite(trms) && rand(rng) < (dE <= 0 ? 1.0 : exp(-dE / T))
        accepted && ((m, rms) = (trial, trms))
        accepted && rms < best_rms && ((best, best_rms) = (copy(m), rms))
        push!(history, (chain = k, iteration = it, temperature = T, trial_rms = trms, rms, best_rms, accepted))
    end
    (; chain = k, best, best_rms, history, acceptance = length(history) > 1 ? mean(x.accepted for x in history[2:end]) : NaN)
end

# chains on threads; the ensemble of their best models in log10 ρ
function _inv1d_vfsa(h, f, rows, m0, ctrl; name = "")
    chains = Vector{Any}(undef, ctrl.chains)
    Threads.@threads :dynamic for k in 1:ctrl.chains
        chains[k] = _inv1d_vfsa_chain(k, h, f, rows, m0, ctrl, MersenneTwister(20260308 + 1000 * (k - 1)))
    end
    L = reduce(hcat, [c.best for c in chains])
    q(p) = [quantile(L[i, :], p) for i in axes(L, 1)]
    ensemble = (mean = vec(mean(L; dims = 2)), median = vec(median(L; dims = 2)),
                std = size(L, 2) > 1 ? vec(std(L; dims = 2)) : zeros(size(L, 1)), p05 = q(0.05), p95 = q(0.95),
                count = size(L, 2))
    b = argmin([c.best_rms for c in chains])
    @printf("VFSA 1D %-10s %d chains  best chain %d RMS %.4f\n", name, length(chains), b, chains[b].best_rms)
    (; chains = identity.(chains), ensemble, best_chain = b, best = chains[b].best, best_rms = chains[b].best_rms)
end

#---------- files ----------

function _inv1d_write_csv(path, history)
    open(path, "w") do io
        println(io, join(string.(keys(first(history))), ","))
        foreach(x -> println(io, join([v isa Bool ? Int(v) : v isa Symbol ? string(v) : @sprintf("%.6g", v) for v in values(x)], ",")),
                history)
    end
    path
end

function _inv1d_open_run(run_dir, data_path, inputs)
    dir = run_dir === nothing ? joinpath(dirname(abspath(data_path)), "run_" * Dates.format(now(), "yyyymmdd_HHMMSS")) :
          String(run_dir)
    mkpath(joinpath(dir, "inputs"))
    foreach(p -> cp(p, joinpath(dir, "inputs", basename(p)); force = true), inputs)
    println("Run directory: ", dir)
    dir
end

_inv1d_write_model(path, h, m) = WriteModel2D(path, [1.0], h, reshape(10 .^ m, :, 1))

#---------- inversion ----------

"""
    Invert1D(data_path, inv_path, meshes=MakeMesh1D(data_path); sites=nothing, run_dir=nothing)

1D inversion of every site of the data file (or those named in `sites`) on its own, on its
`MakeMesh1D` layering, from its background halfspace (also the GN reference model). The
control (`InvCtrl1D`) picks GN or VFSA and the impedance. Writes into `run_dir` (default
`run_YYYYmmdd_HHMMSS/` next to the data) one folder per site with `model.start`,
`model.rho`, `data.pred` and `History.csv` (GN) or `vfsa/` (VFSA: per-chain history, the
ensemble mean, median, 5 and 95% models, `model.rho` = ensemble mean), plus `data.pred`
of all sites, `Summary.txt` and the inputs. Returns the run directory, the per-site
results and the overall rms.
"""
function Invert1D(data_path::AbstractString, inv_path::AbstractString, meshes = nothing;
                  sites = nothing, run_dir::Union{Nothing, AbstractString} = nothing)
    observed = load_data2d(data_path)
    ctrl = ReadInvCtrl1D(inv_path)
    meshes = something(meshes, MakeMesh1D(observed; mode = ctrl.mode))
    picked = sites === nothing ? collect(eachindex(observed.site_names)) :
             [something(findfirst(==(String(s)), observed.site_names), 0) for s in sites]
    all(>(0), picked) || error("unknown sites $(sites[picked .== 0])")
    lookup = Dict(m.site => m for m in meshes)
    missing_mesh = setdiff(observed.site_names[picked], keys(lookup))
    isempty(missing_mesh) || error("no mesh for sites $missing_mesh")
    lo, hi = ctrl.log_bounds
    f = observed.frequencies
    dir = _inv1d_open_run(run_dir, data_path, (data_path, inv_path))

    results = map(picked) do i
        site = mt1d_site_data(observed, i; mode = ctrl.mode)
        name = site.site_names[1]
        h = lookup[name].thicknesses
        rows = _inv1d_rows(site, ctrl.mode)
        ρ0 = log10(lookup[name].background)
        m0 = fill(ctrl.algorithm == :vfsa ? clamp(ρ0, lo, hi) : ρ0, length(h))
        sdir = mkpath(joinpath(dir, name))
        _inv1d_write_model(joinpath(sdir, "model.start"), h, m0)
        history, vfsa = nothing, nothing
        if ctrl.algorithm == :gn
            g = _inv1d_gn(h, f, rows, m0, ctrl; name)
            m, history, reason = g.m, g.history, g.reason
            _inv1d_write_csv(joinpath(sdir, "History.csv"), history)
        else
            vfsa = _inv1d_vfsa(h, f, rows, m0, ctrl; name)
            m = vfsa.ensemble.mean
            vdir = mkpath(joinpath(sdir, "vfsa"))
            for c in vfsa.chains
                _inv1d_write_csv(joinpath(vdir, @sprintf("History_chain_%02d.csv", c.chain)), c.history)
            end
            for k in (:mean, :median, :p05, :p95)
                _inv1d_write_model(joinpath(vdir, "model.$k.rho"), h, getproperty(vfsa.ensemble, k))
            end
            _inv1d_write_model(joinpath(vdir, "model.best.rho"), h, vfsa.best)
            best = _mt1d_predicted(site, reshape(mt1d_impedance(f, 10 .^ vfsa.best, h), :, 1); fractional = false)
            write_data2d(joinpath(vdir, "data.best.pred"), best)
        end
        predicted = _mt1d_predicted(site, reshape(mt1d_impedance(f, 10 .^ m, h), :, 1); fractional = false)
        ctrl.mode == :XY && (predicted.z_yx .= NaN)
        ctrl.mode == :YX && (predicted.z_xy .= NaN)
        ctrl.mode == :DET && (predicted.z_yx .= NaN)
        r = _inv1d_residual(m, h, f, rows)
        rms = _inv1d_rms(r)
        converged = rms <= ctrl.target_rms
        vfsa === nothing || (reason = converged ? :target_rms : :max_iter)
        _inv1d_write_model(joinpath(sdir, "model.rho"), h, m)
        write_data2d(joinpath(sdir, "data.pred"), predicted)
        (; site = name, index = i, run_dir = sdir, thicknesses = h, observed = site, predicted, start = 10 .^ m0,
           final = 10 .^ m, history, vfsa, chi2 = sum(abs2, r), count = length(r), rms, reason, converged, ctrl)
    end

    # all sites in one file, at their own positions
    Z = reduce(hcat, [r.predicted.z_xy for r in results])
    survey = DataFile2D(title = observed.title, periods = observed.periods, frequencies = f,
                        site_names = observed.site_names[picked], receivers = observed.receivers[picked],
                        x_positions = observed.x_positions[picked], z_positions = observed.z_positions[picked],
                        z_xy = zeros(ComplexF64, size(Z)), z_xy_error = reduce(hcat, [r.predicted.z_xy_error for r in results]),
                        z_yx = zeros(ComplexF64, size(Z)), z_yx_error = reduce(hcat, [r.predicted.z_yx_error for r in results]),
                        z_xx = fill(complex(NaN), size(Z)), z_xx_error = fill(NaN, size(Z)),
                        z_yy = fill(complex(NaN), size(Z)), z_yy_error = fill(NaN, size(Z)),
                        rho_xy = zeros(size(Z)), phase_xy = zeros(size(Z)), rho_yx = zeros(size(Z)), phase_yx = zeros(size(Z)),
                        latitudes = isempty(observed.latitudes) ? Float64[] : observed.latitudes[picked],
                        longitudes = isempty(observed.longitudes) ? Float64[] : observed.longitudes[picked],
                        origin = observed.origin, rotation = observed.rotation,
                        rotations = observed.rotations)
    predicted = _mt1d_predicted(survey, reduce(hcat, [mt1d_impedance(f, r.final, r.thicknesses) for r in results]); fractional = false)
    ctrl.mode in (:XY, :DET) && (predicted.z_yx .= NaN)
    ctrl.mode == :YX && (predicted.z_xy .= NaN)
    write_data2d(joinpath(dir, "data.pred"), predicted)
    chi2, ndata = sum(r.chi2 for r in results), sum(r.count for r in results)
    rms = sqrt(chi2 / ndata)
    open(joinpath(dir, "Summary.txt"), "w") do io
        println(io, "Algorithm: ", uppercase(string(ctrl.algorithm)))
        println(io, "Impedance: ", ctrl.mode)
        @printf(io, "RMS: %.6f\n", rms)
        println(io, "Real data count: ", ndata)
        for r in results
            @printf(io, "  %-12s RMS %.6f  %s  %d layers, start %.1f ohm m%s\n", r.site, r.rms, r.reason,
                    length(r.thicknesses), r.start[1], r.history === nothing ?
                    @sprintf("  ensemble mean; best chain %d RMS %.6f", r.vfsa.best_chain, r.vfsa.best_rms) :
                    "  $(r.history[end].iteration) iterations")
        end
    end
    (; run_dir = dir, algorithm = ctrl.algorithm, ctrl, sites = results, predicted, rms)
end

#---------- plots ----------

"""
    plot_mt1d_model(z_cell_sizes, models; output_path, band=nothing, maximum_depth=nothing) -> path

Resistivity-depth steps of one or more 1D models on log axes. `models` holds
`(resistivity, label, color, linestyle)` tuples; `band = (low, high)` shades a range,
e.g. the VFSA 5-95% range. The last layer is drawn down to `maximum_depth` (default
twice its top).
"""
function plot_mt1d_model(z_cell_sizes::AbstractVector{<:Real}, models; output_path::AbstractString,
                         band = nothing, maximum_depth = nothing)
    CairoMakie.activate!()
    tops = vcat(0.0, cumsum(z_cell_sizes[1:end-1]))
    bottom = something(maximum_depth, 2 * max(tops[end], z_cell_sizes[1]))
    ztop = z_cell_sizes[1] / 2
    edges(z) = vcat(max(ztop, 1e-3), repeat(max.(z[2:end], ztop), inner = 2), bottom)
    steps(ρ) = repeat(vec(ρ), inner = 2)
    figure = Figure(size = (600, 750))
    axis = _mt_axis(figure[1, 1]; xlabel = "Resistivity (Ω·m)", ylabel = "Depth (m)", xscale = log10, yscale = log10,
                    yreversed = true)
    if band !== nothing
        z = edges(tops)
        band!(axis, Point2f.(steps(band[1]), z), Point2f.(steps(band[2]), z); color = (:steelblue, 0.25), label = "5-95%")
    end
    for (ρ, label, color, style) in models
        lines!(axis, steps(ρ), edges(tops); color, linestyle = style, linewidth = 2, label)
    end
    ylims!(axis, bottom, ztop)
    axislegend(axis, position = :lb, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

"""
    PlotModel1D(model_path; output_path, true_model_path=nothing, maximum_depth=nothing) -> path

Plot a one-column model file, with the true model when given.
"""
function PlotModel1D(model_path::AbstractString; output_path::AbstractString,
                     true_model_path::Union{Nothing, AbstractString} = nothing, maximum_depth = nothing)
    m = ReadModel2D(model_path)
    models = Any[(m.resistivity[:, 1], "model", :black, :solid)]
    if true_model_path !== nothing
        t = ReadModel2D(true_model_path)
        push!(models, (_mt1d_resample(t, m.z_cell_sizes), "true", :firebrick, :solid))
    end
    plot_mt1d_model(m.z_cell_sizes, models; output_path, maximum_depth)
end

# a 1D model's resistivity at the layer centres of another layering
function _mt1d_resample(model::ModelFile2D, dz::AbstractVector{<:Real})
    tops = cumsum(model.z_cell_sizes)
    centres = cumsum(dz) .- dz ./ 2
    [model.resistivity[min(searchsortedfirst(tops, z), length(tops)), 1] for z in centres]
end

"""
    plot_mt1d_convergence(histories; output_path, target_rms=0.0) -> path

RMS per iteration of a GN history, or current (thin) and best (thick) RMS of each VFSA chain.
"""
function plot_mt1d_convergence(histories::AbstractVector; output_path::AbstractString, target_rms::Real = 0.0)
    CairoMakie.activate!()
    figure = Figure(size = (800, 450))
    axis = _mt_axis(figure[1, 1]; xlabel = "Iteration", ylabel = "RMS", yscale = log10)
    for h in histories
        it = [x.iteration for x in h]
        if hasproperty(first(h), :best_rms)
            lines!(axis, it, [x.rms for x in h]; color = (:gray, 0.4), linewidth = 0.8)
            lines!(axis, it, [x.best_rms for x in h]; color = :navy, linewidth = 1.8)
        else
            scatterlines!(axis, it, [x.rms for x in h]; color = :navy)
        end
    end
    target_rms > 0 && hlines!(axis, [target_rms], color = :gray, linestyle = :dash)
    _mt_save(output_path, figure)
end

"""
    PlotInversion1D(run; true_model_path=nothing, maximum_depth=nothing) -> paths

Standard plots of a run returned by `Invert1D`, into each site's `plots/`: the start,
final and (when given) true models, with the VFSA median and 5-95% range; the data fit;
the convergence.
"""
function PlotInversion1D(run; true_model_path::Union{Nothing, AbstractString} = nothing, maximum_depth = nothing)
    paths = String[]
    truth = true_model_path === nothing ? nothing : ReadModel2D(true_model_path)
    for r in run.sites
        path(n) = joinpath(r.run_dir, "plots", n)
        dz = r.thicknesses
        models = Any[(r.start, "start", :gray, :dash), (r.final, r.vfsa === nothing ? "final" : "mean", :black, :solid)]
        band = nothing
        if r.vfsa !== nothing
            e = r.vfsa.ensemble
            push!(models, (10 .^ e.median, "median", :steelblue, :solid))
            band = (10 .^ e.p05, 10 .^ e.p95)
        end
        truth === nothing || push!(models, (_mt1d_resample(truth, dz), "true", :firebrick, :solid))
        push!(paths, plot_mt1d_model(dz, models; output_path = path("Model.png"), band, maximum_depth))
        push!(paths, plot_mt2d_data_fit(r.observed, r.predicted; output_path = path("DataFit.png"),
                                        names = _mt1d_plot_names(r.ctrl.mode)))
        histories = r.vfsa === nothing ? [r.history] : [c.history for c in r.vfsa.chains]
        push!(paths, plot_mt1d_convergence(histories; output_path = path(r.vfsa === nothing ? "Convergence.png" : "ConvergenceVFSA.png"),
                                           target_rms = r.ctrl.target_rms))
    end
    paths
end

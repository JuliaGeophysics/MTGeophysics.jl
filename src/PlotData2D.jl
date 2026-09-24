# 2D MT data plots
# Author: @pankajkmishra
# Apparent resistivity and phase curves (observed with error bars, predicted as lines), inversion
# convergence, and PlotInversion2D, the standard plots of a run

using CairoMakie

#---------- style ----------

const _MT_COLOURS = (TE = RGBf(0.12, 0.38, 0.72), TM = RGBf(0.84, 0.30, 0.10))
const _MT_MARKER = (marker = :circle, markersize = 10, strokecolor = :black, strokewidth = 1.6)
const _MT_LINEWIDTH = 2.6
const _MT_PX_PER_UNIT = 3

# no grid, no title, regular-weight labels
_mt_axis(pos; kwargs...) = Axis(pos; xgridvisible = false, ygridvisible = false,
                                xlabelfont = :regular, ylabelfont = :regular, kwargs...)

# plain-number labels on log axes, 63 rather than 10^1.8
_mt_log_labels(values) = [v >= 10 ? @sprintf("%.0f", v) : @sprintf("%.2g", v) for v in values]

_mt_save(path, figure) = (mkpath(dirname(abspath(path))); save(path, figure; px_per_unit = _MT_PX_PER_UNIT); String(path))

#---------- curves ----------

# one rho/phase pair per station, observed with error bars, predicted as lines; `names` labels the xy and
# yx curves (TE and TM in 2D, XY/YX or DET in 1D), `nothing` leaves a component out
function _plot_mt2d_curves(observed::Union{Nothing, DataFile2D}, predicted::Union{Nothing, DataFile2D},
                           stations::AbstractVector{<:Integer}, output_path::AbstractString;
                           names = (TE = "TE", TM = "TM"))
    CairoMakie.activate!()
    ref = something(observed, predicted)
    T = ref.periods
    ncol = min(4, length(stations))
    nrow = cld(length(stations), ncol)
    figure = Figure(size = (320 * ncol, 460 * nrow))
    for (n, s) in enumerate(stations)
        r, c = fldmod1(n, ncol)
        cell = figure[r, c] = GridLayout()
        ax_ρ = _mt_axis(cell[1, 1]; xscale = log10, yscale = log10, ytickformat = _mt_log_labels,
                        ylabel = c == 1 ? "Apparent resistivity (Ω·m)" : "", xticklabelsvisible = false)
        ax_φ = _mt_axis(cell[2, 1]; xscale = log10, xlabel = r == nrow ? "Period (s)" : "",
                        ylabel = c == 1 ? "Phase (°)" : "")
        rowgap!(cell, 4)
        linkxaxes!(ax_ρ, ax_φ)
        for (pol, key) in ((:TE, "xy"), (:TM, "yx"))
            names[pol] === nothing && continue
            colour = _MT_COLOURS[pol]
            if predicted !== nothing
                ρp = getproperty(predicted, Symbol("rho_", key))[:, s]
                φp = getproperty(predicted, Symbol("phase_", key))[:, s]
                pol == :TM && (φp = _phase_fold_to_0_90(φp))
                lines!(ax_ρ, T, max.(ρp, eps()), color = colour, linewidth = _MT_LINEWIDTH)
                lines!(ax_φ, T, φp, color = colour, linewidth = _MT_LINEWIDTH)
            end
            observed === nothing && continue
            z = getproperty(observed, Symbol("z_", key))[:, s]
            ρo = getproperty(observed, Symbol("rho_", key))[:, s]
            φo = getproperty(observed, Symbol("phase_", key))[:, s]
            pol == :TM && (φo = _phase_fold_to_0_90(φo))
            ok = isfinite.(ρo) .& (ρo .> 0)
            any(ok) || continue
            # first order: δρ/ρ = 2|δZ|/|Z|, δφ = |δZ|/|Z| rad
            rel = getproperty(observed, Symbol("z_", key, "_error"))[:, s] ./ abs.(z)
            has = ok .& isfinite.(rel)
            if any(has)
                δρ = 2 .* ρo[has] .* rel[has]
                errorbars!(ax_ρ, T[has], ρo[has], min.(δρ, 0.95 .* ρo[has]), δρ,
                           color = :black, linewidth = 1, whiskerwidth = 6)
                errorbars!(ax_φ, T[has], φo[has], rad2deg.(rel[has]), color = :black, linewidth = 1, whiskerwidth = 6)
            end
            scatter!(ax_ρ, T[ok], ρo[ok]; color = colour, _MT_MARKER...)
            scatter!(ax_φ, T[ok], φo[ok]; color = colour, _MT_MARKER...)
        end
        # phases rarely reach the bottom of 0-90, so the site name sits there
        text!(ax_φ, 0.97, 0.04; text = ref.site_names[s], space = :relative, align = (:right, :bottom), fontsize = 13)
        ylims!(ax_φ, 0, 90)
    end
    elements, labels = Any[], String[]
    for pol in (:TE, :TM)
        names[pol] === nothing && continue
        if observed !== nothing
            push!(elements, MarkerElement(; color = _MT_COLOURS[pol], _MT_MARKER...))
            push!(labels, "$(names[pol]) observed")
        end
        if predicted !== nothing
            push!(elements, LineElement(color = _MT_COLOURS[pol], linewidth = _MT_LINEWIDTH))
            push!(labels, "$(names[pol]) predicted")
        end
    end
    Legend(figure[0, :], elements, labels; orientation = :horizontal, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

"""
    plot_mt2d_data_fit(observed, predicted; output_path, station_indices=nothing, names=(TE = "TE", TM = "TM")) -> path

Observed apparent resistivity and phase (circles with error bars) against the
predicted curves (lines), one panel pair per station; all stations by default. `names`
labels the ZXY and ZYX curves (`nothing` leaves one out).
"""
function plot_mt2d_data_fit(observed::DataFile2D, predicted::DataFile2D;
                            output_path::AbstractString,
                            station_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing,
                            names = (TE = "TE", TM = "TM"))
    size(observed.z_xy) == size(predicted.z_xy) || throw(DimensionMismatch("observed and predicted surveys differ"))
    _plot_mt2d_curves(observed, predicted, something(station_indices, 1:length(observed.receivers)), output_path; names)
end

"""
    PlotData2D(data_path; output_path, predicted_path=nothing, station_indices=nothing, names=(TE = "TE", TM = "TM")) -> path

Apparent resistivity and phase of a data file, with error bars, and the curves of
`predicted_path` on top when given. `names` labels the ZXY and ZYX curves, e.g.
`(TE = "XY", TM = "YX")` for 1D.
"""
function PlotData2D(data_path::AbstractString;
                    output_path::AbstractString,
                    predicted_path::Union{Nothing, AbstractString} = nothing,
                    station_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing,
                    names = (TE = "TE", TM = "TM"))
    observed = load_data2d(data_path)
    predicted = predicted_path === nothing ? nothing : load_data2d(predicted_path)
    _plot_mt2d_curves(observed, predicted, something(station_indices, 1:length(observed.receivers)), output_path; names)
end

"""
    plot_mt2d_site_curves(response; station_index, output_path) -> path

Apparent resistivity and phase of one station of a computed response, as lines.
"""
function plot_mt2d_site_curves(response::MT2DResponse;
                               station_index::Integer = cld(length(response.receivers), 2),
                               output_path::AbstractString)
    n = length(response.receivers)
    data = data_from_response2d(response; site_names = [@sprintf("Station %d", i) for i in 1:n])
    _plot_mt2d_curves(nothing, data, [station_index], output_path)
end

#---------- inversion convergence ----------

"""
    plot_inv2d_convergence(history; output_path, target_rms=0.0) -> path

Rms (with the target line when `target_rms > 0`) and the objective terms per iteration.
"""
function plot_inv2d_convergence(history::AbstractVector; output_path::AbstractString, target_rms::Real = 0.0)
    CairoMakie.activate!()

    it = [h.iteration for h in history]
    figure = Figure(size = (1100, 450))
    ax1 = _mt_axis(figure[1, 1]; xlabel = "Iteration", ylabel = "RMS")
    scatterlines!(ax1, it, [h.rms for h in history], color = :navy)
    target_rms > 0 && hlines!(ax1, [target_rms], color = :gray, linestyle = :dash)

    ax2 = _mt_axis(figure[1, 2]; xlabel = "Iteration", ylabel = "Objective terms", yscale = log10)
    scatterlines!(ax2, it, [h.objective for h in history], label = "objective")
    scatterlines!(ax2, it, [h.chi2 / 2 for h in history], label = "χ²/2")
    # regularization is exactly zero at the start model, skip it on a log axis
    reg = [(h.iteration, h.regularization) for h in history if h.regularization > 0]
    isempty(reg) || scatterlines!(ax2, first.(reg), last.(reg), label = "‖R(m - m_ref)‖²/2")
    axislegend(ax2, position = :rt, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

"""
    plot_vfsa2d_convergence(histories; output_path, target_rms=0.0) -> path

Current (thin) and best (thick) rms of each VFSA chain per iteration.
"""
function plot_vfsa2d_convergence(histories::AbstractVector; output_path::AbstractString, target_rms::Real = 0.0)
    CairoMakie.activate!()
    figure = Figure(size = (1100, 500))
    axis = _mt_axis(figure[1, 1]; xlabel = "Iteration", ylabel = "RMS", yscale = log10)
    colors = Makie.wong_colors()
    for (k, h) in enumerate(histories)
        it = [r.iteration for r in h]
        c = colors[mod1(k, length(colors))]
        lines!(axis, it, [r.rms for r in h], color = (c, 0.35), linewidth = 0.8)
        lines!(axis, it, [r.best_rms for r in h], color = c, linewidth = 2, label = "chain $(h[1].chain)")
    end
    target_rms > 0 && hlines!(axis, [target_rms], color = :gray, linestyle = :dash)
    length(histories) <= 12 && Legend(figure[1, 2], axis, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

#---------- inversion run ----------

"""
    PlotInversion2D(run; true_model_path=nothing, maximum_depth_km=10.0,
                    resistivity_log10_range=(0.0, 4.0), background_resistivity=100.0) -> paths

Standard plots of a run returned by the six-file `Invert2D` or the five-file `VFSA2D`, written to `run.run_dir/plots`:
the mesh, the start, final and (when given) true models, the data fit, and the
convergence for GN and NLCG.
"""
function PlotInversion2D(run; true_model_path::Union{Nothing, AbstractString} = nothing,
                         maximum_depth_km::Real = 10.0, resistivity_log10_range = (0.0, 4.0),
                         background_resistivity::Real = 100.0)
    dir = joinpath(run.run_dir, "plots")
    path(name) = joinpath(dir, name)
    water = hasproperty(run, :water) && any(run.water) ? run.water : nothing
    core = (show_padding = false, maximum_depth_km, resistivity_log10_range, water)
    paths = String[
        plot_mt2d_mesh(run.mesh; output_path = path("Mesh.png"), region = :full, background_resistivity),
        plot_mt2d_mesh(run.mesh; output_path = path("MeshCore.png"), region = :core, background_resistivity),
        plot_mt2d_model(run.mesh, run.start; output_path = path("ModelStart.png"), core...),
        plot_mt2d_model(run.mesh, run.final; output_path = path("ModelFinal.png"), core...),
        plot_mt2d_model(run.mesh, run.final; output_path = path("ModelFinalFull.png"), resistivity_log10_range, water),
        plot_mt2d_data_fit(run.observed, run.predicted; output_path = path("DataFit.png")),
    ]
    run.history === nothing || push!(paths, plot_inv2d_convergence(run.history; output_path = path("Convergence.png"),
                                                                    target_rms = run.ctrl.target_rms))
    if hasproperty(run, :vfsa) && run.vfsa !== nothing
        v = run.vfsa
        ens = v.ensemble
        append!(paths, [
            plot_mt2d_model(run.mesh, 10.0 .^ ens.median; output_path = path("ModelMedian.png"), core...),
            plot_mt2d_model(run.mesh, v.best_resistivity; output_path = path("ModelBest.png"), core...),
            plot_mt2d_model(run.mesh, 10.0 .^ ens.p05; output_path = path("ModelP05.png"), core...),
            plot_mt2d_model(run.mesh, 10.0 .^ ens.p95; output_path = path("ModelP95.png"), core...),
            plot_mt2d_model(run.mesh, ifelse.(run.active, ens.std, NaN); output_path = path("ModelStd.png"), core...,
                            log10_values = true, colormap = :viridis, colorbar_label = "σ of log10 ρ",
                            resistivity_log10_range = (0.0, max(0.1, maximum(ens.std[run.active])))),
            plot_mt2d_data_fit(run.observed, _inv2d_predicted(v.best_response, run.observed, run.ctrl.mode);
                               output_path = path("DataFitBest.png")),
            plot_vfsa2d_convergence([c.history for c in v.chains]; output_path = path("ConvergenceVFSA.png"),
                                    target_rms = run.ctrl.target_rms),
        ])
    end
    if true_model_path !== nothing
        truth = ReadModel2D(true_model_path)
        push!(paths, plot_mt2d_model(_mt2d_earth_mesh(truth; receivers = run.mesh.receiver_positions), truth.resistivity;
                                     output_path = path("ModelTrue.png"), core..., water = nothing))
    end
    paths
end

# 2D MT data plots
# Author: @pankajkmishra
# Apparent resistivity and phase curves: observed with error bars, predicted as lines

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

# one rho/phase pair per station, observed with error bars, predicted as lines
function _plot_mt2d_curves(observed::Union{Nothing, DataFile2D}, predicted::Union{Nothing, DataFile2D},
                           stations::AbstractVector{<:Integer}, output_path::AbstractString)
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
        if observed !== nothing
            push!(elements, MarkerElement(; color = _MT_COLOURS[pol], _MT_MARKER...))
            push!(labels, "$pol observed")
        end
        if predicted !== nothing
            push!(elements, LineElement(color = _MT_COLOURS[pol], linewidth = _MT_LINEWIDTH))
            push!(labels, "$pol predicted")
        end
    end
    Legend(figure[0, :], elements, labels; orientation = :horizontal, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

"""
    plot_mt2d_data_fit(observed, predicted; output_path, station_indices=nothing) -> path

Observed apparent resistivity and phase (circles with error bars) against the
predicted curves (lines), TE and TM, one panel pair per station; all stations by default.
"""
function plot_mt2d_data_fit(observed::DataFile2D, predicted::DataFile2D;
                            output_path::AbstractString,
                            station_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing)
    size(observed.z_xy) == size(predicted.z_xy) || throw(DimensionMismatch("observed and predicted surveys differ"))
    _plot_mt2d_curves(observed, predicted, something(station_indices, 1:length(observed.receivers)), output_path)
end

"""
    PlotData2D(data_path; output_path, predicted_path=nothing, station_indices=nothing) -> path

Apparent resistivity and phase of a data file, with error bars, and the curves of
`predicted_path` on top when given.
"""
function PlotData2D(data_path::AbstractString;
                    output_path::AbstractString,
                    predicted_path::Union{Nothing, AbstractString} = nothing,
                    station_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing)
    observed = load_data2d(data_path)
    predicted = predicted_path === nothing ? nothing : load_data2d(predicted_path)
    _plot_mt2d_curves(observed, predicted, something(station_indices, 1:length(observed.receivers)), output_path)
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

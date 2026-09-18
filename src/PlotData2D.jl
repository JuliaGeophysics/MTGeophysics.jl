# 2D MT data plots
# Author: @pankajkmishra
# Response pseudo-sections, station curves, and observed-vs-predicted fit

using CairoMakie

"""
    phase180(values)

Inputs:
- Phase values in degrees.

Output:
- Phase values folded to `[-180, 180]`.

Description:
- Wraps phase angles for the 2D map plot.
"""
phase180(values) = ((values .+ 180) .% 360) .- 180

"""
    plot_mt2d_data_maps(response; output_path)

Inputs:
- 2D MT response and output image path.

Output:
- `String`: Path to the written plot.

Description:
- Writes the standard 2D TE/TM response maps.
"""
function plot_mt2d_data_maps(
    response::MT2DResponse;
    output_path::AbstractString,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    rx_km = response.receivers ./ 1000
    periods = response.periods
    log_periods = log10.(periods)
    tm_phase = _phase_fold_to_0_90(response.phase_yx)
    figure = Figure(size = (1400, 900))

    # Clamp negative apparent resistivities (numerical artifacts) to eps() before log10
    rho_xy_safe = max.(response.rho_xy, eps())
    rho_yx_safe = max.(response.rho_yx, eps())

    ax1 = Axis(figure[1, 1], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TE log10(ρxy)")
    hm1 = heatmap!(ax1, rx_km, log_periods, log10.(rho_xy_safe)', colormap = :Spectral)
    Colorbar(figure[1, 2], hm1)

    ax2 = Axis(figure[1, 3], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TM log10(ρyx)")
    hm2 = heatmap!(ax2, rx_km, log_periods, log10.(rho_yx_safe)', colormap = :Spectral)
    Colorbar(figure[1, 4], hm2)

    ax3 = Axis(figure[2, 1], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TE phase")
    hm3 = heatmap!(ax3, rx_km, log_periods, phase180(response.phase_xy)', colormap = :Spectral, colorrange = (-180, 180))
    Colorbar(figure[2, 2], hm3, ticks = [-180, -90, 0, 90, 180])

    ax4 = Axis(figure[2, 3], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TM phase")
    hm4 = heatmap!(ax4, rx_km, log_periods, tm_phase', colormap = :Spectral, colorrange = (0, 90))
    Colorbar(figure[2, 4], hm4, ticks = [0, 30, 60, 90])

    save(output_path, figure)
    String(output_path)
end

"""
    plot_mt2d_site_curves(response; station_index=cld(length(response.receivers), 2), output_path)

Inputs:
- 2D MT response, station index, and output image path.

Output:
- `String`: Path to the written plot.

Description:
- Writes the representative-station TE/TM response curves.
"""
function plot_mt2d_site_curves(
    response::MT2DResponse;
    station_index::Integer = cld(length(response.receivers), 2),
    output_path::AbstractString,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    tm_phase = _phase_fold_to_0_90(response.phase_yx[:, station_index])
    figure = Figure(size = (900, 700))
    rho_axis = Axis(
        figure[1, 1],
        xlabel = "Period (s)",
        ylabel = "Apparent resistivity (Ω·m)",
        xscale = log10,
        yscale = log10,
        title = "Station $(station_index) response",
    )
    phase_axis = Axis(
        figure[2, 1],
        xlabel = "Period (s)",
        ylabel = "Phase (deg)",
        xscale = log10,
        title = "Phase",
    )

    # Clamp negative apparent resistivities (numerical artifacts) to eps() before log-scale plot
    rho_xy_site = max.(response.rho_xy[:, station_index], eps())
    rho_yx_site = max.(response.rho_yx[:, station_index], eps())
    lines!(rho_axis, response.periods, rho_xy_site, color = :navy, linewidth = 3, label = "TE")
    lines!(rho_axis, response.periods, rho_yx_site, color = :darkorange, linewidth = 3, label = "TM")
    lines!(phase_axis, response.periods, response.phase_xy[:, station_index], color = :navy, linewidth = 3, label = "TE")
    lines!(phase_axis, response.periods, tm_phase, color = :darkorange, linewidth = 3, label = "TM")
    axislegend(rho_axis, position = :rb)
    axislegend(phase_axis, position = :rb)

    save(output_path, figure)
    String(output_path)
end

"""
    PlotData2D(data_path; maps_output_path, curves_output_path, station_index=nothing)

Inputs:
- 2D data path, map image path, curve image path, and optional station index.

Output:
- Named tuple with the written plot paths.

Description:
- Loads a 2D data file and writes both the response maps and a representative-station plot.
"""
function PlotData2D(
    data_path::AbstractString;
    maps_output_path::AbstractString,
    curves_output_path::AbstractString,
    station_index::Union{Nothing, Int} = nothing,
)
    data = load_data2d(data_path)
    response = data_to_response2d(data)
    plot_mt2d_data_maps(response; output_path = maps_output_path)
    plot_mt2d_site_curves(
        response;
        station_index = something(station_index, cld(length(response.receivers), 2)),
        output_path = curves_output_path,
    )
    (maps_output_path = maps_output_path, curves_output_path = curves_output_path)
end

#---------- data fit ----------

"""
    plot_mt2d_data_fit(observed, predicted; output_path, station_indices=nothing)

Inputs:
- `observed`, `predicted`: 2D data objects on the same survey.
- `output_path`: Output image path.
- `station_indices`: Stations to show; `nothing` = first, middle, and last.

Output:
- `String`: Path to the written plot.

Description:
- Observed apparent resistivity and phase (markers, with error bars from the impedance
  errors) against the predicted curves (lines), TE and TM, one column per station.
"""
function plot_mt2d_data_fit(
    observed::DataFile2D,
    predicted::DataFile2D;
    output_path::AbstractString,
    station_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    nr = length(observed.receivers)
    stations = something(station_indices, unique([1, cld(nr, 2), nr]))
    T = observed.periods
    figure = Figure(size = (380 * length(stations), 700))
    for (col, s) in enumerate(stations)
        ax_rho = Axis(figure[1, col], xscale = log10, yscale = log10,
                      ylabel = col == 1 ? "Apparent resistivity (Ω·m)" : "",
                      title = @sprintf("%s (%.1f km)", observed.site_names[s], observed.receivers[s] / 1000))
        ax_pha = Axis(figure[2, col], xscale = log10, xlabel = "Period (s)", ylabel = col == 1 ? "Phase (deg)" : "")
        for (pol, colour) in (("xy", :navy), ("yx", :darkorange))
            z = getproperty(observed, Symbol("z_", pol))[:, s]
            err = getproperty(observed, Symbol("z_", pol, "_error"))[:, s]
            rho_o = max.(getproperty(observed, Symbol("rho_", pol))[:, s], eps())
            rho_p = max.(getproperty(predicted, Symbol("rho_", pol))[:, s], eps())
            pha_o = getproperty(observed, Symbol("phase_", pol))[:, s]
            pha_p = getproperty(predicted, Symbol("phase_", pol))[:, s]
            # TM phase sits in the third quadrant, fold both into (0, 90)
            if pol == "yx"
                pha_o, pha_p = _phase_fold_to_0_90(pha_o), _phase_fold_to_0_90(pha_p)
            end
            # first-order errors: δρ/ρ = 2|δZ|/|Z|, δφ = |δZ|/|Z| rad
            rel = err ./ abs.(z)
            label = pol == "xy" ? "TE" : "TM"
            errorbars!(ax_rho, T, rho_o, rho_o .* min.(2 .* rel, 0.9), rho_o .* 2 .* rel, color = (colour, 0.5))
            scatter!(ax_rho, T, rho_o, color = colour, markersize = 7, label = label)
            lines!(ax_rho, T, rho_p, color = colour, linewidth = 2)
            errorbars!(ax_pha, T, pha_o, rad2deg.(rel), color = (colour, 0.5))
            scatter!(ax_pha, T, pha_o, color = colour, markersize = 7)
            lines!(ax_pha, T, pha_p, color = colour, linewidth = 2)
        end
        col == 1 && axislegend(ax_rho, position = :rb)
    end
    save(output_path, figure)
    String(output_path)
end

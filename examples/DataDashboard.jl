# MT data dashboard: page through the sites of a ModEM data file, observed against an optional ModEM response
# Author: @pankajkmishra
# A directory of EDIs (or one EDI) is converted to <dir>.dat with EDIToModEM first; the window needs GLMakie and a display
# Usage: julia --project=. examples/DataDashboard.jl <data.dat | edi_dir> [predicted.dat] [--png <dir>]
#        --png writes one PNG per site with CairoMakie instead of opening the window
# Residuals and RMS use errors of at least Z_FLOOR·sqrt|Zxy Zyx| and T_FLOOR; PANELS is four of :tipper :ptstrip :beta
# :skew :strike :resid :relerr :nb, or nothing to let the dashboard choose

using MTGeophysics

const Z_FLOOR = 0.05
const T_FLOOR = 0.03
const PANELS = nothing

flags = findall(==("--png"), ARGS)
png_dir = isempty(flags) ? "" : ARGS[flags[1] + 1]
positional = [a for (i, a) in enumerate(ARGS) if !(i in flags) && !(i - 1 in flags)]
isempty(positional) && error("usage: julia --project=. examples/DataDashboard.jl <data.dat | edi_dir> [predicted.dat] [--png <dir>]")

observed = positional[1]
(isdir(observed) || endswith(lowercase(observed), ".edi")) && (observed = EDIToModEM(observed))
predicted = length(positional) >= 2 ? positional[2] : nothing

DataDashboard(observed; predicted, z_floor = Z_FLOOR, t_floor = T_FLOOR, panels = PANELS,
              interactive = isempty(png_dir), snapshot_dir = png_dir)

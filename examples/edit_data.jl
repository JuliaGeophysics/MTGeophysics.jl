# MT data editor: page through the sites of a ModEM data file, mask data and save the mask
# Author: @pankajkmishra
# A directory of EDIs (or one EDI) is converted to <dir>.dat with EDIToModEM first; the window needs GLMakie and a display
# Usage: julia --project=. examples/edit_data.jl <data.dat | edi_dir> [predicted.dat] [--mask <mask.txt>] [--png <dir>]
#        --mask starts from a saved mask, and "Save mask" writes back to it (default <name>_mask.txt here)
#        --png writes one PNG per site with CairoMakie instead of opening the window
# A saved mask applies to the EDIs or to any ModEM file of the survey: apply_data_mask("<edi_dir | data.dat>", "<mask.txt>")
# "Write ModEM" and "Write EDI" write the kept data with their recorded errors (no error floor) beside the data file,
# under a new <name>_<date_time> name each time

using MTGeophysics

option(name) = (i = findfirst(==(name), ARGS); i === nothing ? "" : ARGS[i + 1])
png_dir, mask = option("--png"), option("--mask")
flags = findall(a -> a in ("--png", "--mask"), ARGS)
positional = [a for (i, a) in enumerate(ARGS) if !(i in flags) && !(i - 1 in flags)]
isempty(positional) && error("usage: julia --project=. examples/edit_data.jl <data.dat | edi_dir> [predicted.dat] [--mask <mask.txt>] [--png <dir>]")

observed = positional[1]
if isdir(observed)
    edis = filter(f -> endswith(lowercase(f), ".edi"), readdir(observed))
    println("Found $(length(edis)) EDI files in $observed")
    println("Converting to ModEM: writing $(rstrip(abspath(observed), '/')).dat")
    observed = EDIToModEM(observed)
elseif endswith(lowercase(observed), ".edi")
    println("Found one EDI file, $observed; converting to ModEM")
    observed = EDIToModEM(observed)
else
    println("Found ModEM data file $observed")
end
predicted = length(positional) >= 2 ? positional[2] : nothing

DataDashboard(observed; predicted, mask = isempty(mask) ? nothing : mask,
              interactive = isempty(png_dir), snapshot_dir = png_dir)

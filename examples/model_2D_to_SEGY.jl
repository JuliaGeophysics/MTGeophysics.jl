# 2D resistivity model to SEG-Y, for 3D scenes (TerraScope.jl) and seismic section importers
# Author: @pankajkmishra
# A SEG-Y rev 1 depth section written with SegyIO.jl: log10 resistivity (ohm m) as IEEE float32 samples,
# one trace per uniform step along the profile. Trace X/Y follow the stations' lat/lon, piecewise linear
# between them and straight past the ends, so a curved profile stays curved in a 3D scene
# Usage: julia --project=. examples/model_2D_to_SEGY.jl [model.rho data.dat [model.sgy]]
#        without arguments it converts model.rho of the newest run in examples/data/2D-III

using MTGeophysics
using Printf
import SegyIO

#---------- text header ----------

# printable ascii 32:126 in EBCDIC (code page 037), the rev 1 text header encoding
const EBCDIC = UInt8[
    0x40, 0x5a, 0x7f, 0x7b, 0x5b, 0x6c, 0x50, 0x7d, 0x4d, 0x5d, 0x5c, 0x4e, 0x6b, 0x60, 0x4b, 0x61,
    0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0x7a, 0x5e, 0x4c, 0x7e, 0x6e, 0x6f,
    0x7c, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6,
    0xd7, 0xd8, 0xd9, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xba, 0xe0, 0xbb, 0xb0, 0x6d,
    0x79, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96,
    0x97, 0x98, 0x99, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xc0, 0x4f, 0xd0, 0xa1,
]

# 40 cards of 80 characters, "C 1 ..." to "C40 ..."
function segy_text_header(lines)
    card(i) = rpad(@sprintf("C%2d %s", i, i <= length(lines) ? uppercase(lines[i]) : ""), 80)[1:80]
    String([32 <= Int(c) <= 126 ? EBCDIC[Int(c) - 31] : 0x40 for i in 1:40 for c in card(i)])
end

#---------- profile geometry ----------

# utm zone of a WGS84 point, as an EPSG code
utm_epsg(lat, lon) = "EPSG:" * string((lat >= 0 ? 32600 : 32700) + clamp(floor(Int, (lon + 180) / 6) + 1, 1, 60))

# stations in map coordinates, sorted and unique in profile y; local (y, 0) without lat/lon
function station_map(data::DataFile2D, crs)
    order = sortperm(data.receivers)
    keep = order[[true; diff(data.receivers[order]) .> 0]]
    ys = data.receivers[keep]
    if !(any(!iszero, data.latitudes) && any(!iszero, data.longitudes))
        crs === nothing || throw(ArgumentError("crs given, but the data file has no station lat/lon"))
        return ys, ys, zeros(length(ys)), "LOCAL"
    end
    crs = something(crs, utm_epsg(sum(data.latitudes) / length(data.latitudes), sum(data.longitudes) / length(data.longitudes)))
    to_map = MTGeophysics.Proj.Transformation("EPSG:4326", crs; always_xy = true)
    en = [to_map(data.longitudes[i], data.latitudes[i]) for i in keep]
    ys, first.(en), last.(en), crs
end

# map position of profile coordinate y: piecewise linear between the stations, and past the
# ends along the end segment
function profile_xy(ys, es, ns, y)
    length(ys) == 1 && return es[1] + (y - ys[1]), ns[1]
    if ys[1] <= y <= ys[end]
        k = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
        t = (y - ys[k]) / (ys[k+1] - ys[k])
        return es[k] + t * (es[k+1] - es[k]), ns[k] + t * (ns[k+1] - ns[k])
    end
    a, b, e = y < ys[1] ? (1, 2, 1) : (length(ys) - 1, length(ys), length(ys))
    L = hypot(es[b] - es[a], ns[b] - ns[a])
    s = y - ys[e]
    es[e] + s * (es[b] - es[a]) / L, ns[e] + s * (ns[b] - ns[a]) / L
end

#---------- converter ----------

"""
    model_2D_to_SEGY(model_path, data_path, segy_path=model.sgy next to the model;
                     y_range=:core, dy=nothing, dz=nothing, max_depth=Inf, crs=nothing) -> segy_path

Write a ModEM-layout 2D model as a SEG-Y rev 1 depth section through SegyIO.jl. Samples
are log10 resistivity (ohm m), IEEE float32, taken from the cell under each trace and
sample, so nothing is smoothed. The data file gives the station lat/lon and profile y.
- `y_range`: `:core` (uniform cells, no padding), `:all`, or `(y1, y2)` in metres
- `dy`: trace spacing, default the narrowest cell in `y_range`
- `dz`: sample interval in whole metres, default the thinnest layer; sample 1 is the model top
- `max_depth`: deepest sample, default the model bottom
- `crs`: projected CRS of the trace X/Y, default the UTM zone of the survey

Trace X/Y are in source, group and CDP (scalar -100, metres), the profile y in the shot
point (byte 197, scalar -100) and the sample interval in metres in `dt`. The text header
records the CRS, units and sampling.
"""
function model_2D_to_SEGY(model_path::AbstractString, data_path::AbstractString,
                          segy_path::AbstractString = joinpath(dirname(abspath(model_path)), "model.sgy");
                          y_range = :core, dy = nothing, dz = nothing, max_depth::Real = Inf, crs = nothing)
    model, data = ReadModel2D(model_path), load_data2d(data_path)
    ρ = model.resistivity
    y_nodes = model.origin[2] .+ vcat(0.0, cumsum(model.y_cell_sizes))
    z_nodes = vcat(0.0, cumsum(model.z_cell_sizes))

    #---------- sampling ----------
    y_range isa Symbol && y_range ∉ (:core, :all) && throw(ArgumentError("y_range must be :core, :all or (y1, y2)"))
    cols = y_range == :core ? MTGeophysics._core_range(model.y_cell_sizes) : 1:length(model.y_cell_sizes)
    y1, y2 = y_range isa Symbol ? (y_nodes[first(cols)], y_nodes[last(cols)+1]) : Float64.(y_range)
    y_nodes[1] <= y1 < y2 <= y_nodes[end] || throw(ArgumentError("y_range must lie inside the model"))
    inside = [i for i in eachindex(model.y_cell_sizes) if y_nodes[i+1] > y1 && y_nodes[i] < y2]
    Δy = Float64(something(dy, minimum(model.y_cell_sizes[inside])))
    Δy > 0 || throw(ArgumentError("dy must be positive"))
    depth = min(z_nodes[end], Float64(max_depth))
    depth > 0 || throw(ArgumentError("max_depth must be positive"))
    Δz = something(dz, max(1, round(Int, minimum(model.z_cell_sizes)), ceil(Int, depth / typemax(Int16))))
    isinteger(Δz) && Δz >= 1 || throw(ArgumentError("dz must be a whole number of metres, at least 1"))
    Δz = Int(Δz)
    n_samples = floor(Int, depth / Δz - 1e-9) + 1
    n_samples <= typemax(Int16) || throw(ArgumentError("$n_samples samples per trace, over the SEG-Y limit; raise dz or lower max_depth"))
    n_traces = floor(Int, (y2 - y1) / Δy + 1e-9)
    n_traces >= 2 || throw(ArgumentError("dy leaves fewer than two traces"))
    ytrace = y1 .+ Δy .* ((1:n_traces) .- 0.5)
    iz = [clamp(searchsortedlast(z_nodes, (k - 1) * Δz), 1, size(ρ, 1)) for k in 1:n_samples]
    iy = [clamp(searchsortedlast(y_nodes, y), 1, size(ρ, 2)) for y in ytrace]

    #---------- samples and headers ----------
    ys, es, ns, crs = station_map(data, crs)
    xy = [profile_xy(ys, es, ns, y) for y in ytrace]
    block = SegyIO.SeisBlock(Float32.(log10.(ρ[iz, iy])))
    text = segy_text_header([
        "MTGeophysics.jl 2D resistivity model, $(basename(model_path))",
        "Samples: log10 resistivity (ohm m), IEEE float32, big-endian",
        @sprintf("Depth section: sample 1 at the model top, interval %d m (dt), %d samples", Δz, n_samples),
        @sprintf("Traces: %d, spacing %.6g m along the profile, cdp = trace number", n_traces, Δy),
        "Coordinates: $crs, source = group = cdp x y (m), scalar -100",
        "Profile y (m) of each trace in the shot point, byte 197, scalar -100",
        @sprintf("Profile y from %.6g to %.6g m, %d stations", y1, y2, length(ys)),
        "SEG-Y rev 1",
    ])
    bfh = block.fileheader.bfh
    for (name, value) in ((:Job, 1), (:Line, 1), (:Reel, 1), (:DataTracePerEnsemble, 1), (:dt, Δz), (:dtOrig, Δz),
                          (:nsOrig, n_samples), (:EnsembleFold, 1), (:TraceSorting, 4), (:MeasurementSystem, 1),
                          (:SegyFormatRevisionNumber, 0x0100), (:FixedLengthTraceFlag, 1))
        SegyIO.set_fileheader!(bfh, name, value)
    end
    block.fileheader = SegyIO.FileHeader(text, bfh)
    cm(v) = round.(Int32, 100 .* v)
    trace, unit = collect(1:n_traces), ones(Int, n_traces)
    for (name, value) in ((:TraceNumWithinLine, trace), (:TraceNumWithinFile, trace), (:FieldRecord, unit),
                          (:TraceNumber, trace), (:EnergySourcePoint, unit), (:CDP, trace), (:CDPTrace, unit),
                          (:TraceIDCode, unit), (:DataUse, unit), (:CoordUnits, unit), (:Inline3D, unit),
                          (:Crossline3D, trace), (:dt, fill(Δz, n_traces)), (:ElevationScalar, fill(-100, n_traces)),
                          (:RecSourceScalar, fill(-100, n_traces)), (:ShotPointScalar, fill(-100, n_traces)),
                          (:SourceX, cm(first.(xy))), (:SourceY, cm(last.(xy))), (:GroupX, cm(first.(xy))),
                          (:GroupY, cm(last.(xy))), (:CDPX, cm(first.(xy))), (:CDPY, cm(last.(xy))),
                          (:ShotPoint, cm(ytrace)))
        SegyIO.set_traceheader!(block.traceheaders, name, value)
    end
    mkpath(dirname(abspath(segy_path)))
    SegyIO.segy_write(String(segy_path), block)
    String(segy_path)
end

#---------- script entry ----------

if abspath(PROGRAM_FILE) == @__FILE__
    case_dir = joinpath(@__DIR__, "data", "2D-III")
    args = if isempty(ARGS)
        runs = sort(filter(d -> startswith(d, "run_") && isfile(joinpath(case_dir, d, "model.rho")),
                           isdir(case_dir) ? readdir(case_dir) : String[]))
        isempty(runs) && error("no inversion run in $case_dir; run julia --project=. examples/run_inv2D.jl first")
        [joinpath(case_dir, runs[end], "model.rho"), joinpath(case_dir, "data.dat")]
    elseif length(ARGS) in (2, 3)
        ARGS
    else
        error("usage: julia --project=. examples/model_2D_to_SEGY.jl [model.rho data.dat [model.sgy]]")
    end
    path = model_2D_to_SEGY(args...)
    block = SegyIO.segy_read(path)
    @printf("SEG-Y : %s\n", path)
    @printf("        %d traces x %d samples, dz = %d m, log10 rho %.2f to %.2f\n", size(block.data, 2),
            size(block.data, 1), block.fileheader.bfh.dt, extrema(block.data)...)
end

# Data masks
# Author: @pankajkmishra
# A mask lists the data to leave out, one line per site, period and component, so a mask drawn in DataDashboard on one
# file applies to any copy of the same survey: a ModEM data file (the masked lines are dropped and the rest of the file
# kept as it is) or EDIs (written again without the masked values, together with their ModEM file). Sites match by
# name (case, surrounding blanks and spaces against underscores ignored) and periods within rtol; a period or a
# component of * stands for all of them. Components belong to the frame the mask was drawn in, which the file records
# as "# zrot"; a target in another frame is flagged

const _MASK_COMPONENTS = ("Zxx", "Zxy", "Zyx", "Zyy", "Tzx", "Tzy")
const _MASK_MODEM = Dict("ZXX" => 1, "ZXY" => 2, "ZYX" => 3, "ZYY" => 4, "TX" => 5, "TY" => 6,
                         "RHOXY" => 2, "PHSXY" => 2, "RHOYX" => 3, "PHSYX" => 3)
const _MaskEntry = @NamedTuple{site::String, T::Float64, component::Int}

_mask_site(s::AbstractString) = uppercase(replace(strip(s), r"\s+" => "_"))

function _mask_component(s::AbstractString)
    s == "*" && return 0
    i = findfirst(c -> uppercase(c) == uppercase(s), _MASK_COMPONENTS)
    i === nothing || return i
    i = get(_MASK_MODEM, uppercase(s), 0)
    i > 0 || error("unknown component \"$s\"; use one of $(join(_MASK_COMPONENTS, " ")) or *")
    i
end

# the fields of a "> ..." header line, nothing for any other line
_mask_header(l::AbstractString) = (s = lstrip(l); startswith(s, '>') ? split(s[2:end]) : nothing)

_mask_hit(e::_MaskEntry, T::Real, ic::Integer, rtol::Real) =
    (isnan(e.T) || abs(T / e.T - 1) <= rtol) && (e.component == 0 || e.component == ic)

function _mask_by_site(entries::AbstractVector{_MaskEntry})
    by = Dict{String, Vector{Int}}()
    for (k, e) in enumerate(entries)
        push!(get!(by, e.site, Int[]), k)
    end
    by
end

function _mask_check_frame(mask, zrot::Real, target::AbstractString)
    isfinite(mask.zrot) && abs(mask.zrot - zrot) > 0.01 &&
        @warn @sprintf("the mask was drawn in a frame at %.2f°, %s is at %.2f°; components may not correspond", mask.zrot, target, zrot)
end

"""
    write_data_mask(path, d::Data, keep; source=d.name) -> path

Write the data that `keep` (`nf × 6 × ns`, components Zxx Zxy Zyx Zyy Tzx Tzy) leaves out
as a mask file: one `site period component` line per masked datum, `period *` when every
component at a period is masked and `* *` for a whole site. Entries without data are
not listed.
"""
function write_data_mask(path::AbstractString, d::Data, keep::AbstractArray{Bool, 3}; source::AbstractString = d.name)
    size(keep) == (d.nf, 6, d.ns) || throw(DimensionMismatch("keep must be $(d.nf)×6×$(d.ns), got $(size(keep))"))
    present = cat(isfinite.(d.Z), isfinite.(d.tip); dims = 2)
    kept = present .& keep
    masked = present .& .!keep
    lines = 0
    open(path, "w") do io
        println(io, "# MTGeophysics data mask, ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), isempty(source) ? "" : ", from $source")
        println(io, "# one masked datum per line; * is every period or every component; components ", join(_MASK_COMPONENTS, " "))
        isempty(d.zrot) || !allequal(d.zrot) || @printf(io, "# zrot %.2f\n", d.zrot[1])
        println(io, "# site  period_s  component")
        for is in 1:d.ns
            any(masked[:, :, is]) || continue
            if !any(kept[:, :, is])
                println(io, d.site[is], "  *  *"); lines += 1
                continue
            end
            for ip in 1:d.nf
                any(masked[ip, :, is]) || continue
                if !any(kept[ip, :, is])
                    @printf(io, "%s  %.6e  *\n", d.site[is], d.T[ip]); lines += 1
                else
                    for ic in findall(masked[ip, :, is])
                        @printf(io, "%s  %.6e  %s\n", d.site[is], d.T[ip], _MASK_COMPONENTS[ic]); lines += 1
                    end
                end
            end
        end
    end
    @printf("write_data_mask: %d masked data at %d sites in %d lines -> %s\n", count(masked),
            count(is -> any(masked[:, :, is]), 1:d.ns), lines, path)
    String(path)
end

"""
    read_data_mask(path) -> (entries, zrot)

Read a mask file written by `write_data_mask` or by hand: `site period component` per
line, `*` for every period or component, `#` starts a comment. Components are Zxx Zxy
Zyx Zyy Tzx Tzy (ModEM's ZXY, TX, ... also work). `zrot` is the frame the mask was
drawn in, NaN when the file does not say.
"""
function read_data_mask(path::AbstractString)
    entries, zrot = _MaskEntry[], NaN
    for (k, line) in enumerate(eachline(path))
        s = strip(line)
        m = match(r"^#\s*zrot\s+(\S+)", s)
        m === nothing || (zrot = parse(Float64, m[1]); continue)
        (isempty(s) || startswith(s, '#')) && continue
        f = split(s)
        length(f) == 3 || error("$path, line $k: expected \"site period component\", got \"$s\"")
        T = f[2] == "*" ? NaN : tryparse(Float64, f[2])
        (T === nothing || T <= 0) && error("$path, line $k: period \"$(f[2])\" is not a positive number or *")
        push!(entries, (site = _mask_site(f[1]), T = Float64(T), component = _mask_component(f[3])))
    end
    (entries = entries, zrot = zrot)
end

_mask_read(mask) = mask isa AbstractString ? read_data_mask(mask) : mask

"""
    mask_keep(d::Data, mask; rtol=0.01) -> BitArray

The `nf × 6 × ns` array of the data `mask` (a path or the result of `read_data_mask`)
keeps: true where `d` holds a datum and no mask line covers it.
"""
function mask_keep(d::Data, mask; rtol::Real = 0.01)
    m = _mask_read(mask)
    keep = cat(isfinite.(d.Z), isfinite.(d.tip); dims = 2)
    by, used = _mask_by_site(m.entries), falses(length(m.entries))
    for is in 1:d.ns
        ks = get(by, _mask_site(d.site[is]), nothing)
        ks === nothing && continue
        for ip in 1:d.nf, ic in 1:6, k in ks
            _mask_hit(m.entries[k], d.T[ip], ic, rtol) || continue
            used[k] = true
            keep[ip, ic, is] = false
        end
    end
    any(!, used) && @warn "$(count(!, used)) mask line(s) match no site and period: $(join(unique(e.site for e in m.entries[.!used]), ", "))"
    keep
end

"""
    apply_data_mask!(d::Data, mask; rtol=0.01) -> Int

Set the data `mask` covers (and their errors) to NaN in `d` and return how many were
set. Periods match within `rtol`.
"""
function apply_data_mask!(d::Data, mask; rtol::Real = 0.01)
    keep = mask_keep(d, mask; rtol)
    n = 0
    for is in 1:d.ns, ip in 1:d.nf, ic in 1:6
        keep[ip, ic, is] && continue
        if ic <= 4
            isfinite(d.Z[ip, ic, is]) || continue
            d.Z[ip, ic, is] = d.Zerr[ip, ic, is] = complex(NaN, NaN)
        else
            isfinite(d.tip[ip, ic - 4, is]) || continue
            d.tip[ip, ic - 4, is] = d.tiperr[ip, ic - 4, is] = complex(NaN, NaN)
        end
        n += 1
    end
    d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)
    n
end

# drop the masked lines of a ModEM file and rewrite each block's "> nperiods nsites" line; blocks left empty go
function _mask_modem_file(src::AbstractString, dst::AbstractString, m; rtol::Real)
    blocks = Tuple{Vector{String}, Vector{String}}[]
    for line in eachline(src)
        s = lstrip(line)
        if isempty(s) || startswith(s, '#') || startswith(s, '>')
            (isempty(blocks) || !isempty(blocks[end][2])) && push!(blocks, (String[], String[]))
            push!(blocks[end][1], line)
        else
            isempty(blocks) && push!(blocks, (String[], String[]))
            push!(blocks[end][2], line)
        end
    end
    by = _mask_by_site(m.entries)
    frames, dropped = Set{Float64}(), 0
    open(dst, "w") do io
        for (header, data) in blocks
            if isempty(data)
                foreach(l -> println(io, l), header)
                continue
            end
            kept = filter(data) do line
                f = split(line)
                length(f) >= 8 || return true
                ks = get(by, _mask_site(f[2]), nothing)
                ic = get(_MASK_MODEM, uppercase(f[8]), 0)
                T = tryparse(Float64, f[1])
                (ks === nothing || ic == 0 || T === nothing) && return true
                !any(k -> _mask_hit(m.entries[k], T, ic, rtol), ks)
            end
            dropped += length(data) - length(kept)
            isempty(kept) && continue
            fields = _mask_header.(header)
            counts = findlast(f -> f !== nothing && length(f) == 2 && all(x -> tryparse(Int, x) !== nothing, f), fields)
            for (i, l) in enumerate(header)
                f = fields[i]
                f !== nothing && length(f) == 1 && (r = tryparse(Float64, f[1])) !== nothing && push!(frames, r)
                if i == counts
                    println(io, "> ", length(unique(split(k)[1] for k in kept)), " ", length(unique(split(k)[2] for k in kept)))
                else
                    println(io, l)
                end
            end
            foreach(l -> println(io, l), kept)
        end
    end
    foreach(r -> _mask_check_frame(m, r, src), frames)
    dropped
end

"""
    apply_data_mask(target, mask; output="", rtol=0.01) -> path or (edi_dir, modem)

Apply a mask file to a copy of the survey without touching the original.

  - `target` a ModEM data file: the masked lines are dropped, everything else is kept
    as written, and the result goes to `output` (default `<name>_masked_<date_time>.dat`
    next to it).
  - `target` a directory of EDIs (or one EDI): the EDIs are read as `EDIToModEM` reads
    them (sites in different frames turned to north), the masked values removed, and
    each site with data left written to `<dir>_masked_<date_time>/` (or `output`), with
    the ModEM file of them as `<that directory>.dat`.

Sites match by name, periods within `rtol`.
"""
function apply_data_mask(target::AbstractString, mask; output::AbstractString = "", rtol::Real = 0.01)
    m = _mask_read(mask)
    stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    src = rstrip(abspath(target), '/')
    if isfile(src) && !occursin(r"\.edi$"i, src)
        dst = isempty(output) ? splitext(src)[1] * "_masked_$stamp.dat" : String(output)
        n = _mask_modem_file(src, dst, m; rtol)
        println("apply_data_mask: $n data lines dropped -> $dst")
        return dst
    end
    (isdir(src) || isfile(src)) || error("not found: $target")
    d = load_data_edi(src)
    _edi_to_one_frame!(d)
    _mask_check_frame(m, d.zrot[1], target)
    d.site = replace.(strip.(d.site), r"\s+" => "_")
    n = apply_data_mask!(d, m; rtol)
    dir = isempty(output) ? (isdir(src) ? src : splitext(src)[1]) * "_masked_$stamp" : rstrip(String(output), '/')
    mkpath(dir)
    left = [is for is in 1:d.ns if any(isfinite, d.Z[:, :, is]) || any(isfinite, d.tip[:, :, is])]
    isempty(left) && error("every datum is masked")
    for is in left
        _write_edi(joinpath(dir, d.site[is] * ".edi"), d, is, basename(src))
    end
    println("apply_data_mask: $n data masked, $(length(left)) of $(d.ns) EDIs -> $dir")
    modem = EDIToModEM(dir, dir * ".dat")
    (edi_dir = dir, modem = modem)
end

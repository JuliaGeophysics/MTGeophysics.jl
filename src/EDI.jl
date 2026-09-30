# EDI reader and writer
# Author: @pankajkmishra
# EDIToModEM reads SEG EDI files, one or a directory on a shared period axis, into a ModEM data file; ModEMToEDI writes
# one EDI per site back. Three layouts are read: transfer functions (>ZXXR ... >TYVAR.EXP, errors the square root of
# the variances), apparent resistivity and phase only (>RHOXY/>PHSXY, |Z| from ρ and δ|Z|/|Z| = δρ/2ρ), and cross
# spectra (>=SPECTRASECT). A spectra block is the real matrix of the listed channels: autopowers on the diagonal, the
# real part of <A B*> below it and minus its imaginary part above. Z = <E R*><H R*>⁻¹ and T = <Hz R*><H R*>⁻¹, R the
# two channels after the local five (the local H when there are none); errors are the 68 % bounds of the local
# regression, 4·F(4, AVGT-4)/(AVGT-4) times the residual power times diag<H H*>⁻¹. EMPTY values, EMPTY cross powers
# and all-zero placeholders (1e-10 in the MTNet files) are NaN. Impedances in [mV/km]/[nT] go to Ohm, exp(+iωt), as in
# load_data_modem; a file one unit factor off (SI impedances, or the MTNet BC87 double conversion) is caught by its
# apparent resistivities. A DATAID with spaces, or shared by several files, gives way to SECTID, then to the file
# name. Station x, y are the local transverse Mercator of MakeMesh3D about the sites' centroid, and sites recorded in
# different ZROT frames are turned to north, since a ModEM file holds one rotation.
# References: SEG MT/EMAP data interchange standard (1987); Bendat & Piersol, Random Data (2010)

const _EDI_ZCOMP = ("ZXX", "ZXY", "ZYX", "ZYY")
const _EDI_TCOMP = ("TX", "TY")
const _EDI_PLACEHOLDER = 1e-9
const _EDI_FIELD = 4π * 1e-7 * 1000

function _edi_sections(path::AbstractString)
    sections = Vector{NamedTuple{(:name, :header, :lines), Tuple{String, String, Vector{String}}}}()
    for raw in eachline(path)
        line = strip(replace(raw, '\r' => ""))
        if startswith(line, ">")
            startswith(line, ">!") && continue
            name = uppercase(first(split(line[2:end] * " ")))
            push!(sections, (name = name, header = String(line), lines = String[]))
        elseif !isempty(sections) && !isempty(line)
            push!(sections[end].lines, String(line))
        end
    end
    sections
end

_edi_numbers(lines) = Float64[x for l in lines for t in split(l) for x in (tryparse(Float64, replace(t, r"[dD]" => "E")),) if x !== nothing]

function _edi_keys(sections)
    keys = Dict{String, String}()
    for s in sections, l in (s.name in ("HEAD", "INFO", "=DEFINEMEAS", "=MTSECT", "=SPECTRASECT") ? s.lines : String[])
        m = match(r"^\s*([A-Za-z][A-Za-z0-9_.]*)\s*[=:]\s*(.*?)\s*$", l)
        m === nothing && continue
        k = uppercase(m.captures[1])
        haskey(keys, k) || (keys[k] = strip(m.captures[2], ['"', ' ']))
    end
    keys
end

function _edi_degrees(s::AbstractString)
    t = strip(s)
    occursin(':', t) || return something(tryparse(Float64, t), NaN)
    parts = something.(tryparse.(Float64, strip.(split(t, ':'))), NaN)
    v = abs(parts[1]) + sum(parts[k] / 60.0^(k - 1) for k in 2:length(parts); init = 0.0)
    startswith(t, "-") ? -v : v
end

_edi_first(keys, names...) = (for n in names; haskey(keys, n) && !isempty(keys[n]) && return keys[n]; end; "")

function _f4_quantile(p::Real, m::Real)
    m > 0 || return NaN
    b = m / 2
    lo, hi = 0.0, 1.0
    for _ in 1:100
        x = (lo + hi) / 2
        1 - (1 - x)^b * (1 + b * x) < p ? (lo = x) : (hi = x)
    end
    x = (lo + hi) / 2
    m * x / (4 * (1 - x))
end

function _edi_from_spectra(sections, path, empty)
    chtype = Dict{Float64, String}()
    for s in sections
        s.name in ("HMEAS", "EMEAS") || continue
        id, ty = match(r"\bID=\s*([0-9.]+)"i, s.header), match(r"\bCHTYPE=\s*(\w+)"i, s.header)
        id === nothing || ty === nothing || (chtype[parse(Float64, id[1])] = uppercase(ty[1]))
    end
    sect = sections[findfirst(s -> s.name == "=SPECTRASECT", sections)]
    k = findfirst(l -> startswith(l, "//"), sect.lines)
    k === nothing && error("no channel list under >=SPECTRASECT in $path")
    ids = _edi_numbers(vcat(replace(sect.lines[k], r"^//\s*\d*" => ""), sect.lines[k+1:end]))
    n = length(ids)
    types = [get(chtype, id, "") for id in ids]
    at(t) = findfirst(==(t), types)
    hx, hy, hz, ex, ey = at("HX"), at("HY"), at("HZ"), at("EX"), at("EY")
    any(isnothing, (hx, hy, ex, ey)) && error("the spectra of $path lack one of HX, HY, EX, EY")
    rest = setdiff(1:n, filter(!isnothing, [hx, hy, hz, ex, ey]))
    rx, ry = length(rest) >= 2 ? (rest[1], rest[2]) : (hx, hy)

    blocks = [s for s in sections if s.name == "SPECTRA"]
    nf = length(blocks)
    nf == 0 && error("no >SPECTRA blocks in $path")
    freq, zrot = fill(NaN, nf), zeros(nf)
    Z, Zerr = fill(complex(NaN, NaN), nf, 4), fill(NaN, nf, 4)
    T, Terr = fill(complex(NaN, NaN), nf, 2), fill(NaN, nf, 2)
    for (i, s) in enumerate(blocks)
        value(key, default) = (m = match(Regex("\\b$key=\\s*([-+0-9.eEdD]+)", "i"), s.header);
                               m === nothing ? default : something(tryparse(Float64, replace(m[1], r"[dD]" => "E")), default))
        freq[i], zrot[i], ν = value("FREQ", NaN), value("ROTSPEC", 0.0), value("AVGT", NaN)
        M = _edi_numbers(s.lines)
        length(M) >= n^2 || continue
        M = permutedims(reshape([abs(x) >= 0.1 * abs(empty) ? NaN : x for x in M[1:n^2]], n, n))
        S = [a == b ? complex(M[a, a]) : a < b ? complex(M[b, a], -M[a, b]) : complex(M[a, b], M[b, a]) for a in 1:n, b in 1:n]
        H, R = [hx, hy], [rx, ry]
        SHR = S[H, R]
        abs(det(SHR)) > 0 || continue
        Z[i, :] .= vec(permutedims(S[[ex, ey], R] / SHR))
        hz === nothing || (T[i, :] .= vec(S[[hz], R] / SHR))
        SHH = S[H, H]
        dH = real(det(SHH))
        (dH > 0 && ν > 4) || continue
        w = real.([SHH[2, 2], SHH[1, 1]]) ./ dH
        k = 4 * _f4_quantile(0.68, ν - 4) / (ν - 4)
        function σ(o)
            r = real(S[o, o] - (S[[o], H] / SHH * S[H, [o]])[1])
            [v >= 0 ? sqrt(v) : NaN for v in k * r .* w]
        end
        Zerr[i, 1:2] .= σ(ex)
        Zerr[i, 3:4] .= σ(ey)
        hz === nothing || (Terr[i, :] .= σ(hz))
    end
    (freq = freq, Z = Z, Zerr = Zerr, tip = T, tiperr = Terr, zrot = zrot)
end

"""
    read_edi(path; fix_units=true) -> NamedTuple

One EDI file: site, lat/lon/elevation, periods, impedance and tipper with their
standard errors, in Ohm and exp(+iωt), and the rotation angles. Missing entries are NaN.

Three layouts are read: transfer functions (>ZXXR ... >TYVAR.EXP, errors the square
root of the variances), apparent resistivity and phase only (>RHOXY, >PHSXY), and
cross spectra (>=SPECTRASECT, Z and T through the remote reference channels, errors
from the local coherence). With `fix_units`, a file whose apparent resistivities land
beyond 10⁶ Ω·m, or below 10⁻³ Ω·m, and come back in range by one [mV/km]/[nT] ↔ Ohm
factor is rescaled with a warning: the MTNet BC87 EDIs, or impedances written in SI.
"""
function read_edi(path::AbstractString; fix_units::Bool = true)
    sections = _edi_sections(path)
    names = Set(s.name for s in sections)
    keys = _edi_keys(sections)
    empty = something(tryparse(Float64, _edi_first(keys, "EMPTY")), 1e32)
    block(name) = (i = findfirst(s -> s.name == name, sections); i === nothing ? nothing : _edi_numbers(sections[i].lines))
    anyblock(names...) = (for n in names; b = block(n); b === nothing || return b; end; nothing)

    if "=SPECTRASECT" in names
        sp = _edi_from_spectra(sections, path, empty)
        freq, Z, Zerr, T, Terr, zrot = sp.freq, sp.Z, sp.Zerr, sp.tip, sp.tiperr, sp.zrot
        trot = copy(zrot)
        nf = length(freq)
    else
        freq = block("FREQ")
        freq === nothing && error("no >FREQ block (nor spectra) in $path")
        nf = length(freq)
        clean(v) = v === nothing ? fill(NaN, nf) :
                   [k <= length(v) && isfinite(v[k]) && abs(v[k]) < 0.1 * abs(empty) ? v[k] : NaN for k in 1:nf]

        Z, Zerr = fill(complex(NaN, NaN), nf, 4), fill(NaN, nf, 4)
        for (ic, c) in enumerate(_EDI_ZCOMP)
            Z[:, ic] .= complex.(clean(block(c * "R")), clean(block(c * "I")))
            Zerr[:, ic] .= sqrt.(abs.(clean(block(c * ".VAR"))))
        end
        zrot = clean(block("ZROT"))
        if all(isnan, real.(Z)) && "RHOXY" in names
            for (ic, c) in enumerate(("XX", "XY", "YX", "YY"))
                ρ, φ, δρ = clean(block("RHO" * c)), clean(block("PHS" * c)), clean(block("RHO" * c * ".ERR"))
                absZ = sqrt.(ρ .* 2π .* freq .* 4π * 1e-7) ./ _EDI_FIELD
                Z[:, ic] .= absZ .* cis.(deg2rad.(φ))
                Zerr[:, ic] .= absZ .* δρ ./ (2 .* ρ)
            end
            zrot = clean(block("RHOROT"))
        end
        all(isnan, zrot) && (zrot .= 0.0)
        T, Terr = fill(complex(NaN, NaN), nf, 2), fill(NaN, nf, 2)
        for (ic, c) in enumerate(_EDI_TCOMP)
            T[:, ic] .= complex.(clean(anyblock(c * "R.EXP", c * "R")), clean(anyblock(c * "I.EXP", c * "I")))
            Terr[:, ic] .= sqrt.(abs.(clean(anyblock(c * "VAR.EXP", c * "VAR", c * ".VAR", c * "RVAR.EXP"))))
        end
        trot = clean(anyblock("TROT.EXP", "TROT"))
        all(isnan, trot) && (trot .= zrot)
    end

    placeholder(z) = abs(real(z)) <= _EDI_PLACEHOLDER && abs(imag(z)) <= _EDI_PLACEHOLDER
    for A in (Z, T), i in eachindex(A)
        placeholder(A[i]) && (A[i] = complex(NaN, NaN))
    end
    Zerr[isnan.(real.(Z))] .= NaN
    Terr[isnan.(real.(T))] .= NaN

    if occursin(r"exp\(\s*-", lowercase(_edi_first(keys, "SIGNCONVENTION")))
        Z .= conj.(Z); T .= conj.(T)
    end

    scale = _EDI_FIELD
    if fix_units
        ρ = [abs2(Z[i, ic] * scale) / (2π * freq[i] * 4π * 1e-7) for i in 1:nf, ic in 2:3 if isfinite(Z[i, ic])]
        m = isempty(ρ) ? NaN : median(ρ)
        if m > 1e6 && m * _EDI_FIELD^2 < 1e6
            @warn "$(basename(path)): apparent resistivities near $(round(m, sigdigits = 2)) Ω·m; impedances rescaled by μ0·1000 once more (as the MTNet BC87 EDIs need)"
            scale *= _EDI_FIELD
        elseif m < 1e-3 && m / _EDI_FIELD^2 > 1e-3
            @warn "$(basename(path)): apparent resistivities near $(round(m, sigdigits = 2)) Ω·m; impedances read as Ohm (SI), not [mV/km]/[nT]"
            scale = 1.0
        end
    end
    Z .*= scale
    Zerr .*= scale

    lat = _edi_degrees(_edi_first(keys, "LAT", "LATITUDE", "REFLAT"))
    lon = _edi_degrees(_edi_first(keys, "LONG", "LON", "LONGITUDE", "REFLONG"))
    elev = something(tryparse(Float64, _edi_first(keys, "ELEV", "ELEVATION")), 0.0)
    elev == 0 && (elev = something(tryparse(Float64, _edi_first(keys, "REFELEV")), 0.0))
    stem = splitext(basename(path))[1]
    site = _edi_first(keys, "DATAID", "SITE", "STATION", "SECTID")
    isempty(site) && (site = stem)
    sectid = strip(_edi_first(keys, "SECTID"))
    occursin(' ', strip(site)) && !isempty(sectid) && !occursin(' ', sectid) && (site = sectid)

    ok = findall(f -> isfinite(f) && f > 0, freq)
    order = ok[sortperm(1 ./ freq[ok])]
    (site = String(site), sectid = String(sectid), stem = String(stem), lat = lat, lon = lon, elev = elev,
     T = (1 ./ freq)[order], Z = Z[order, :], Zerr = Zerr[order, :], tip = T[order, :], tiperr = Terr[order, :],
     zrot = zrot[order], trot = trot[order], path = String(path))
end

function _shared_periods(sets::AbstractVector; rtol::Real = 2e-3)
    axis = Float64[]
    index = [zeros(Int, length(s)) for s in sets]
    members = Set{Int}()
    for (t, is, k) in sort([(t, is, k) for (is, s) in enumerate(sets) for (k, t) in enumerate(s)])
        if isempty(axis) || t > axis[end] * (1 + rtol) || is in members
            push!(axis, t)
            empty!(members)
        end
        push!(members, is)
        index[is][k] = length(axis)
    end
    axis, index
end

_period_index(axis, T; rtol = 2e-3) = (i = argmin(abs.(log.(axis ./ T))); abs(axis[i] / T - 1) <= rtol ? i : nothing)

"""
    load_data_edi(path; fix_units=true, pattern=r"\\.edi\$"i) -> Data

Read one EDI file or every EDI in a directory into a `Data` object, on the union of
the sites' periods: periods of different sites within 0.2 % share an index, those of one
site never do. Each file is read by `read_edi` (transfer
functions, ρ/φ only or cross spectra); in a directory, a file that cannot be read or
holds no periods is skipped with a warning. Impedances come out in Ohm and exp(+iωt),
and `zrot` keeps each site's ZROT. Station x (north) and y (east) are metres in the
transverse Mercator about the sites' centroid, which is also the origin. Site names
come from DATAID; a DATAID shared by several files falls back to SECTID, then to the
file name.
"""
function load_data_edi(path::AbstractString; fix_units::Bool = true, pattern::Regex = r"\.edi$"i)
    files = isdir(path) ? sort(filter(f -> occursin(pattern, f), readdir(path; join = true))) : [String(path)]
    isempty(files) && error("no EDI files in $path")
    edis = Any[]
    for f in files
        e = try
            read_edi(f; fix_units)
        catch err
            length(files) == 1 && rethrow()
            @warn "skipping $(basename(f)): $(sprint(showerror, err))"
            continue
        end
        isempty(e.T) ? @warn("skipping $(basename(f)): no periods") : push!(edis, e)
    end
    isempty(edis) && error("no readable EDI files in $path")
    names = [e.site for e in edis]
    for (i, e) in enumerate(edis)
        count(==(e.site), names) > 1 || continue
        names[i] = !isempty(e.sectid) && count(x -> x.sectid == e.sectid, edis) == 1 ? e.sectid : e.stem
    end

    T, index = _shared_periods([e.T for e in edis])
    nf, ns = length(T), length(edis)
    d = make_nan_data()
    d.name, d.niter = String(path), ""
    d.T, d.f, d.nf, d.ns = T, 1 ./ T, nf, ns
    d.site = names
    d.Z, d.Zerr = fill(complex(NaN, NaN), nf, 4, ns), fill(complex(NaN, NaN), nf, 4, ns)
    d.tip, d.tiperr = fill(complex(NaN, NaN), nf, 2, ns), fill(complex(NaN, NaN), nf, 2, ns)
    d.zrot, d.trot = zeros(nf, ns), zeros(nf, ns)
    for (is, e) in enumerate(edis)
        d.zrot[:, is] .= e.zrot[1]; d.trot[:, is] .= e.trot[1]
        for (k, ip) in enumerate(index[is])
            d.Z[ip, :, is] .= e.Z[k, :]
            d.Zerr[ip, :, is] .= complex.(e.Zerr[k, :], 0.0)
            d.tip[ip, :, is] .= e.tip[k, :]
            d.tiperr[ip, :, is] .= complex.(e.tiperr[k, :], 0.0)
            d.zrot[ip, is], d.trot[ip, is] = e.zrot[k], e.trot[k]
        end
    end
    present = [any(isfinite, real.(d.Z[:, ic, :])) for ic in 1:4]
    tipper = [any(isfinite, real.(d.tip[:, ic, :])) for ic in 1:2]
    d.responses = vcat(collect(_EDI_ZCOMP[present]), collect(_EDI_TCOMP[tipper]))
    d.nr = length(d.responses)

    d.loc = hcat([e.lat for e in edis], [e.lon for e in edis], [e.elev for e in edis])
    ok = isfinite.(d.loc[:, 1]) .& isfinite.(d.loc[:, 2])
    lat0, lon0 = any(ok) ? (mean(d.loc[ok, 1]), mean(d.loc[ok, 2])) : (0.0, 0.0)
    d.origin = [lat0, lon0, 0.0]
    tm = Proj.Transformation("EPSG:4326", _local_tm_proj_string(lat0, lon0); always_xy = true)
    d.x, d.y = fill(NaN, ns), fill(NaN, ns)
    for is in findall(ok)
        e, n = tm((d.loc[is, 2], d.loc[is, 1]))
        d.x[is], d.y[is] = n, e - 500000.0
    end
    d.z = copy(d.loc[:, 3])
    d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)
    any(!, ok) && @warn "$(count(!, ok)) site(s) without coordinates: $(join(names[.!ok], ", "))"
    d
end

function _error_floors(d::Data, z_floor::Real, t_floor::Real)
    ez, et = abs.(d.Zerr), abs.(d.tiperr)
    for is in 1:d.ns, ip in 1:d.nf
        g = z_floor * sqrt(abs(d.Z[ip, 2, is] * d.Z[ip, 3, is]))
        for ic in 1:4
            e = ez[ip, ic, is]
            ez[ip, ic, is] = isfinite(e) ? (isfinite(g) ? max(e, g) : e) : (isfinite(g) && g > 0 ? g : NaN)
        end
        for ic in 1:2
            e = et[ip, ic, is]
            et[ip, ic, is] = isfinite(e) ? max(e, t_floor) : (t_floor > 0 ? t_floor : NaN)
        end
    end
    ez, et
end

function _edi_to_one_frame!(d::Data)
    all(==(d.zrot[1]), d.zrot) && all(==(d.zrot[1]), d.trot) && return d
    dropped = 0
    for is in 1:d.ns, ip in 1:d.nf
        θ = -d.zrot[ip, is]
        if θ != 0 && any(isfinite, d.Z[ip, :, is])
            r = _rotate_impedance(permutedims(reshape(d.Z[ip, :, is], 2, 2)), permutedims(reshape(abs.(d.Zerr[ip, :, is]), 2, 2)), θ)
            if r === nothing
                d.Z[ip, :, is] .= complex(NaN, NaN); dropped += 1
            else
                d.Z[ip, :, is] .= vec(permutedims(r.Z)); d.Zerr[ip, :, is] .= vec(permutedims(r.σ))
            end
        end
        θ = -d.trot[ip, is]
        if θ != 0 && any(isfinite, d.tip[ip, :, is])
            r = _rotate_tipper(d.tip[ip, :, is], abs.(d.tiperr[ip, :, is]), θ)
            if r === nothing
                d.tip[ip, :, is] .= complex(NaN, NaN); dropped += 1
            else
                d.tip[ip, :, is] .= r.T; d.tiperr[ip, :, is] .= r.σ
            end
        end
    end
    @warn "the EDIs are in different frames (ZROT/TROT); every site turned to north" * (dropped > 0 ? ", $dropped incomplete site-period block(s) dropped" : "")
    d.zrot .= 0.0; d.trot .= 0.0
    d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)
    d
end

"""
    EDIToModEM(edi_path, output_path=""; z_floor=0.0, t_floor=0.0, fix_units=true) -> path

Convert one EDI file or a directory of them to a ModEM data file (Full_Impedance and,
where present, Full_Vertical_Components; [mV/km]/[nT], exp(+iωt)). The default output
is `<directory>.dat` next to the directory, or `<file>.dat` next to the file.

The EDIs may be transfer functions, ρ/φ only or cross spectra (see `read_edi`); a file
in a directory that cannot be read is skipped with a warning. Periods are shared
across sites within 0.2 % (two periods of one site stay apart), station x/y are metres in the transverse Mercator about the
sites' centroid (the file's origin), and site names lose any spaces. Sites in
different frames (ZROT/TROT) are turned to north, since ModEM holds one rotation per
file. Errors are the EDI ones; `z_floor` (times √|Zxy·Zyx|) and `t_floor` raise them,
and an entry without an error is written with ModEM's 1e12 unless a floor fills it.
"""
function EDIToModEM(edi_path::AbstractString, output_path::AbstractString = "";
                    z_floor::Real = 0.0, t_floor::Real = 0.0, fix_units::Bool = true)
    d = load_data_edi(edi_path; fix_units)
    _edi_to_one_frame!(d)
    ez, et = _error_floors(d, z_floor, t_floor)
    d.Zerr, d.tiperr = complex.(ez, 0.0), complex.(et, 0.0)
    d.site = replace.(strip.(d.site), r"\s+" => "_")
    d.name = ""
    src = rstrip(abspath(edi_path), '/')
    path = isempty(output_path) ? (isdir(src) ? src : splitext(src)[1]) * ".dat" : String(output_path)
    missing_err = count(i -> isfinite(d.Z[i]) && !isfinite(ez[i]), eachindex(d.Z)) +
                  count(i -> isfinite(d.tip[i]) && !isfinite(et[i]), eachindex(d.tip))
    redirect_stdout(devnull) do
        write_data_modem(path, d; sign = 1, units = "[mV/km]/[nT]", include_tipper = has_tipper_data(d),
                         description = "EDIToModEM from $(basename(src))" *
                                       (z_floor > 0 || t_floor > 0 ? @sprintf(", error floors %.3g sqrt|ZxyZyx| and %.3g", z_floor, t_floor) : ""))
    end
    @printf("EDIToModEM: %d sites, %d periods (%.4g .. %.4g s), %s -> %s\n", d.ns, d.nf, minimum(d.T), maximum(d.T),
            join(d.responses, " "), path)
    missing_err > 0 && @warn "$missing_err data without an error written with 1e12; set z_floor/t_floor to fill them"
    path
end

function _edi_sexagesimal(v::Real)
    isfinite(v) || return "0:00:00.00"
    a = abs(v)
    d = floor(Int, a); m = floor(Int, (a - d) * 60); s = (a - d - m / 60) * 3600
    s >= 59.995 && (s = 0.0; m += 1)
    m == 60 && (m = 0; d += 1)
    @sprintf("%s%d:%02d:%05.2f", v < 0 ? "-" : "", d, m, s)
end

function _edi_write_block(io, header::AbstractString, v::AbstractVector)
    println(io, header, " // ", length(v))
    for chunk in Iterators.partition(v, 6)
        println(io, join((@sprintf("%15.6E", x) for x in chunk), ""))
    end
    println(io)
end

function _write_edi(path::AbstractString, d::Data, is::Integer, source::AbstractString)
    keep = [ip for ip in d.nf:-1:1 if any(isfinite, d.Z[ip, :, is]) || any(isfinite, d.tip[ip, :, is])]
    n = length(keep)
    empty = 1.0e32
    val(x) = isfinite(x) ? x : empty
    lat, lon, elev = d.loc[is, 1], d.loc[is, 2], d.loc[is, 3]
    rot = isempty(d.zrot) ? zeros(n) : d.zrot[keep, is]
    trot = isempty(d.trot) ? rot : d.trot[keep, is]
    tipper = any(isfinite, d.tip[keep, :, is])
    open(path, "w") do io
        println(io, ">HEAD\n  DATAID=\"", d.site[is], "\"\n  ACQBY=\"\"\n  FILEBY=\"MTGeophysics.jl\"")
        println(io, "  FILEDATE=", Dates.format(Dates.today(), "yyyy/mm/dd"))
        println(io, "  LAT=", _edi_sexagesimal(lat), "\n  LONG=", _edi_sexagesimal(lon), "\n  ELEV=", @sprintf("%.2f", elev))
        println(io, "  STDVERS=\"SEG 1.0\"\n  PROGVERS=\"MTGeophysics.jl ModEMToEDI\"\n  EMPTY=1.0E+32\n")
        println(io, ">INFO\n  MAXINFO=999\n  SOURCE=", source, "\n  SIGNCONVENTION=exp(+i\\omega t)\n  UNITS=[mV/km]/[nT]\n")
        println(io, ">=DEFINEMEAS\n  MAXCHAN=", tipper ? 5 : 4, "\n  MAXRUN=999\n  MAXMEAS=9999\n  UNITS=M\n  REFTYPE=CART")
        println(io, "  REFLAT=", _edi_sexagesimal(lat), "\n  REFLONG=", _edi_sexagesimal(lon), "\n  REFELEV=", @sprintf("%.2f", elev), "\n")
        a = isempty(rot) ? 0.0 : rot[1]
        @printf(io, ">HMEAS ID=1001.001 CHTYPE=HX X=0.0 Y=0.0 Z=0.0 AZM=%.2f\n", a)
        @printf(io, ">HMEAS ID=1002.001 CHTYPE=HY X=0.0 Y=0.0 Z=0.0 AZM=%.2f\n", a + 90)
        tipper && println(io, ">HMEAS ID=1003.001 CHTYPE=HZ X=0.0 Y=0.0 Z=0.0 AZM=0.00")
        @printf(io, ">EMEAS ID=1004.001 CHTYPE=EX X=-50.0 Y=0.0 Z=0.0 X2=50.0 Y2=0.0 AZM=%.2f\n", a)
        @printf(io, ">EMEAS ID=1005.001 CHTYPE=EY X=0.0 Y=-50.0 Z=0.0 X2=0.0 Y2=50.0 AZM=%.2f\n\n", a + 90)
        println(io, ">=MTSECT\n  SECTID=\"", d.site[is], "\"\n  NFREQ=", n)
        println(io, "  HX=1001.001\n  HY=1002.001", tipper ? "\n  HZ=1003.001" : "", "\n  EX=1004.001\n  EY=1005.001\n")
        _edi_write_block(io, ">FREQ ORDER=DEC", 1 ./ d.T[keep])
        _edi_write_block(io, ">ZROT", rot)
        for (ic, c) in enumerate(_EDI_ZCOMP)
            z = d.Z[keep, ic, is] ./ _EDI_FIELD
            e = abs.(d.Zerr[keep, ic, is]) ./ _EDI_FIELD
            _edi_write_block(io, ">$(c)R ROT=ZROT", val.(real.(z)))
            _edi_write_block(io, ">$(c)I ROT=ZROT", val.(imag.(z)))
            _edi_write_block(io, ">$(c).VAR ROT=ZROT", [isfinite(z[k]) && isfinite(e[k]) ? e[k]^2 : empty for k in 1:n])
        end
        if tipper
            _edi_write_block(io, ">TROT", trot)
            for (ic, c) in enumerate(_EDI_TCOMP)
                t, e = d.tip[keep, ic, is], abs.(d.tiperr[keep, ic, is])
                _edi_write_block(io, ">$(c)R.EXP ROT=TROT", val.(real.(t)))
                _edi_write_block(io, ">$(c)I.EXP ROT=TROT", val.(imag.(t)))
                _edi_write_block(io, ">$(c)VAR.EXP ROT=TROT", [isfinite(t[k]) && isfinite(e[k]) ? e[k]^2 : empty for k in 1:n])
            end
        end
        println(io, ">END")
    end
    path
end

"""
    ModEMToEDI(data_path, output_dir=""; sites=nothing) -> Vector{String}

Write one EDI per site of a ModEM data file (impedance and tipper sections,
[mV/km]/[nT], exp(+iωt), variances the squared errors, frequencies decreasing). Only
the periods a site holds are written, missing components as EMPTY; the header rotation
becomes ZROT/TROT. ModEM's 1e12 "no error" entries come out without a variance. The
default `output_dir` is `<data>_EDI` next to the data file; `sites` limits the output.
"""
function ModEMToEDI(data_path::AbstractString, output_dir::AbstractString = ""; sites = nothing)
    d = redirect_stdout(() -> load_data_modem(data_path; warn_rotation = false), devnull)
    for A in (d.Zerr, d.tiperr), i in eachindex(A)
        abs(A[i]) > 1e10 && (A[i] = complex(NaN, NaN))
    end
    dir = isempty(output_dir) ? splitext(abspath(data_path))[1] * "_EDI" : String(output_dir)
    mkpath(dir)
    which = sites === nothing ? (1:d.ns) : [something(findfirst(==(s), d.site), 0) for s in sites]
    any(==(0), which) && error("sites not in $data_path: $(join(sites[which .== 0], ", "))")
    paths = [_write_edi(joinpath(dir, d.site[is] * ".edi"), d, is, basename(data_path)) for is in which]
    println("ModEMToEDI: $(length(paths)) EDI files in $dir")
    paths
end

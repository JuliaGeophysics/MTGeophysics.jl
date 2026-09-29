# MT data dashboard and editor
# Author: @pankajkmishra
# One window to page through the sites of a ModEM data file (EDIs go through EDIToModEM first) and mask data: apparent
# resistivity above phase of Zxy (red) and Zyx (blue), with Zxx and Zyy (paler) joining them under "Full tensor", every
# phase as recorded (Zyx in the third quadrant) on a fixed -200..200°; the real (red) and imaginary (blue) tipper sit in a
# panel that collapses to the right, a small site map in one that collapses to the left, and a ModEM response can be
# drawn against the data. Masking is by drag only, and the full tensor switch only changes the view: a drag selects a
# band of periods across the full height of ρa and φ and masks all four impedances in it (across the tipper, Tzx and
# Tzy), or restores them with Shift held; a click on a data panel does nothing. The mask is saved as a text file that apply_data_mask applies to any copy of the survey
# (DataMask.jl), and the kept data are written with their own errors, under a new time-stamped name each time, as a
# ModEM file or as EDIs. ModEM's 1e12 "no error" is read as missing. The window is built on any Makie backend, so
# CairoMakie renders it headless; DataDashboard shows it with GLMakie

const _DB_ZCOLOURS = (RGBf(0.62, 0.84, 0.60), RGBf(0.84, 0.30, 0.10), RGBf(0.12, 0.38, 0.72), RGBf(0.78, 0.66, 0.88))
const _DB_ZNAMES = ("Zxx", "Zxy", "Zyx", "Zyy")
const _DB_TNAMES = ("Tzx", "Tzy")
const _DB_RE = (colour = RGBf(0.84, 0.30, 0.10), marker = :circle, linestyle = :solid)
const _DB_IM = (colour = RGBf(0.12, 0.38, 0.72), marker = :circle, linestyle = :dash)
const _DB_MASKED = (markersize = 9, color = :transparent, strokecolor = :grey55, strokewidth = 1.2)
const _DB_PHASE_RANGE = (-200.0, 200.0)
const _DB_CLICK_PX = 12
const _DB_DRAG_PX = 5

struct _Decades end
Makie.get_tickvalues(::_Decades, vmin, vmax) = collect(Float64, ceil(vmin):floor(vmax))
const _DB_DECADES = LogTicks(_Decades())

# whole decades as 10ⁿ when at least two fall in the range, else 1-2-5 steps in plain numbers
function _db_log_ticks(lo, hi)
    decades = ceil(Int, log10(lo)):floor(Int, log10(hi))
    length(decades) >= 2 && return (10.0 .^ decades, [rich("10", superscript(string(n))) for n in decades])
    v = [m * 10.0^n for n in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 5) if lo <= m * 10.0^n <= hi]
    (v, _mt_log_labels(v))
end

_db_rho(z, T) = abs2(z) / (2π / T * 4π * 1e-7)
_db_phase(z) = rad2deg(angle(z))
_db_present(st, ip, ic, is) = ic <= 4 ? isfinite(st.obs.Z[ip, ic, is]) : isfinite(st.obs.tip[ip, ic - 4, is])

function _db_break_wraps(x::AbstractVector, φ::AbstractVector)
    xs, ys = Float64[], Float64[]
    for k in eachindex(φ)
        k > firstindex(φ) && abs(φ[k] - φ[k-1]) > 180 && (push!(xs, NaN); push!(ys, NaN))
        push!(xs, x[k]); push!(ys, φ[k])
    end
    xs, ys
end

function _db_align_predicted(obs::Data, pred::Data)
    Z = fill(complex(NaN, NaN), obs.nf, 4, obs.ns)
    tip = fill(complex(NaN, NaN), obs.nf, 2, obs.ns)
    index = Dict(uppercase(strip(s)) => j for (j, s) in enumerate(pred.site))
    matched, turned = 0, 0
    for is in 1:obs.ns
        j = get(index, uppercase(strip(obs.site[is])), nothing)
        j === nothing && continue
        matched += 1
        for (k, t) in enumerate(pred.T)
            ip = _period_index(obs.T, t; rtol = 0.02)
            ip === nothing && continue
            z, tz = pred.Z[k, :, j], pred.tip[k, :, j]
            θ = obs.zrot[ip, is] - (isempty(pred.zrot) ? 0.0 : pred.zrot[k, j])
            if abs(θ) > 1e-3
                r = _rotate_impedance(permutedims(reshape(z, 2, 2)), zeros(2, 2), θ)
                z = r === nothing ? fill(complex(NaN, NaN), 4) : vec(permutedims(r.Z))
                r = _rotate_tipper(tz, zeros(2), θ)
                tz = r === nothing ? fill(complex(NaN, NaN), 2) : r.T
                turned += 1
            end
            Z[ip, :, is] .= z
            tip[ip, :, is] .= tz
        end
    end
    matched == 0 && @warn "no predicted site matches an observed site by name"
    println("  predicted: $matched of $(obs.ns) sites matched", turned > 0 ? ", $turned site-periods rotated to the observed frame" : "")
    (Z = Z, tip = tip)
end

function _db_site_rms(st, is)
    st.pred === nothing && return (rms = NaN, n = 0)
    num, n = 0.0, 0
    for ic in 1:6, ip in 1:st.obs.nf
        st.keep[ip, ic, is] || continue
        o, p, e = ic <= 4 ? (st.obs.Z[ip, ic, is], st.pred.Z[ip, ic, is], abs(st.obs.Zerr[ip, ic, is])) :
                            (st.obs.tip[ip, ic - 4, is], st.pred.tip[ip, ic - 4, is], abs(st.obs.tiperr[ip, ic - 4, is]))
        (isfinite(o) && isfinite(p) && isfinite(e) && e > 0) || continue
        num += abs2(real(p - o) / e) + abs2(imag(p - o) / e)
        n += 2
    end
    n > 0 ? (rms = sqrt(num / n), n = n) : (rms = NaN, n = 0)
end

function _db_total_rms(st)
    num, n = 0.0, 0
    for is in 1:st.obs.ns
        r = _db_site_rms(st, is)
        r.n > 0 && (num += r.rms^2 * r.n; n += r.n)
    end
    n > 0 ? sqrt(num / n) : NaN
end

# the kept data with their own errors, on the sites and periods that keep any
function _db_export_data(st)
    d = deepcopy(st.obs)
    for is in 1:d.ns, ip in 1:d.nf
        for ic in 1:4
            st.keep[ip, ic, is] || (d.Z[ip, ic, is] = d.Zerr[ip, ic, is] = complex(NaN, NaN))
        end
        for ic in 1:2
            st.keep[ip, ic + 4, is] || (d.tip[ip, ic, is] = d.tiperr[ip, ic, is] = complex(NaN, NaN))
        end
    end
    sites = [is for is in 1:d.ns if any(st.keep[:, :, is])]
    pers = [ip for ip in 1:d.nf if any(st.keep[ip, :, sites])]
    isempty(sites) && error("every data point is masked")
    d.site, d.loc, d.x, d.y, d.z = d.site[sites], d.loc[sites, :], d.x[sites], d.y[sites], d.z[sites]
    d.T, d.f = d.T[pers], d.f[pers]
    d.Z, d.Zerr, d.tip, d.tiperr = d.Z[pers, :, sites], d.Zerr[pers, :, sites], d.tip[pers, :, sites], d.tiperr[pers, :, sites]
    d.zrot, d.trot = d.zrot[pers, sites], d.trot[pers, sites]
    d.ns, d.nf = length(sites), length(pers)
    d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)
    d.name = ""
    length(unique(d.zrot)) > 1 && @warn "sites carry different ZROT; the ModEM header takes $(d.zrot[1])°"
    d
end

_db_stamp() = Dates.format(now(), "yyyymmdd_HHMMSS")

# the kept data as <stem>_<date_time>.dat in dir; a datum without an error gets ModEM's 1e12 back
function _db_write_modem(st, dir::AbstractString, stem::AbstractString)
    path = joinpath(dir, "$(stem)_$(_db_stamp()).dat")
    redirect_stdout(devnull) do
        write_data_modem(path, _db_export_data(st); sign = 1, units = "[mV/km]/[nT]", include_tipper = any(st.keep[:, 5:6, :]),
                         description = "DataDashboard export of $stem, recorded errors")
    end
    println("Wrote $path")
    path
end

# the kept data as one EDI per site in <stem>_EDI_<date_time>/ in dir
function _db_write_edi(st, dir::AbstractString, stem::AbstractString)
    d = _db_export_data(st)
    out = mkpath(joinpath(dir, "$(stem)_EDI_$(_db_stamp())"))
    foreach(is -> _write_edi(joinpath(out, d.site[is] * ".edi"), d, is, stem), 1:d.ns)
    println("Wrote $(d.ns) EDI files to $out")
    out
end

_db_period_axis(pos; kw...) = _mt_axis(pos; xscale = log10, xticks = _DB_DECADES, xlabel = "Period (s)", kw...)

function _db_limits!(ax, xl, yl = nothing)
    ax.limits[] = (xl, yl)
    reset_limits!(ax)
end

_db_note!(ax, s) = text!(ax, 0.02, 0.96; text = s, space = :relative, align = (:left, :top), fontsize = 12, color = :grey30)

#---------- pixels and hits ----------

_db_frac(v, lo, hi, lg) = lg ? (log10(v) - log10(lo)) / (log10(hi) - log10(lo)) : (v - lo) / (hi - lo)
_db_unfrac(f, lo, hi, lg) = lg ? 10^(log10(lo) + f * (log10(hi) - log10(lo))) : lo + f * (hi - lo)

# window pixels of a data point, and back; finallimits are in data units on log axes too
function _db_px(ax, x, y)
    vp, lim = ax.scene.viewport[], ax.finallimits[]
    lo, hi = lim.origin, lim.origin .+ lim.widths
    Point2f(vp.origin[1] + vp.widths[1] * _db_frac(x, lo[1], hi[1], ax.xscale[] === log10),
            vp.origin[2] + vp.widths[2] * _db_frac(y, lo[2], hi[2], ax.yscale[] === log10))
end

function _db_data(ax, p)
    vp, lim = ax.scene.viewport[], ax.finallimits[]
    lo, hi = lim.origin, lim.origin .+ lim.widths
    f = (p .- vp.origin) ./ vp.widths
    Point2f(_db_unfrac(f[1], lo[1], hi[1], ax.xscale[] === log10), _db_unfrac(f[2], lo[2], hi[2], ax.yscale[] === log10))
end

# every (component, period) drawn on the given axes with a period between lo and hi
function _db_in_band(hits, axes, lo, hi)
    lo, hi = minmax(lo, hi)
    out = Tuple{Int, Int}[]
    for h in hits
        any(ax -> ax === h.ax, axes) || continue
        for k in eachindex(h.ip)
            lo <= h.x[k] <= hi && push!(out, (h.ic, h.ip[k]))
        end
    end
    unique(out)
end

# keep or mask (component, period) pairs at one site; comps, when given, replaces the picked component by these at
# each picked period
function _db_set!(st, is, picks, value::Bool, comps = nothing)
    n = 0
    for (ic, ip) in picks, c in something(comps, (ic,))
        _db_present(st, ip, c, is) || continue
        n += st.keep[ip, c, is] != value
        st.keep[ip, c, is] = value
    end
    n
end

#---------- panels ----------

# the points an axis range covers: the kept ones, so masking an outlier zooms in; all of them when none is kept
_db_scaled(valid, keep) = (k = filter(i -> keep[i], valid); isempty(k) ? valid : k)

# one observed series: error bars and filled markers where kept, hollow markers where masked
function _db_series!(ax, T, y, lo, hi, keep, colour, marker)
    valid = findall(isfinite, y)
    kept, dropped = filter(i -> keep[i], valid), filter(i -> !keep[i], valid)
    e = filter(i -> isfinite(lo[i]) && isfinite(hi[i]), kept)
    isempty(e) || errorbars!(ax, T[e], y[e], lo[e], hi[e]; color = :black, linewidth = 1, whiskerwidth = 6)
    scatter!(ax, T[kept], y[kept]; _MT_MARKER..., color = colour, marker = marker)
    scatter!(ax, T[dropped], y[dropped]; _DB_MASKED..., marker = marker)
    valid
end

function _db_impedance!(axρ, axφ, st, is, comps, hits)
    empty!(axρ); empty!(axφ)
    T = st.obs.T
    ρseen, φseen = Float64[], Float64[]
    for ic in comps
        colour = _DB_ZCOLOURS[ic]
        if st.pred !== nothing
            zp = st.pred.Z[:, ic, is]
            ok = findall(isfinite, zp)
            if !isempty(ok)
                ρp = _db_rho.(zp[ok], T[ok])
                lines!(axρ, T[ok], ρp; color = colour, linewidth = _MT_LINEWIDTH)
                lines!(axφ, _db_break_wraps(T[ok], _db_phase.(zp[ok]))...; color = colour, linewidth = _MT_LINEWIDTH)
                append!(ρseen, ρp)
            end
        end
        z = st.obs.Z[:, ic, is]
        ρ = [isfinite(v) && abs(v) > 0 ? _db_rho(v, t) : NaN for (v, t) in zip(z, T)]
        φ = [isfinite(r) ? _db_phase(v) : NaN for (v, r) in zip(z, ρ)]
        rel = abs.(st.obs.Zerr[:, ic, is]) ./ abs.(z)
        δρ, δφ = 2 .* ρ .* rel, min.(rad2deg.(rel), 90)
        keep = st.keep[:, ic, is]
        valid = _db_series!(axρ, T, ρ, min.(δρ, 0.95 .* ρ), δρ, keep, colour, :circle)
        _db_series!(axφ, T, φ, δφ, δφ, keep, colour, :circle)
        push!(hits, (ax = axρ, ic = ic, ip = valid, x = T[valid], y = ρ[valid]))
        push!(hits, (ax = axφ, ic = ic, ip = valid, x = T[valid], y = φ[valid]))
        shown = _db_scaled(valid, keep)
        append!(ρseen, ρ[shown]); append!(φseen, φ[shown])
    end
    isempty(φseen) && _db_note!(axρ, "no $(join(_DB_ZNAMES[collect(comps)], ", ")) at this site")
    ρl = isempty(ρseen) ? (1.0, 1e4) : (10^(log10(minimum(ρseen)) - 0.5), 10^(log10(maximum(ρseen)) + 0.5))
    axρ.yticks = _db_log_ticks(ρl...)
    _db_limits!(axρ, st.Trange, ρl)
    _db_limits!(axφ, st.Trange, _DB_PHASE_RANGE)
end

function _db_tipper!(ax, st, is, j, hits)
    empty!(ax)
    T = st.obs.T
    t = st.obs.tip[:, j, is]
    if !any(isfinite, t)
        _db_note!(ax, "no tipper at this site")
        return _db_limits!(ax, st.Trange, (-0.5, 0.5))
    end
    hlines!(ax, [0.0]; color = :grey70, linewidth = 1)
    e = abs.(st.obs.tiperr[:, j, is])
    keep = st.keep[:, j + 4, is]
    shown = _db_scaled(findall(isfinite, t), keep)
    top = maximum(abs, vcat(real.(t[shown]), imag.(t[shown])); init = 0.0)
    for (part, sty) in ((real, _DB_RE), (imag, _DB_IM))
        if st.pred !== nothing
            tp = st.pred.tip[:, j, is]
            ok = findall(isfinite, tp)
            isempty(ok) || lines!(ax, T[ok], part.(tp[ok]); color = sty.colour, linewidth = _MT_LINEWIDTH, linestyle = sty.linestyle)
        end
        y = [isfinite(v) ? part(v) : NaN for v in t]
        valid = _db_series!(ax, T, y, e, e, keep, sty.colour, sty.marker)
        push!(hits, (ax = ax, ic = j + 4, ip = valid, x = T[valid], y = y[valid]))
    end
    _db_limits!(ax, st.Trange, 1.1 .* (-max(top, 0.2), max(top, 0.2)))
end

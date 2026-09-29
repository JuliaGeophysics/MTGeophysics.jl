# MT data dashboard
# Author: @pankajkmishra
# One window to page through the sites of a ModEM data file (EDIs go through EDIToModEM first): apparent resistivity and
# phase of all four components, the yx phase turned by 180° into the first quadrant, observed against a ModEM response,
# and four panels to choose from (tipper, phase tensor ellipses and induction arrows along period, PT skew, ellipticity
# and Swift/Bahr skew, strike, normalised residuals, relative errors, Niblett-Bostick). The Swift strike takes the
# branch that maximises |Zxy|² + |Zyx|²; log axes are labelled at whole decades. A site map shows the per-site RMS and
# the phase tensors at one period. Points can be masked by clicking and the kept data written as a ModEM file with
# error floors; ModEM's 1e12 "no error" is read as missing. The window is built on any Makie backend, so CairoMakie
# renders it headless; DataDashboard shows it with GLMakie
# References: Swift (1967); Bahr (1988); Caldwell et al. (2004); Niblett & Sayn-Wittgenstein (1960), Bostick (1977)

const _DB_ZCOLOURS = (RGBf(0.20, 0.60, 0.30), RGBf(0.12, 0.38, 0.72), RGBf(0.84, 0.30, 0.10), RGBf(0.55, 0.30, 0.65))
const _DB_ZNAMES = ("XX", "XY", "YX", "YY")
const _DB_TCOLOURS = (RGBf(0.12, 0.38, 0.72), RGBf(0.84, 0.30, 0.10))
const _DB_MASKED = (marker = :circle, markersize = 9, color = :transparent, strokecolor = :grey55, strokewidth = 1.2)
const _DB_PANELS = [
    (:tipper, "Tipper"),
    (:ptstrip, "Phase tensors & induction arrows"),
    (:beta, "PT skew β"),
    (:skew, "Ellipticity, Swift & Bahr skew"),
    (:strike, "Strike"),
    (:resid, "Normalised residuals"),
    (:relerr, "Relative errors"),
    (:nb, "Niblett–Bostick"),
]
const _DB_PT_COLOURMAP = Reverse(:Spectral)

struct _Decades end
Makie.get_tickvalues(::_Decades, vmin, vmax) = collect(Float64, ceil(vmin):floor(vmax))
const _DB_DECADES = LogTicks(_Decades())

_db_wrap180(a) = mod(a + 180, 360) - 180
_db_rho(z, T) = abs2(z) / (2π / T * 4π * 1e-7)
_db_phase(z, ic) = ic == 3 ? _db_wrap180(rad2deg(angle(z)) + 180) : rad2deg(angle(z))

function _db_break_wraps(x::AbstractVector, φ::AbstractVector)
    xs, ys = Float64[], Float64[]
    for k in eachindex(φ)
        k > firstindex(φ) && abs(φ[k] - φ[k-1]) > 180 && (push!(xs, NaN); push!(ys, NaN))
        push!(xs, x[k]); push!(ys, φ[k])
    end
    xs, ys
end

function _db_swift_bahr(zxx, zxy, zyx, zyy)
    all(isfinite, (zxx, zxy, zyx, zyy)) || return (κ = NaN, η = NaN, strike = NaN)
    S1, S2, D1, D2 = zxx + zyy, zxy + zyx, zxx - zyy, zxy - zyx
    comm(a, b) = real(a) * imag(b) - real(b) * imag(a)
    (κ = abs(S1) / abs(D2), η = sqrt(abs(comm(D1, S2) - comm(S1, D2))) / abs(D2),
     strike = atand(-2 * real(D1 * conj(S2)), abs2(S2) - abs2(D1)) / 4)
end

function _db_bostick(z::AbstractVector, T::AbstractVector)
    ρa = _db_rho.(z, T)
    φ = mod.(angle.(z), π / 2)
    depth = sqrt.(ρa .* T ./ (2π * 4π * 1e-7))
    ρb = ρa .* (π ./ (2 .* φ) .- 1)
    ok = isfinite.(ρb) .& (ρb .> 0) .& isfinite.(depth)
    depth[ok], ρb[ok]
end

function _db_zdet(Z::AbstractMatrix)
    map(axes(Z, 1)) do i
        zxx, zxy, zyx, zyy = Z[i, :]
        isfinite(zxx) && isfinite(zyy) || (zxx = zyy = 0.0im)
        z = sqrt(zxx * zyy - zxy * zyx)
        real(z) < 0 ? -z : z
    end
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

function _db_residuals(st, is)
    st.pred === nothing && return nothing
    re, im = fill(NaN, st.obs.nf, 6), fill(NaN, st.obs.nf, 6)
    for ic in 1:6, ip in 1:st.obs.nf
        st.keep[ip, ic, is] || continue
        o, p, e = ic <= 4 ? (st.obs.Z[ip, ic, is], st.pred.Z[ip, ic, is], st.ez[ip, ic, is]) :
                            (st.obs.tip[ip, ic - 4, is], st.pred.tip[ip, ic - 4, is], st.et[ip, ic - 4, is])
        (isfinite(p) && isfinite(e) && e > 0) || continue
        re[ip, ic], im[ip, ic] = real(p - o) / e, imag(p - o) / e
    end
    (re = re, im = im)
end

function _db_site_rms(st, is)
    r = _db_residuals(st, is)
    r === nothing && return (rms = NaN, n = 0)
    v = filter(isfinite, vcat(vec(r.re), vec(r.im)))
    isempty(v) ? (rms = NaN, n = 0) : (rms = sqrt(sum(abs2, v) / length(v)), n = length(v))
end

function _db_total_rms(st)
    num, n = 0.0, 0
    for is in 1:st.obs.ns
        r = _db_site_rms(st, is)
        r.n > 0 && (num += r.rms^2 * r.n; n += r.n)
    end
    n > 0 ? sqrt(num / n) : NaN
end

function _db_pt(st, is; predicted::Bool = false)
    Z = predicted ? st.pred.Z : st.obs.Z
    map(1:st.obs.nf) do ip
        all(st.keep[ip, 1:4, is]) || return nothing
        phase_tensor(Z[ip, 1, is], Z[ip, 2, is], Z[ip, 3, is], Z[ip, 4, is])
    end
end

function _db_export_data(st)
    d = deepcopy(st.obs)
    for is in 1:d.ns, ip in 1:d.nf
        for ic in 1:4
            st.keep[ip, ic, is] || (d.Z[ip, ic, is] = complex(NaN, NaN))
            d.Zerr[ip, ic, is] = complex(st.ez[ip, ic, is], 0.0)
        end
        for ic in 1:2
            st.keep[ip, ic + 4, is] || (d.tip[ip, ic, is] = complex(NaN, NaN))
            d.tiperr[ip, ic, is] = complex(st.et[ip, ic, is], 0.0)
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

_db_period_axis(pos; kw...) = _mt_axis(pos; xscale = log10, xticks = _DB_DECADES, xlabel = "Period (s)", kw...)

function _db_decades(v, default)
    v = filter(x -> isfinite(x) && x > 0, v)
    isempty(v) && return default
    lo = floor(log10(minimum(v)))
    (10.0^lo, 10.0^max(ceil(log10(maximum(v))), lo + 1))
end

function _db_limits!(ax, xl, yl = nothing)
    ax.limits[] = (xl, yl)
    reset_limits!(ax)
end

_db_note!(ax, s) = text!(ax, 0.02, 0.96; text = s, space = :relative, align = (:left, :top), fontsize = 12, color = :grey30)

function _db_rhophi!(axρ, axφ, st, is, diagonals::Bool, pickmap)
    empty!(axρ); empty!(axφ)
    T = st.obs.T
    φlo, φhi, ρseen = 0.0, 90.0, Float64[]
    for ic in (diagonals ? (1, 2, 3, 4) : (2, 3))
        colour = _DB_ZCOLOURS[ic]
        if st.pred !== nothing
            zp = st.pred.Z[:, ic, is]
            ok = findall(isfinite.(zp))
            if !isempty(ok)
                lines!(axρ, T[ok], _db_rho.(zp[ok], T[ok]); color = colour, linewidth = _MT_LINEWIDTH)
                append!(ρseen, _db_rho.(zp[ok], T[ok]))
                lines!(axφ, _db_break_wraps(T[ok], _db_phase.(zp[ok], ic))...; color = colour, linewidth = _MT_LINEWIDTH)
            end
        end
        z = st.obs.Z[:, ic, is]
        ρ, φ = st.obs.ρ[:, ic, is], _db_phase.(z, ic)
        valid = isfinite.(z) .& (ρ .> 0)
        kept = findall(valid .& st.keep[:, ic, is])
        dropped = findall(valid .& .!st.keep[:, ic, is])
        rel = abs.(st.obs.Zerr[kept, ic, is]) ./ abs.(z[kept])
        e = findall(isfinite, rel)
        if !isempty(e)
            δρ = 2 .* ρ[kept[e]] .* rel[e]
            errorbars!(axρ, T[kept[e]], ρ[kept[e]], min.(δρ, 0.95 .* ρ[kept[e]]), δρ; color = :black, linewidth = 1, whiskerwidth = 6)
            append!(ρseen, ρ[kept[e]] .+ δρ, ρ[kept[e]] .- min.(δρ, 0.95 .* ρ[kept[e]]))
            errorbars!(axφ, T[kept[e]], φ[kept[e]], min.(rad2deg.(rel[e]), 90); color = :black, linewidth = 1, whiskerwidth = 6)
        end
        for (ax, y) in ((axρ, ρ), (axφ, φ))
            p = scatter!(ax, T[kept], y[kept]; color = colour, _MT_MARKER...)
            pickmap[objectid(p)] = (ic, kept)
            q = scatter!(ax, T[dropped], y[dropped]; _DB_MASKED...)
            pickmap[objectid(q)] = (ic, dropped)
        end
        append!(ρseen, ρ[kept], ρ[dropped])
        fin = filter(isfinite, φ[valid])
        isempty(fin) || (φlo = min(φlo, minimum(fin)); φhi = max(φhi, maximum(fin)))
    end
    _db_limits!(axρ, st.Trange, _db_decades(ρseen, (1.0, 1e4)))
    _db_limits!(axφ, st.Trange, (max(φlo, -180) - 5, min(φhi, 180) + 5))
end

function _db_tipper!(ax, st, is, pickmap)
    T = st.obs.T
    any(isfinite, st.obs.tip[:, :, is]) || return _db_note!(ax, "no tipper at this site")
    hlines!(ax, [0.0]; color = :grey70, linewidth = 1)
    for ic in 1:2, (part, marker, style) in ((real, :circle, :solid), (imag, :utriangle, :dash))
        colour = _DB_TCOLOURS[ic]
        if st.pred !== nothing
            tp = st.pred.tip[:, ic, is]
            ok = findall(isfinite.(tp))
            isempty(ok) || lines!(ax, T[ok], part.(tp[ok]); color = colour, linewidth = _MT_LINEWIDTH, linestyle = style)
        end
        t = st.obs.tip[:, ic, is]
        valid = isfinite.(t)
        kept = findall(valid .& st.keep[:, ic + 4, is])
        dropped = findall(valid .& .!st.keep[:, ic + 4, is])
        e = abs.(st.obs.tiperr[kept, ic, is])
        has = findall(isfinite, e)
        isempty(has) || errorbars!(ax, T[kept[has]], part.(t[kept[has]]), e[has]; color = :black, linewidth = 1, whiskerwidth = 6)
        p = scatter!(ax, T[kept], part.(t[kept]); color = colour, _MT_MARKER..., marker = marker)
        pickmap[objectid(p)] = (ic + 4, kept)
        q = scatter!(ax, T[dropped], part.(t[dropped]); _DB_MASKED..., marker = marker)
        pickmap[objectid(q)] = (ic + 4, dropped)
    end
    ax.ylabel = "Tipper"
    axislegend(ax, [MarkerElement(; color = _DB_TCOLOURS[1], _MT_MARKER...), MarkerElement(; color = _DB_TCOLOURS[2], _MT_MARKER...),
                    MarkerElement(; color = :white, _MT_MARKER..., marker = :utriangle)],
               ["Tzx", "Tzy", "imaginary"]; position = :lt, framevisible = false, labelsize = 11, patchsize = (16, 12))
    _db_limits!(ax, st.Trange)
end

function _db_ptstrip!(ax, st, is, colour_slot)
    T = st.obs.T
    have = findall(ip -> any(isfinite, st.obs.Z[ip, :, is]) || any(isfinite, st.obs.tip[ip, :, is]), 1:st.obs.nf)
    isempty(have) && return _db_note!(ax, "no data at this site")
    lt = log10.(T[have])
    shown = Set(have[1:cld(length(have), 30):end])
    Δ = length(shown) > 1 ? (lt[end] - lt[1]) / (length(shown) - 1) : 1.0
    xs(ip) = log10(T[ip]) / Δ
    rows = st.pred === nothing ? [(:pt, false, 1.1, "PT observed"), (:iv, false, 0.0, "Arrows observed")] :
           [(:pt, true, 3.3, "PT predicted"), (:pt, false, 2.2, "PT observed"), (:iv, false, 1.1, "Arrows observed"), (:iv, true, 0.0, "Arrows predicted")]
    polys, cols = Vector{Point2f}[], Float32[]
    lines_re, lines_im, heads_re, heads_im = Point2f[], Point2f[], Vector{Point2f}[], Vector{Point2f}[]
    for (kind, predicted, y, _) in rows
        if kind == :pt
            for (ip, pt) in enumerate(_db_pt(st, is; predicted))
                (pt === nothing || !(ip in shown)) && continue
                a, b = _ellipse_semiaxes(pt, 0.9, 1.0)
                a > 0 || continue
                ex, ey = _ellipse_ring(xs(ip), y, a, b, pt.azimuth + st.obs.zrot[ip, is])
                push!(polys, Point2f.(ex, ey)); push!(cols, Float32(atand(pt.phimin)))
            end
        else
            tip = predicted ? st.pred.tip : st.obs.tip
            for ip in 1:st.obs.nf
                (ip in shown && all(st.keep[ip, 5:6, is])) || continue
                iv = induction_vector(tip[ip, 1, is], tip[ip, 2, is]; convention = st.convention)
                iv === nothing && continue
                for (vec, store, heads) in ((iv.re, lines_re, heads_re), (iv.im, lines_im, heads_im))
                    (sx, sy), (hx, hy) = _arrow_parts(xs(ip), y, 1.6 * vec[1], 1.6 * vec[2])
                    isempty(sx) && continue
                    append!(store, Point2f.(sx, sy)); push!(store, Point2f(NaN, NaN))
                    push!(heads, Point2f.(hx, hy))
                end
            end
        end
    end
    isempty(polys) && isempty(lines_re) && _db_note!(ax, "needs the full impedance tensor or both tipper components")
    isempty(polys) || poly!(ax, polys; color = cols, colormap = _DB_PT_COLOURMAP, colorrange = (0, 90), strokecolor = :black, strokewidth = 0.8)
    lines!(ax, lines_im; color = :grey55, linewidth = 1.4); poly!(ax, heads_im; color = :grey55)
    lines!(ax, lines_re; color = :black, linewidth = 1.4); poly!(ax, heads_re; color = :black)
    decades = floor(Int, minimum(lt)):ceil(Int, maximum(lt))
    ax.xticks = (collect(decades) ./ Δ, [rich("10", superscript(string(n))) for n in decades])
    ax.yticks = ([r[3] for r in rows], [r[4] for r in rows])
    ax.aspect = DataAspect()
    ax.xlabel = "Period (s)"
    _db_limits!(ax, (log10(st.Trange[1]) / Δ, log10(st.Trange[2]) / Δ), (-1.2, rows[1][3] + 0.8))
    colour_slot(Colorbar, (; colormap = _DB_PT_COLOURMAP, limits = (0, 90), label = "Φmin (°)", width = 10))
end

function _db_beta!(ax, st, is)
    hspan!(ax, -3, 3; color = (:grey, 0.15))
    hlines!(ax, [0.0]; color = :grey70, linewidth = 1)
    T = st.obs.T
    if st.pred !== nothing
        b = [pt === nothing ? NaN : pt.beta for pt in _db_pt(st, is; predicted = true)]
        ok = findall(isfinite, b)
        lines!(ax, T[ok], b[ok]; color = :black, linewidth = _MT_LINEWIDTH)
    end
    b = [pt === nothing ? NaN : pt.beta for pt in _db_pt(st, is)]
    ok = findall(isfinite, b)
    scatter!(ax, T[ok], b[ok]; color = _DB_ZCOLOURS[4], _MT_MARKER...)
    isempty(ok) && _db_note!(ax, "needs the full impedance tensor")
    ax.ylabel = "β (°)"
    _db_limits!(ax, st.Trange)
end

function _db_skews!(ax, st, is)
    T = st.obs.T
    hlines!(ax, [0.3]; color = :grey60, linestyle = :dash, linewidth = 1)
    series(Z, pts) = (ellipticity = [pt === nothing ? NaN : pt.ellipticity for pt in pts],
                      κ = [all(st.keep[ip, 1:4, is]) ? _db_swift_bahr(Z[ip, :, is]...).κ : NaN for ip in 1:st.obs.nf],
                      η = [all(st.keep[ip, 1:4, is]) ? _db_swift_bahr(Z[ip, :, is]...).η : NaN for ip in 1:st.obs.nf])
    names = (ellipticity = "PT ellipticity", κ = "Swift skew κ", η = "Bahr skew η")
    colours = (ellipticity = _DB_ZCOLOURS[4], κ = _DB_ZCOLOURS[1], η = RGBf(0.85, 0.65, 0.10))
    p = st.pred === nothing ? nothing : series(st.pred.Z, _db_pt(st, is; predicted = true))
    o = series(st.obs.Z, _db_pt(st, is))
    for k in keys(names)
        if p !== nothing
            ok = findall(isfinite, p[k])
            lines!(ax, T[ok], p[k][ok]; color = colours[k], linewidth = _MT_LINEWIDTH)
        end
        ok = findall(isfinite, o[k])
        scatter!(ax, T[ok], o[k][ok]; color = colours[k], _MT_MARKER..., label = names[k])
    end
    any(isfinite, o.κ) ? axislegend(ax; position = :lt, framevisible = false, labelsize = 11, patchsize = (16, 12)) :
                         _db_note!(ax, "needs the full impedance tensor")
    ax.ylabel = "Skew, ellipticity"
    _db_limits!(ax, st.Trange, (0, max(1.0, maximum(filter(isfinite, vcat(o.κ, o.η, [0.0]))) * 1.05)))
end

function _db_strike!(ax, st, is)
    T = st.obs.T
    fold(a) = mod(a, 90)
    for (predicted, Z) in (st.pred === nothing ? ((false, st.obs.Z),) : ((true, st.pred.Z), (false, st.obs.Z)))
        pts = _db_pt(st, is; predicted)
        αpt = [pt === nothing ? NaN : fold(pt.azimuth + st.obs.zrot[ip, is]) for (ip, pt) in enumerate(pts)]
        αsw = [all(st.keep[ip, 1:4, is]) ? fold(_db_swift_bahr(Z[ip, :, is]...).strike + st.obs.zrot[ip, is]) : NaN for ip in 1:st.obs.nf]
        for (α, colour, label) in ((αpt, _DB_ZCOLOURS[4], "PT azimuth"), (αsw, _DB_ZCOLOURS[1], "Swift"))
            ok = findall(isfinite, α)
            if predicted
                scatter!(ax, T[ok], α[ok]; color = colour, marker = :xcross, markersize = 9)
            else
                scatter!(ax, T[ok], α[ok]; color = colour, _MT_MARKER..., label = label)
            end
        end
    end
    any(p -> p isa Scatter && !isempty(p[1][]), ax.scene.plots) ? axislegend(ax; position = :lt, framevisible = false, labelsize = 11, patchsize = (16, 12)) :
        _db_note!(ax, "needs the full impedance tensor")
    ax.ylabel = "Strike (°, modulo 90)"
    ax.yticks = 0:15:90
    _db_limits!(ax, st.Trange, (0, 90))
end

function _db_resid!(ax, st, is)
    r = _db_residuals(st, is)
    r === nothing && return _db_note!(ax, "no predicted data; pass predicted = \"<ModEM response>.dat\"")
    T = st.obs.T
    hspan!(ax, -1, 1; color = (:grey, 0.15))
    hlines!(ax, [0.0]; color = :grey70, linewidth = 1)
    colours = (_DB_ZCOLOURS..., RGBf(0.35, 0.35, 0.35), RGBf(0.65, 0.65, 0.65))
    labels = ("Zxx", "Zxy", "Zyx", "Zyy", "Tzx", "Tzy")
    for ic in 1:6
        ok = findall(isfinite, r.re[:, ic])
        isempty(ok) && continue
        scatter!(ax, T[ok], r.re[ok, ic]; color = colours[ic], _MT_MARKER..., label = labels[ic])
        scatter!(ax, T[ok], r.im[ok, ic]; color = colours[ic], _MT_MARKER..., marker = :utriangle)
    end
    s = _db_site_rms(st, is)
    text!(ax, 0.98, 0.96; text = @sprintf("RMS %.2f (N = %d)", s.rms, s.n), space = :relative, align = (:right, :top), fontsize = 12)
    s.n > 0 && axislegend(ax; position = :lb, framevisible = false, labelsize = 10, patchsize = (14, 10), nbanks = 3)
    ax.ylabel = "(pred − obs) / error"
    _db_limits!(ax, st.Trange)
end

function _db_relerr!(ax, st, is)
    T = st.obs.T
    st.z_floor > 0 && hlines!(ax, [100 * st.z_floor]; color = :grey60, linestyle = :dash, linewidth = 1)
    any_data = false
    for ic in 1:4
        z = st.obs.Z[:, ic, is]
        e = 100 .* abs.(st.obs.Zerr[:, ic, is]) ./ abs.(z)
        ok = findall(i -> isfinite(e[i]) && e[i] > 0 && st.keep[i, ic, is], eachindex(e))
        isempty(ok) && continue
        any_data = true
        scatter!(ax, T[ok], e[ok]; color = _DB_ZCOLOURS[ic], _MT_MARKER..., label = "Z" * lowercase(_DB_ZNAMES[ic]))
    end
    for ic in 1:2
        e = 100 .* abs.(st.obs.tiperr[:, ic, is])
        ok = findall(i -> isfinite(e[i]) && e[i] > 0 && st.keep[i, ic + 4, is], eachindex(e))
        isempty(ok) && continue
        any_data = true
        scatter!(ax, T[ok], e[ok]; color = _DB_TCOLOURS[ic], _MT_MARKER..., marker = :diamond, label = ic == 1 ? "100·δTzx" : "100·δTzy")
    end
    ax.yscale = log10
    ax.yticks = _DB_DECADES
    ax.ylabel = "Relative error (%)"
    any_data ? axislegend(ax; position = :lt, framevisible = false, labelsize = 10, patchsize = (14, 10), nbanks = 2) :
               _db_note!(ax, "no error estimates")
    _db_limits!(ax, st.Trange, any_data ? nothing : (0.1, 100.0))
end

function _db_nb!(ax, st, is)
    T = st.obs.T
    series = ((2, "xy", _DB_ZCOLOURS[2]), (3, "yx", _DB_ZCOLOURS[3]), (0, "det", :black))
    pick(Z, ic) = ic == 0 ? _db_zdet(Z[:, :, is]) : Z[:, ic, is]
    anything, ρs, hs = false, Float64[], Float64[]
    for (ic, label, colour) in series
        if st.pred !== nothing
            h, ρ = _db_bostick(pick(st.pred.Z, ic), T)
            isempty(h) || lines!(ax, ρ, h ./ 1000; color = colour, linewidth = _MT_LINEWIDTH)
            append!(ρs, ρ); append!(hs, h ./ 1000)
        end
        z = pick(st.obs.Z, ic)
        kept = ic == 0 ? [all(st.keep[ip, 2:3, is]) for ip in 1:st.obs.nf] : st.keep[:, ic, is]
        h, ρ = _db_bostick(ifelse.(kept, z, complex(NaN, NaN)), T)
        isempty(h) && continue
        anything = true
        append!(ρs, ρ); append!(hs, h ./ 1000)
        scatter!(ax, ρ, h ./ 1000; color = colour, _MT_MARKER..., label = label)
    end
    ax.xscale = log10; ax.yscale = log10
    ax.xticks = _DB_DECADES; ax.yticks = _DB_DECADES
    ax.yreversed = true
    ax.xlabel = "Bostick resistivity (Ω·m)"; ax.ylabel = "Depth (km)"
    anything ? axislegend(ax; position = :rt, framevisible = false, labelsize = 11, patchsize = (16, 12)) :
               _db_note!(ax, "no impedance at this site")
    _db_limits!(ax, _db_decades(ρs, (1.0, 1e4)), _db_decades(hs, (0.1, 100.0)))
end

function _dashboard_window(st, stem::AbstractString; panels, export_path::AbstractString, figsize)
    d = st.obs
    fig = Figure(size = figsize)

    header = fig[1, 1:3] = GridLayout()
    b_first = Button(header[1, 1], label = "|<")
    b_prev = Button(header[1, 2], label = "< Prev")
    b_next = Button(header[1, 3], label = "Next >")
    b_last = Button(header[1, 4], label = ">|")
    menu_sites = Menu(header[1, 5], options = [(s, i) for (i, s) in enumerate(d.site)], default = d.site[1], width = 170)
    info = Label(header[1, 6], "", fontsize = 13, halign = :left, tellwidth = false)
    diag_toggle = Toggle(header[1, 7], active = any(isfinite, d.Z[:, [1, 4], :]))
    Label(header[1, 8], "Diagonals", fontsize = 12)
    legend_items = Any[MarkerElement(; color = _DB_ZCOLOURS[ic], _MT_MARKER...) for ic in 1:4]
    legend_labels = ["Z" * lowercase(n) for n in _DB_ZNAMES]
    if st.pred !== nothing
        push!(legend_items, LineElement(color = :black, linewidth = _MT_LINEWIDTH)); push!(legend_labels, "predicted")
    end
    push!(legend_items, MarkerElement(; _DB_MASKED...)); push!(legend_labels, "masked")
    Legend(header[1, 9], legend_items, legend_labels; orientation = :horizontal, framevisible = false, labelsize = 12, patchsize = (16, 12))

    left = fig[2:4, 1] = GridLayout()
    fr = _ptiv_frame(d, "EPSG:4326"; pt_scale = 0.8, iv_scale = 0.0)
    x0, x1, y0, y1 = fr.limits
    width = (x1 - x0) * fr.lon_scale
    if y1 - y0 < 0.5 * width
        c = (y0 + y1) / 2
        y0, y1 = c - width / 4, c + width / 4
    end
    map_ax = _mt_axis(left[1, 1:2]; xlabel = "Longitude (°)", ylabel = "Latitude (°)", aspect = AxisAspect(width / (y1 - y0)),
                      limits = (x0, x1, y0, y1), xtickformat = _plain_tickformat, ytickformat = _plain_tickformat)
    map_period = Observable(cld(d.nf, 2))
    map_polys, map_cols = Observable(Vector{Point2f}[]), Observable(Float32[])
    map_pt = Observable(true)
    poly!(map_ax, map_polys; color = map_cols, colormap = _DB_PT_COLOURMAP, colorrange = (0, 90), strokecolor = :black, strokewidth = 0.6, visible = map_pt)
    site_colour = Observable(fill(NaN, d.ns))
    rms_range = Observable((0.0, 3.0))
    site_plot = scatter!(map_ax, fr.site_x, fr.site_y; color = site_colour, colormap = :viridis, colorrange = rms_range,
                         nan_color = :grey60, strokecolor = :black, strokewidth = 0.8, markersize = 8)
    current = Observable(Point2f(fr.site_x[1], fr.site_y[1]))
    scatter!(map_ax, current; color = :transparent, strokecolor = :red, strokewidth = 2.5, markersize = 22)
    st.pred === nothing || Colorbar(left[1, 3], colormap = :viridis, limits = rms_range, label = "Site RMS", width = 10)

    map_controls = left[2, 1:3] = GridLayout()
    pt_toggle = Toggle(map_controls[1, 1], active = true)
    Label(map_controls[1, 2], "PT map", fontsize = 12)
    period_slider = Slider(map_controls[1, 3], range = 1:d.nf, startvalue = map_period[])
    period_label = Label(map_controls[1, 4], "", fontsize = 12, width = 95)

    edit_controls = left[3, 1:3] = GridLayout()
    edit_toggle = Toggle(edit_controls[1, 1], active = false)
    Label(edit_controls[1, 2], "Edit mask", fontsize = 12)
    b_site = Button(edit_controls[1, 3], label = "Drop site")
    b_reset = Button(edit_controls[1, 4], label = "Reset site")
    b_export = Button(edit_controls[2, 1:2], label = "Export ModEM")
    b_png = Button(edit_controls[2, 3:4], label = "Save PNG")
    status = Label(left[4, 1:3], "", fontsize = 11, halign = :left, tellwidth = false, word_wrap = true)
    colsize!(fig.layout, 1, Fixed(430))

    axρ = _db_period_axis(fig[2, 2]; yscale = log10, yticks = _DB_DECADES, ylabel = "Apparent resistivity (Ω·m)")
    axφ = _db_period_axis(fig[2, 3]; ylabel = "Phase (°)")
    linkxaxes!(axρ, axφ)
    slots = map(enumerate(((3, 2), (3, 3), (4, 2), (4, 3)))) do (k, (r, c))
        layout = fig[r, c] = GridLayout()
        menu = Menu(layout[1, 1], options = [(l, p) for (p, l) in _DB_PANELS],
                    default = Dict(_DB_PANELS)[panels[k]], width = 260, halign = :left, tellwidth = false)
        (layout = layout, menu = menu, kind = Ref(panels[k]), blocks = Any[])
    end
    rowsize!(fig.layout, 2, Relative(0.34))

    cur = Ref(1)
    pickmap = Dict{UInt64, Any}()
    syncing = Ref(false)

    function draw_slot!(slot)
        foreach(delete!, slot.blocks); empty!(slot.blocks)
        trim!(slot.layout)
        kind = slot.kind[]
        before = length(fig.content)
        ax = kind == :ptstrip ? _mt_axis(slot.layout[2, 1]) : _db_period_axis(slot.layout[2, 1])
        add!(B, kw) = B(slot.layout[2, 2]; kw...)
        is = cur[]
        kind == :tipper  ? _db_tipper!(ax, st, is, pickmap) :
        kind == :ptstrip ? _db_ptstrip!(ax, st, is, add!) :
        kind == :beta    ? _db_beta!(ax, st, is) :
        kind == :skew    ? _db_skews!(ax, st, is) :
        kind == :strike  ? _db_strike!(ax, st, is) :
        kind == :resid   ? _db_resid!(ax, st, is) :
        kind == :relerr  ? _db_relerr!(ax, st, is) : _db_nb!(ax, st, is)
        append!(slot.blocks, fig.content[before+1:end])
        nothing
    end

    function refresh_map!()
        ip = map_period[]
        polys, cols = Vector{Point2f}[], Float32[]
        for is in 1:d.ns
            pt = fr.PT[ip, is]
            (pt === nothing || !all(st.keep[ip, 1:4, is])) && continue
            a, b = _ellipse_semiaxes(pt, 0.8, fr.Lref)
            a > 0 || continue
            xs, ys = _ellipse_ring(fr.site_x[is], fr.site_y[is], a, b, pt.azimuth + d.zrot[ip, is]; kx = fr.lon_stretch)
            push!(polys, Point2f.(xs, ys)); push!(cols, Float32(atand(pt.phimin)))
        end
        map_polys[], map_cols[] = polys, cols
        period_label.text[] = @sprintf("T = %.4g s", d.T[ip])
        if st.pred !== nothing
            r = [_db_site_rms(st, is).rms for is in 1:d.ns]
            fin = filter(isfinite, r)
            isempty(fin) || (rms_range[] = (0.0, max(1.0, quantile(fin, 0.95))))
            site_colour[] = r
        end
    end

    function refresh_info!()
        is = cur[]
        lat, lon, elev = d.loc[is, :]
        nT = count(ip -> any(isfinite, d.Z[ip, :, is]) || any(isfinite, d.tip[ip, :, is]), 1:d.nf)
        nkept, nall = count(st.keep[:, :, is]), count(isfinite, d.Z[:, :, is]) + count(isfinite, d.tip[:, :, is])
        s = @sprintf("%s  (%d/%d)   %.4f°, %.4f°, %.0f m   %d periods   kept %d/%d", d.site[is], is, d.ns, lat, lon, elev, nT, nkept, nall)
        if st.pred !== nothing
            s *= @sprintf("   RMS %.2f (all %.2f)", _db_site_rms(st, is).rms, _db_total_rms(st))
        end
        info.text[] = s
    end

    function redraw!()
        empty!(pickmap)
        _db_rhophi!(axρ, axφ, st, cur[], diag_toggle.active[], pickmap)
        foreach(draw_slot!, slots)
        refresh_info!()
    end

    function goto!(is::Integer)
        cur[] = clamp(is, 1, d.ns)
        current[] = Point2f(fr.site_x[cur[]], fr.site_y[cur[]])
        syncing[] = true
        menu_sites.i_selected[] = cur[]
        syncing[] = false
        redraw!()
        cur[]
    end

    on(b_first.clicks) do _; goto!(1); end
    on(b_prev.clicks) do _; goto!(cur[] == 1 ? d.ns : cur[] - 1); end
    on(b_next.clicks) do _; goto!(cur[] % d.ns + 1); end
    on(b_last.clicks) do _; goto!(d.ns); end
    on(menu_sites.selection) do i
        syncing[] || i === nothing || i == cur[] || goto!(i)
    end
    on(diag_toggle.active) do _; redraw!(); end
    for slot in slots
        on(slot.menu.selection) do kind
            (kind === nothing || kind == slot.kind[]) && return
            slot.kind[] = kind
            empty!(pickmap)
            _db_rhophi!(axρ, axφ, st, cur[], diag_toggle.active[], pickmap)
            foreach(draw_slot!, slots)
        end
    end
    on(pt_toggle.active) do v; map_pt[] = v; end
    on(period_slider.value) do ip; map_period[] = ip; refresh_map!(); end

    function set_site!(value::Bool)
        is = cur[]
        for ic in 1:6, ip in 1:d.nf
            fin = ic <= 4 ? isfinite(d.Z[ip, ic, is]) : isfinite(d.tip[ip, ic - 4, is])
            st.keep[ip, ic, is] = value && fin
        end
        refresh_map!(); redraw!()
    end
    on(b_site.clicks) do _; set_site!(false); status.text[] = "Dropped $(d.site[cur[]])"; end
    on(b_reset.clicks) do _; set_site!(true); status.text[] = "Restored $(d.site[cur[]])"; end
    on(b_export.clicks) do _
        try
            path = isempty(export_path) ? joinpath(pwd(), "$(stem)Edited.dat") : export_path
            redirect_stdout(devnull) do
                write_data_modem(path, _db_export_data(st); sign = 1, units = "[mV/km]/[nT]",
                                 include_tipper = any(st.keep[:, 5:6, :]),
                                 description = @sprintf("DataDashboard export of %s, error floors %.3g sqrt|ZxyZyx| and %.3g", stem, st.z_floor, st.t_floor))
            end
            status.text[] = "Wrote $path"
            println("Wrote $path")
        catch e
            status.text[] = "Export failed: $(sprint(showerror, e))"
        end
    end
    on(b_png.clicks) do _
        path = joinpath(pwd(), "$(stem)-$(d.site[cur[]]).png")
        try
            save(path, fig; px_per_unit = 2, backend = CairoMakie)
            status.text[] = "Saved $path"
        catch e
            status.text[] = "Save failed: $(sprint(showerror, e))"
        end
    end

    on(events(fig).mousebutton, priority = 10) do ev
        (ev.button == Mouse.left && ev.action == Mouse.press) || return Consume(false)
        if is_mouseinside(map_ax.scene)
            plt, i = pick(map_ax.scene, events(map_ax.scene).mouseposition[], 12)
            if plt === site_plot && i > 0
                goto!(i)
                return Consume(true)
            end
            return Consume(false)
        end
        edit_toggle.active[] || return Consume(false)
        axes_now = vcat([axρ, axφ], [s.blocks[1] for s in slots if s.kind[] == :tipper && !isempty(s.blocks)])
        for ax in axes_now
            is_mouseinside(ax.scene) || continue
            plt, i = pick(ax.scene, events(ax.scene).mouseposition[], 12)
            hit = plt === nothing ? nothing : get(pickmap, objectid(plt), nothing)
            (hit === nothing || i < 1 || i > length(hit[2])) && return Consume(false)
            ic, idx = hit
            ip = idx[i]
            st.keep[ip, ic, cur[]] = !st.keep[ip, ic, cur[]]
            status.text[] = @sprintf("%s %s T = %.4g s %s", d.site[cur[]], ic <= 4 ? "Z" * lowercase(_DB_ZNAMES[ic]) : (ic == 5 ? "Tzx" : "Tzy"),
                                     d.T[ip], st.keep[ip, ic, cur[]] ? "restored" : "masked")
            refresh_map!(); redraw!()
            return Consume(true)
        end
        Consume(false)
    end
    on(events(fig).keyboardbutton) do ev
        ev.action in (Keyboard.press, Keyboard.repeat) || return Consume(false)
        ev.key == Keyboard.right && (goto!(cur[] % d.ns + 1); return Consume(true))
        ev.key == Keyboard.left && (goto!(cur[] == 1 ? d.ns : cur[] - 1); return Consume(true))
        ev.key == Keyboard.up && (set_close_to!(period_slider, min(d.nf, map_period[] + 1)); return Consume(true))
        ev.key == Keyboard.down && (set_close_to!(period_slider, max(1, map_period[] - 1)); return Consume(true))
        Consume(false)
    end

    refresh_map!()
    goto!(1)
    (fig = fig, goto! = goto!, state = st, map_axis = map_ax, rho_axis = axρ, phase_axis = axφ, slots = slots)
end

"""
    DataDashboard(observed; predicted=nothing, z_floor=0.05, t_floor=0.03, panels=nothing,
                  iv_convention=:parkinson, export_path="", figsize=(1850, 1080),
                  interactive=true, block=!isinteractive(), snapshot_dir="", snapshot_sites=nothing)

Page through the sites of a ModEM data file (convert EDIs with `EDIToModEM` first).
`predicted` is a ModEM response, matched to the observed sites by name and to their
periods within 2 %, and turned into the observed frame when the rotations differ.

Each site shows apparent resistivity and phase (all four components, with error bars),
and four panels, each with its own chooser: tipper, phase tensor ellipses and
induction arrows along period, PT skew β, ellipticity with Swift and Bahr skew, PT and
Swift strike, normalised residuals, relative errors, and Niblett–Bostick depth.
Residuals and RMS use the errors floored at `z_floor`·√|Zxy·Zyx| and `t_floor`; the
error bars are the recorded ones. The map shows the sites, coloured by RMS, and the
phase tensors at one period.

Keys: ←/→ change site, ↑/↓ change the map period; a click on the map picks a site. With
"Edit mask" on, a click on a ρ, φ or tipper point masks or restores it; "Export ModEM"
writes the kept data with floored errors to `export_path` (default
`<name>Edited.dat` in the working directory).

`interactive = false` needs no display: it builds the window with CairoMakie and, when
`snapshot_dir` is given, writes one PNG per site (or per name in `snapshot_sites`).
Returns the window's NamedTuple (`fig`, `goto!`, `state`, ...).
"""
function DataDashboard(observed::AbstractString;
    predicted::Union{Nothing, AbstractString} = nothing,
    z_floor::Real = 0.05,
    t_floor::Real = 0.03,
    panels = nothing,
    iv_convention::Symbol = :parkinson,
    export_path::AbstractString = "",
    figsize = (1850, 1080),
    interactive::Bool = true,
    block::Bool = !isinteractive(),
    snapshot_dir::AbstractString = "",
    snapshot_sites = nothing)

    (isdir(observed) || occursin(r"\.edi$"i, observed)) &&
        error("DataDashboard reads ModEM data; convert the EDIs first: EDIToModEM(\"$observed\")")
    isfile(observed) || error("not found: $observed")
    obs = redirect_stdout(() -> load_data_modem(observed; warn_rotation = false), devnull)
    for A in (obs.Zerr, obs.tiperr), i in eachindex(A)
        abs(A[i]) > 1e10 && (A[i] = complex(NaN, NaN))
    end
    stem = replace(splitext(basename(observed))[1], r"[^A-Za-z0-9]" => "")
    println("Observed: $observed")
    @printf("  %d sites, %d periods (%.4g .. %.4g s), components %s\n", obs.ns, obs.nf, minimum(obs.T), maximum(obs.T), join(obs.responses, " "))
    pred = if predicted === nothing
        nothing
    else
        isfile(predicted) || error("predicted file not found: $predicted")
        println("Predicted: $predicted")
        _db_align_predicted(obs, redirect_stdout(() -> load_data_modem(predicted; warn_rotation = false), devnull))
    end

    ez, et = _error_floors(obs, z_floor, t_floor)
    keep = BitArray(undef, obs.nf, 6, obs.ns)
    keep[:, 1:4, :] .= isfinite.(obs.Z)
    keep[:, 5:6, :] .= isfinite.(obs.tip)
    Trange = (minimum(obs.T) / 1.5, maximum(obs.T) * 1.5)
    st = (obs = obs, pred = pred, ez = ez, et = et, keep = keep, Trange = Trange,
          z_floor = z_floor, t_floor = t_floor, convention = iv_convention)
    pred === nothing || @printf("  RMS %.3f over all sites (errors floored at %.3g·sqrt|ZxyZyx| and %.3g)\n", _db_total_rms(st), z_floor, t_floor)

    tipper = has_tipper_data(obs)
    panels = something(panels, pred === nothing ? [tipper ? :tipper : :nb, :ptstrip, :beta, :strike] :
                                                  [tipper ? :tipper : :nb, :ptstrip, :resid, :beta])
    length(panels) == 4 && all(p -> p in first.(_DB_PANELS), panels) ||
        throw(ArgumentError("panels must be four of $(first.(_DB_PANELS))"))

    if interactive
        isdefined(@__MODULE__, :GLMakie) || error("the dashboard window needs GLMakie and a display; interactive = false renders with CairoMakie")
        GLMakie.activate!(title = "MTGeophysics data dashboard: $(basename(observed))")
    else
        CairoMakie.activate!()
    end
    w = _dashboard_window(st, stem; panels = collect(panels), export_path = export_path, figsize = figsize)

    if !isempty(snapshot_dir)
        mkpath(snapshot_dir)
        which = snapshot_sites === nothing ? (1:obs.ns) : [findfirst(==(s), obs.site) for s in snapshot_sites]
        for is in which
            is === nothing && continue
            w.goto!(is)
            save(joinpath(snapshot_dir, "$(stem)-$(obs.site[is]).png"), w.fig; px_per_unit = 1.5, backend = CairoMakie)
        end
        println("  $(length(which)) snapshot(s) in $snapshot_dir")
        w.goto!(1)
    end
    if interactive
        screen = display(w.fig)
        println("Dashboard ready: ←/→ site, ↑/↓ map period, click the map to pick a site")
        block && wait(screen)
    end
    w
end

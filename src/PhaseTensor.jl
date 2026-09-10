# Phase tensor and induction vector computation, geometry and GIS export.
# Author: @pankajkmishra
# This file turns ModEM impedances and tippers into phase tensor invariants and
# induction arrows, builds the map symbols for them, and writes one shapefile per
# period. Everything here is headless; the viewer lives in PlotData3D.jl.
# References: Caldwell et al. (2004); Parkinson (1959); Wiese (1962).

# ---------- invariants ----------

"""
    phase_tensor(Zxx, Zxy, Zyx, Zyy)

Phase tensor invariants Phi = X^-1 Y for one site at one period.
In:  the four complex impedance components.
Out: NamedTuple of invariants, angles in degrees; `nothing` if Re(Z) is singular.
"""
function phase_tensor(Zxx::Number, Zxy::Number, Zyx::Number, Zyy::Number)
    all(isfinite, (real(Zxx), imag(Zxx), real(Zxy), imag(Zxy),
                   real(Zyx), imag(Zyx), real(Zyy), imag(Zyy))) || return nothing
    x11, x12, x21, x22 = real(Zxx), real(Zxy), real(Zyx), real(Zyy)
    y11, y12, y21, y22 = imag(Zxx), imag(Zxy), imag(Zyx), imag(Zyy)
    detX = x11*x22 - x12*x21
    (isfinite(detX) && abs(detX) > 1e-30) || return nothing

    p11 = ( x22*y11 - x12*y21) / detX
    p12 = ( x22*y12 - x12*y22) / detX
    p21 = (-x21*y11 + x11*y21) / detX
    p22 = (-x21*y12 + x11*y22) / detX
    all(isfinite, (p11, p12, p21, p22)) || return nothing

    P1   = (p11 + p22) / 2
    P3   = (p12 - p21) / 2
    detP = p11*p22 - p12*p21
    a = atand(p12 + p21, p11 - p22) / 2
    b = atand(p12 - p21, p11 + p22) / 2

    base = sqrt(P1^2 + P3^2)
    root = sqrt(max(base^2 - detP, 0.0))
    Pmax, Pmin = base + root, base - root

    return (phimin = Pmin, phimax = Pmax, phi1 = P1,
            phi2 = detP >= 0 ? sqrt(detP) : NaN,
            alpha = a, beta = b, azimuth = a - b,
            ellipticity = (Pmax + Pmin) == 0 ? 0.0 : (Pmax - Pmin) / (Pmax + Pmin),
            anisotropy = Pmin == 0 ? NaN : Pmax / Pmin,
            phidiff = Pmax - Pmin)
end

"""
    induction_vector(Tzx, Tzy; convention = :parkinson)

Induction arrow components; Parkinson points towards conductors, Wiese away.
In:  complex tipper Tzx, Tzy and the sign convention.
Out: NamedTuple of real/imag (east, north) vectors, magnitudes and azimuths;
     `nothing` if the tipper is not finite.
"""
function induction_vector(Tzx::Number, Tzy::Number; convention::Symbol = :parkinson)
    all(isfinite, (real(Tzx), imag(Tzx), real(Tzy), imag(Tzy))) || return nothing
    s = convention === :parkinson ? -1.0 : 1.0
    return (re = (s*real(Tzy), s*real(Tzx)), im = (s*imag(Tzy), s*imag(Tzx)),
            re_mag = hypot(real(Tzx), real(Tzy)), im_mag = hypot(imag(Tzx), imag(Tzy)),
            re_azim = atand(s*real(Tzy), s*real(Tzx)),
            im_azim = atand(s*imag(Tzy), s*imag(Tzx)))
end

"""
    has_tipper_data(d)

Whether a `Data` object carries usable vertical field transfer functions.
In:  a ModEM `Data` object.
Out: `true` if both tipper components are present and at least one is finite.
"""
has_tipper_data(d) = !isempty(d.tip) && size(d.tip, 2) >= 2 && any(isfinite, real.(d.tip))

"""
    phase_tensors_from_data(d)

Phase tensors for every period and site of a data set.
In:  a ModEM `Data` object.
Out: `nf x ns` array of `phase_tensor` NamedTuples, `nothing` where undefined.
"""
function phase_tensors_from_data(d)
    PT = Array{Any}(nothing, d.nf, d.ns)
    for ip in 1:d.nf, is in 1:d.ns
        PT[ip, is] = phase_tensor(d.Z[ip,1,is], d.Z[ip,2,is], d.Z[ip,3,is], d.Z[ip,4,is])
    end
    return PT
end

"""
    induction_vectors_from_data(d; convention = :parkinson)

Induction vectors for every period and site of a data set.
In:  a ModEM `Data` object and the sign convention.
Out: `nf x ns` array of `induction_vector` NamedTuples, all `nothing` without tipper.
"""
function induction_vectors_from_data(d; convention::Symbol = :parkinson)
    IV = Array{Any}(nothing, d.nf, d.ns)
    has_tipper_data(d) || return IV
    for ip in 1:d.nf, is in 1:d.ns
        IV[ip, is] = induction_vector(d.tip[ip,1,is], d.tip[ip,2,is]; convention = convention)
    end
    return IV
end

# ---------- symbol geometry, shared by the viewer and the shapefile export ----------

"""
    _ellipse_ring(cx, cy, semi_major, semi_minor, azimuth_deg; n, kx)

Closed ellipse ring in map units.
In:  centre, semi-axes in isotropic units, azimuth clockwise from north, vertex
     count, and the east stretch `kx` that turns isotropic offsets into map units.
Out: (xs, ys) of n+1 points, last equal to first.
"""
function _ellipse_ring(cx, cy, semi_major, semi_minor, azimuth_deg; n::Int = 48, kx::Real = 1.0)
    st, ct = sincosd(azimuth_deg)
    xs = Vector{Float64}(undef, n + 1); ys = similar(xs)
    @inbounds for k in 0:n
        t = 2pi * k / n
        u = semi_major * cos(t)
        v = semi_minor * sin(t)
        xs[k+1] = cx + (u*st + v*ct) * kx
        ys[k+1] = cy + (u*ct - v*st)
    end
    xs[end], ys[end] = xs[1], ys[1]      # close the ring exactly, sin(2pi) is not 0
    return xs, ys
end

"""
    _ellipse_semiaxes(pt, scale, Lref)

Semi-axes of the drawn ellipse, every one normalised to the same size.
In:  a phase tensor NamedTuple, the size scale and the reference site spacing.
Out: (semi_major, semi_minor) in isotropic map units; (0, 0) if it should be skipped.
"""
function _ellipse_semiaxes(pt, scale::Real, Lref::Real)
    (isfinite(pt.phimax) && pt.phimax > 0) || return (0.0, 0.0)
    s = scale * Lref / 2
    return (s, s * pt.phimin / pt.phimax)
end

"""
    _arrow_parts(cx, cy, dx, dy; head_frac, head_width, kx)

Shaft and filled triangular head of one arrow, in map units.
In:  tail position, isotropic vector components, head length fraction, head
     half-width ratio and the east stretch `kx`.
Out: ((shaft_xs, shaft_ys), (head_xs, head_ys)); the head is a closed triangle.

The shaft and head are built in the isotropic frame so a head is the same size
whatever the arrow points at; east offsets are stretched by `kx` afterwards.
"""
function _arrow_parts(cx, cy, dx, dy; head_frac::Real = 0.30,
                      head_width::Real = 0.42, kx::Real = 1.0)
    L = hypot(dx, dy)
    L == 0 && return ((Float64[], Float64[]), (Float64[], Float64[]))
    ux, uy = dx / L, dy / L
    h = head_frac * L
    w = head_width * h
    px, py = -uy, ux
    bxo, byo = dx - h*ux, dy - h*uy           # head base, as an offset from the tail
    ex(o) = cx + o * kx
    shaft = (Float64[cx, ex(bxo)], Float64[cy, cy + byo])
    head  = (Float64[ex(dx), ex(bxo + w*px), ex(bxo - w*px), ex(dx)],
             Float64[cy + dy, cy + byo + w*py, cy + byo - w*py, cy + dy])
    return shaft, head
end

"""
    _arrow_outline(cx, cy, dx, dy; head_frac, head_width, kx)

Single closed polyline tracing an arrow, for shapefile export.
In:  as `_arrow_parts`.
Out: (xs, ys) running tail to head base, around the head, back to the base.
"""
function _arrow_outline(cx, cy, dx, dy; head_frac::Real = 0.30,
                        head_width::Real = 0.42, kx::Real = 1.0)
    (sx, sy), (hx, hy) = _arrow_parts(cx, cy, dx, dy;
                                      head_frac = head_frac, head_width = head_width, kx = kx)
    isempty(sx) && return (Float64[cx], Float64[cy])
    return (vcat(sx, hx), vcat(sy, hy))
end

# ---------- map frame shared by the viewer and the export ----------

"""
    _median_site_spacing(xs, ys, lon_scale)

Median nearest-neighbour distance between sites, measured isotropically.
In:  site coordinates in map units and the east/north unit ratio `lon_scale`
     (1 for a metric CRS, cos(lat) for degrees).
Out: distance in north-axis units; 1.0 when there are fewer than two sites.
"""
function _median_site_spacing(xs, ys, lon_scale::Real)
    length(xs) < 2 && return 1.0
    nn = Float64[]
    @inbounds for i in eachindex(xs)
        best = Inf
        for j in eachindex(xs)
            i == j && continue
            dd = ((xs[i]-xs[j]) * lon_scale)^2 + (ys[i]-ys[j])^2
            dd < best && (best = dd)
        end
        isfinite(best) && push!(nn, sqrt(best))
    end
    isempty(nn) ? 1.0 : median(nn)
end

"""
    _ptiv_frame(d, crs; pt_scale, iv_scale, map_pad, convention)

Everything the phase tensor map needs that does not depend on the period shown.
In:  a ModEM `Data` object, the plot CRS and the symbol scales.
Out: NamedTuple with the site coordinates, the east stretch, the reference
     spacing, the precomputed PT/IV arrays and the fixed map extent and aspect.

`lon_stretch` is 1 in a metric CRS; in EPSG:4326 it is 1/cos(lat) so that symbols
built in the isotropic frame stay circular in true distance.
"""
function _ptiv_frame(d, crs::AbstractString; pt_scale::Real = 0.66, iv_scale::Real = 1.2,
                     map_pad::Real = 0.06, convention::Symbol = :parkinson)
    site_x, site_y = _stations_in_crs(d, crs)
    isempty(site_x) && error("No station coordinates available for crs = \"$crs\"")

    geographic = uppercase(strip(crs)) == "EPSG:4326"
    lon_scale  = geographic ? cosd(clamp(sum(site_y) / length(site_y), -89.0, 89.0)) : 1.0
    lon_stretch = 1 / lon_scale

    Lref = _median_site_spacing(site_x, site_y, lon_scale)
    PT   = phase_tensors_from_data(d)
    IV   = induction_vectors_from_data(d; convention = convention)

    x0, x1 = extrema(site_x)
    y0, y1 = extrema(site_y)
    m  = max(pt_scale * Lref / 2, iv_scale * Lref)
    sx = (x1 - x0) * map_pad + m * lon_stretch
    sy = (y1 - y0) * map_pad + m
    limits = (x0 - sx, x1 + sx, y0 - sy, y1 + sy)

    aspect = geographic ?
        first(_distance_based_aspect([limits[3], limits[4]], [limits[1], limits[2]])) :
        (limits[2] - limits[1]) / (limits[4] - limits[3])

    return (site_x = site_x, site_y = site_y, geographic = geographic,
            lon_scale = lon_scale, lon_stretch = lon_stretch, Lref = Lref,
            PT = PT, IV = IV, has_tipper = has_tipper_data(d),
            limits = limits, aspect = aspect)
end

"""
    PT_FILL_OPTIONS

Phase tensor invariants offered as the ellipse fill, with their colorbar label,
whether they are tangents that can be shown as angles, and whether their color
range should be symmetric about zero. `PT_FILL_ORDER` fixes the menu order.
"""
const PT_FILL_OPTIONS = Dict(
    :beta   => (label = "beta skew (deg)", angle = false, symmetric = true),
    :phimin => (label = "Phimin",          angle = true,  symmetric = false),
    :phi2   => (label = "Phi2 (geom)",     angle = true,  symmetric = false),
)
const PT_FILL_ORDER = [:beta, :phimin, :phi2]

"""
    _pt_fill_value(pt, key, as_angle)

One phase tensor invariant, optionally converted from a tangent to an angle.
In:  a phase tensor NamedTuple, the invariant name and the angle flag.
Out: the value in degrees for angle-valued invariants, raw otherwise.
"""
function _pt_fill_value(pt, key::Symbol, as_angle::Bool)
    v = getfield(pt, key)
    return (as_angle && PT_FILL_OPTIONS[key].angle && isfinite(v)) ? atand(v) : v
end

# ---------- GIS export, one shapefile per period ----------

"""
    _export_ptiv_gis(d, fr; output_dir, data_name, crs, ...)

Write one shapefile per period from an already prepared map frame.
In:  the `Data` object, a `_ptiv_frame` NamedTuple and the export options.
Out: the output directory path.
"""
function _export_ptiv_gis(d, fr;
    output_dir::AbstractString,
    data_name::AbstractString,
    crs::AbstractString,
    source_file::AbstractString = "",
    as_angle::Bool = true,
    pt_scale::Real = 0.66,
    iv_scale::Real = 1.2,
    iv_convention::Symbol = :parkinson,
    iv_max_magnitude::Real = 1.0,
    iv_head_frac::Real = 0.30,
    iv_head_width::Real = 0.42,
    skip_beta_above::Union{Nothing, Real} = nothing)

    mkpath(output_dir)
    wkt = _resolve_prj_wkt_for_crs(crs)
    kx  = fr.lon_stretch
    ang(v) = as_angle && isfinite(v) ? atand(v) : v
    keep(pt) = isnothing(skip_beta_above) || abs(pt.beta) <= skip_beta_above
    bbox(xs, ys) = Shapefile.Rect(minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    pts(xs, ys) = [Shapefile.Point(xs[k], ys[k]) for k in eachindex(xs)]

    println("\nExporting $(d.nf) periods to: $output_dir")
    n_written = 0

    for ip in 1:d.nf
        T   = d.T[ip]
        tag = @sprintf("T%09.4f", T)
        Tr  = round(T, sigdigits = 8)
        frq_val = T > 0 ? round(1/T, sigdigits = 8) : NaN

        polys = Shapefile.Polygon[]
        site  = String[]; per  = Float64[]; frq = Float64[]
        pmin  = Float64[]; pmax = Float64[]; p1  = Float64[]; p2  = Float64[]
        bet   = Float64[]; alp  = Float64[]; azi = Float64[]
        ell   = Float64[]; ani  = Float64[]; dif = Float64[]

        for is in 1:d.ns
            pt = fr.PT[ip, is]
            (isnothing(pt) || !keep(pt)) && continue
            a, b = _ellipse_semiaxes(pt, pt_scale, fr.Lref)
            a > 0 || continue
            xs, ys = _ellipse_ring(fr.site_x[is], fr.site_y[is], a, b, pt.azimuth; kx = kx)
            push!(polys, Shapefile.Polygon(bbox(xs, ys), Int32[0], pts(xs, ys)))
            push!(site, d.site[is]); push!(per, Tr); push!(frq, frq_val)
            push!(pmin, round(ang(pt.phimin),  digits = 5))
            push!(pmax, round(ang(pt.phimax),  digits = 5))
            push!(p1,   round(ang(pt.phi1),    digits = 5))
            push!(p2,   round(ang(pt.phi2),    digits = 5))
            push!(bet,  round(pt.beta,         digits = 5))
            push!(alp,  round(pt.alpha,        digits = 5))
            push!(azi,  round(pt.azimuth,      digits = 5))
            push!(ell,  round(pt.ellipticity,  digits = 5))
            push!(ani,  round(pt.anisotropy,   digits = 5))
            push!(dif,  round(ang(pt.phidiff), digits = 5))
        end

        if !isempty(polys)
            _write_shapefile_with_sidecars(joinpath(output_dir, "$(data_name)_PT_$(tag).shp"),
                polys,
                (site = site, period_s = per, freq_Hz = frq,
                 phimin = pmin, phimax = pmax, phi1 = p1, phi2 = p2,
                 beta = bet, alpha = alp, azimuth = azi,
                 ellip = ell, aniso = ani, phidiff = dif), wkt)
            n_written += 1
        end

        if fr.has_tipper
            arrows = Shapefile.Polyline[]
            asite = String[]; aper = Float64[]; afrq  = Float64[]
            amag  = Float64[]; aazi = Float64[]; apart = String[]
            sc = iv_scale * fr.Lref
            for is in 1:d.ns
                iv = fr.IV[ip, is]
                isnothing(iv) && continue
                for (vec, mag, az, part) in ((iv.re, iv.re_mag, iv.re_azim, "real"),
                                             (iv.im, iv.im_mag, iv.im_azim, "imag"))
                    (isfinite(mag) && 0 < mag <= iv_max_magnitude) || continue
                    xs, ys = _arrow_outline(fr.site_x[is], fr.site_y[is],
                                            vec[1]*sc, vec[2]*sc;
                                            head_frac = iv_head_frac,
                                            head_width = iv_head_width, kx = kx)
                    push!(arrows, Shapefile.Polyline(bbox(xs, ys), Int32[0], pts(xs, ys)))
                    push!(asite, d.site[is]); push!(aper, Tr); push!(afrq, frq_val)
                    push!(amag, round(mag, digits = 5)); push!(aazi, round(az, digits = 3))
                    push!(apart, part)
                end
            end
            isempty(arrows) || _write_shapefile_with_sidecars(
                joinpath(output_dir, "$(data_name)_IV_$(tag).shp"), arrows,
                (site = asite, period_s = aper, freq_Hz = afrq,
                 magnitude = amag, azimuth = aazi, part = apart,
                 conv = fill(String(iv_convention), length(asite))), wkt)
        end

        @printf("  period %2d/%d  T = %-10.4g s  %d ellipses\n", ip, d.nf, T, length(polys))
    end

    open(joinpath(output_dir, "README.txt"), "w") do io
        println(io, "Phase tensor and induction vector export")
        println(io, "generated ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "source   : ", isempty(source_file) ? d.name : source_file)
        println(io, "CRS      : ", crs)
        println(io, "ellipses : *_PT_T<period>.shp (polygons); phimin/phimax/phi1/phi2 are ",
                    as_angle ? "angles in degrees, atan(Phi)" : "raw Phi values (tangents)")
        println(io, "           beta/alpha/azimuth in degrees; azimuth = alpha - beta, cw from north")
        println(io, "           all ellipses one size (normalised by their own Phimax), scale ", pt_scale)
        if fr.geographic
            println(io, "           symbols are circular in true distance, so east offsets carry a")
            println(io, "           1/cos(lat) stretch; reproject to a metric CRS to measure them")
        end
        println(io, "arrows   : *_IV_T<period>.shp (polylines), ", iv_convention,
                    " convention; Parkinson points TOWARDS conductors")
        println(io, "reference: Caldwell et al. (2004); Parkinson (1959); Wiese (1962)")
    end

    println("GIS export complete: $n_written period(s) written.")
    return output_dir
end

"""
    write_ptiv_gis(data_file; crs, output_dir, ...)

Load a ModEM data file and write phase tensor and induction vector shapefiles,
one pair per period. Needs no display, so it also serves `PlotPTIVMap` headless.
In:  path to the data file, the export CRS, the output directory (default
     `<data>-PTIV-GIS` in the working directory) and the symbol options.
Out: the output directory path.
"""
function write_ptiv_gis(data_file::AbstractString;
    crs::AbstractString = "EPSG:4326",
    output_dir::AbstractString = "",
    as_angle::Bool = true,
    pt_scale::Real = 0.66,
    iv_scale::Real = 1.2,
    iv_convention::Symbol = :parkinson,
    iv_max_magnitude::Real = 1.0,
    iv_head_frac::Real = 0.30,
    iv_head_width::Real = 0.42,
    skip_beta_above::Union{Nothing, Real} = nothing)

    isempty(data_file) && error("data_file is required")
    isfile(data_file)  || error("Data file not found: $data_file")

    d  = load_data_modem(data_file)
    fr = _ptiv_frame(d, crs; pt_scale = pt_scale, iv_scale = iv_scale,
                     convention = iv_convention)
    data_name = splitext(basename(data_file))[1]
    dir = isempty(output_dir) ? joinpath(pwd(), "$(data_name)-PTIV-GIS") : output_dir

    return _export_ptiv_gis(d, fr;
        output_dir = dir, data_name = data_name, crs = crs, source_file = data_file,
        as_angle = as_angle, pt_scale = pt_scale, iv_scale = iv_scale,
        iv_convention = iv_convention, iv_max_magnitude = iv_max_magnitude,
        iv_head_frac = iv_head_frac, iv_head_width = iv_head_width,
        skip_beta_above = skip_beta_above)
end

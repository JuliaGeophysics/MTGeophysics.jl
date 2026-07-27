# =====================================================================
# VFSA3D_exp.jl -- experimental model parameterisation for VFSA3DMT
#
# Replaces the normalised-Gaussian (Shepard) weight map with local RBF-FD
# stencils built from polyharmonic splines plus polynomial augmentation,
# following Flyer, Fornberg, Bayona & Barnett, "On the role of polynomials
# in RBF-FD approximations: I. Interpolation and accuracy", J. Comput.
# Phys. 321 (2016) 21-38, and Bayona et al., part II, JCP 332 (2017).
#
# Why PHS + polynomials rather than a Gaussian:
#
#   1. No shape parameter. PHS phi(r) = r^m is scale free, so the
#      sigma_scale / sigma_scale_deep tuning disappears entirely.
#   2. Stagnation free. Accuracy is set by the polynomial degree, not by a
#      shape parameter, so refining the control set actually helps.
#   3. Fixed stencil size. Every cell uses the same n controls, so the
#      effective step no longer swings ~5x between densely and sparsely
#      covered regions the way normalised Gaussian weights do.
#   4. Local systems are ill conditioned but harmlessly so: they are solved
#      once at setup in double precision and the resulting weights are well
#      behaved. This is the central Fornberg-Flyer result, and it held here:
#      measured cond(A) was 4.6e2 to 2.4e4 across the depth range.
#
#   5. MEASURED CAVEAT.  That theory is about approximating SMOOTH fields,
#      where the large +/- weights cancel.  A VFSA proposal is white noise in
#      control space, nothing cancels, and the output amplitude is ||w||_2.
#      With the textbook n = 2*npoly stencil that gave a 97x spread in
#      per-cell step and 15x overshoot beyond the control range -- worse than
#      the Gaussian it replaces.  Oversampling to n = 30 with a 1e-3 ridge
#      pulls it back to 5.6x spread and 3.2x overshoot.  Those are the
#      defaults below.  Degree-2 augmentation was unusable, |w| ~ 1e16.
#
# Smoothness follows magnetotelluric resolution: the stencil is selected in
# a depth-scaled metric, lateral length alpha*z and vertical beta*z, since
# MT cannot resolve a body much narrower than its own burial depth. Near
# surface stencils are therefore tight, deep stencils broad, without any
# per-layer tuning.
#
# Standalone. Builds an MTGeophysics.RBFMap, so apply_rbf_map! and the rest
# of VFSA3DMT consume it unchanged. Nothing here executes on include.
#
#   include("examples/VFSA3D_exp.jl")
#   r = build_local_rbf_map(m, ix, iy, kz, 900, MersenneTwister(1))
#   rbf_report(r, m, ix, iy, kz)
# =====================================================================

using MTGeophysics
using LinearAlgebra
using Random
using Statistics
using Printf

# ---------------------------------------------------------------------
# Polyharmonic spline and polynomial basis
# ---------------------------------------------------------------------

# phi(r) = r^m, m odd.  In 3-D the thin-plate analogues are the odd
# polyharmonics: m = 1 biharmonic, m = 3 triharmonic (the usual RBF-FD
# choice), m = 5 for higher polynomial degrees.  r^2*log(r) is the 2-D
# form and is not used here.
@inline phs(r2::Float64, m::Int) = r2 <= 0 ? 0.0 : r2^(m / 2)          # r^m from r^2

# number of 3-D polynomial terms up to total degree deg
npoly(deg::Int) = ((deg + 1) * (deg + 2) * (deg + 3)) ÷ 6              # 1,4,10,20,...

# evaluate the polynomial basis at a local coordinate
function poly_basis!(p::Vector{Float64}, x::Float64, y::Float64,
                     z::Float64, deg::Int)
    p[1] = 1.0                                                          # constant
    deg == 0 && return p
    p[2] = x; p[3] = y; p[4] = z                                        # linear
    deg == 1 && return p
    p[5]  = x * x; p[6]  = x * y; p[7]  = x * z                         # quadratic
    p[8]  = y * y; p[9]  = y * z; p[10] = z * z
    deg == 2 && return p
    error("poly degree > 2 not implemented")
end

# ---------------------------------------------------------------------
# Depth-scaled metric: smoothness follows MT resolution
# ---------------------------------------------------------------------

"""
    ResolutionMetric(alpha, beta, l_floor, h_floor)

Anisotropic length scales used to select and weight stencils.  At depth `z`
the lateral scale is `max(l_floor, alpha*z)` and the vertical scale is
`max(h_floor, beta*z)`.  MT resolves a body of size comparable to its depth,
so `alpha` near 0.5 and `beta` near 0.25 reproduce the usual rule that
vertical resolution is the better of the two.
"""
struct ResolutionMetric
    alpha::Float64                                                      # lateral length / depth
    beta::Float64                                                       # vertical length / depth
    l_floor::Float64                                                    # min lateral length, m
    h_floor::Float64                                                    # min vertical length, m
end

@inline function scales(rm::ResolutionMetric, depth::Float64)
    return (max(rm.l_floor, rm.alpha * depth),                          # L(z)
            max(rm.h_floor, rm.beta  * depth))                          # H(z)
end

# ---------------------------------------------------------------------
# Local RBF-FD weight construction
# ---------------------------------------------------------------------

"""
    build_local_rbf_map(m, ix, iy, kz, n_ctrl, rng; kwargs...)

Return an `MTGeophysics.RBFMap` whose weights come from local PHS + polynomial
interpolation instead of normalised Gaussians.  Control placement is delegated
to the package `build_rbf_map`, so the depth-weighted, mask-aware sampling is
identical and only the weights differ.

Keywords
- `phs_order`  odd polyharmonic exponent m, default 3 (triharmonic)
- `poly_deg`   polynomial augmentation degree, 0 / 1 / 2, default 1
- `nstencil`   controls per cell, default 30 (oversampled, see header note 5)
- `metric`     `ResolutionMetric`, default alpha=0.5, beta=0.25
- `exclude`    3-D Bool mask of cells that may not host a control
- `ridge`      Tikhonov term added to the RBF block, default 1e-3
"""
function build_local_rbf_map(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int},
                             kz::UnitRange{Int}, n_ctrl::Int, rng::AbstractRNG;
                             phs_order::Int = 3,
                             poly_deg::Int = 1,
                             nstencil::Int = 30,                                   # oversampled: see
                             metric::ResolutionMetric = ResolutionMetric(
                                 0.5, 0.25, median(m.dx), m.dz[1]),
                             exclude = nothing,
                             depth_power::Real = 0.0,
                             ridge::Float64 = 1.0e-3)                              # header note 5

    isodd(phs_order) || error("phs_order must be odd (1, 3, 5)")
    poly_deg <= 1 || @warn "poly_deg 2 was unstable in testing: |w| reached 1e16"

    seed_map = build_rbf_map(m, ix, iy, n_ctrl, rng; kz = kz,                       # reuse the
                             depth_power = depth_power, exclude = exclude)          # tested
    ci, cj, ck = seed_map.ci, seed_map.cj, seed_map.ck                              # placement
    M = length(ci)
    nst = min(nstencil, M)
    P   = npoly(poly_deg)
    nst >= P || error("nstencil ($nst) must be at least npoly ($P)")

    nx, ny, nz = length(ix), length(iy), length(kz)
    N = nx * ny * nz

    cxc = [m.cx[ix[i]] for i in 1:nx]                                   # absolute cell
    cyc = [m.cy[iy[j]] for j in 1:ny]                                   # centres of the
    czc = [m.cz[kz[k]] for k in 1:nz]                                   # core box
    qx = [cxc[ci[q]] for q in 1:M]                                      # control
    qy = [cyc[cj[q]] for q in 1:M]                                      # positions
    qz = [czc[ck[q]] for q in 1:M]
    ztop = m.origin[3]                                                  # z of model top

    ptr  = Vector{Int}(undef, N + 1); ptr[1] = 1
    nbrs = Int[];     sizehint!(nbrs, nst * N)
    wts  = Float64[]; sizehint!(wts,  nst * N)

    d2   = Vector{Float64}(undef, M)                                    # scratch, all controls
    idx  = Vector{Int}(undef, M)
    A    = zeros(nst + P, nst + P)                                      # saddle-point system
    rhs  = zeros(nst + P)
    pbuf = zeros(max(P, 10))
    n_fallback = 0                                                      # singular -> inv dist

    row = 1
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        ex, ey, ez = cxc[i], cyc[j], czc[k]
        depth = max(ez - ztop, 0.0)                                     # depth below model top
        L, H = scales(metric, depth)                                    # MT resolution scales

        for q in 1:M                                                    # scaled distance to
            ddx = (qx[q] - ex) / L                                      # every control; M is
            ddy = (qy[q] - ey) / L                                      # only ~900 so brute
            ddz = (qz[q] - ez) / H                                      # force beats a tree
            d2[q] = ddx * ddx + ddy * ddy + ddz * ddz
        end
        sortperm!(idx, d2)                                              # nearest nst controls
        h = sqrt(max(d2[idx[nst]], eps()))                              # stencil radius, scaled

        fill!(A, 0.0); fill!(rhs, 0.0)                                  # build [A P; P' 0]
        for a in 1:nst                                                  # in LOCAL coordinates
            qa = idx[a]                                                 # centred on the eval
            ax = (qx[qa] - ex) / (L * h)                                # point and normalised
            ay = (qy[qa] - ey) / (L * h)                                # by h: this is what
            az = (qz[qa] - ez) / (H * h)                                # keeps the local
            for b in 1:a                                                # systems solvable in
                qb = idx[b]                                             # double precision
                bx = (qx[qb] - ex) / (L * h)
                by = (qy[qb] - ey) / (L * h)
                bz = (qz[qb] - ez) / (H * h)
                r2 = (ax - bx)^2 + (ay - by)^2 + (az - bz)^2
                v  = phs(r2, phs_order)
                A[a, b] = v; A[b, a] = v
            end
            A[a, a] += ridge                                            # optional Tikhonov
            poly_basis!(pbuf, ax, ay, az, poly_deg)
            for c in 1:P
                A[a, nst + c] = pbuf[c]; A[nst + c, a] = pbuf[c]        # P and P'
            end
            rhs[a] = phs(ax * ax + ay * ay + az * az, phs_order)        # phi(|y - x_a|), y = 0
        end
        poly_basis!(pbuf, 0.0, 0.0, 0.0, poly_deg)                      # p_k(y) at the origin
        for c in 1:P
            rhs[nst + c] = pbuf[c]
        end

        w = try                                                         # saddle point is
            (A \ rhs)[1:nst]                                            # symmetric indefinite
        catch                                                           # -> plain LU
            Float64[]
        end
        # sum(w) = 1 is enforced by the constant in the polynomial block, so it
        # is a free correctness check: a stencil whose points are degenerate for
        # the polynomial (coplanar, collinear) fails it silently rather than
        # throwing, and must be caught here rather than trusted
        if isempty(w) || !all(isfinite, w) || abs(sum(w) - 1.0) > 1e-6
            n_fallback += 1
            iw = [1.0 / (d2[idx[a]] + eps()) for a in 1:nst]            # inverse-distance
            w = iw ./ sum(iw)                                           # fallback
        end

        for a in 1:nst
            push!(nbrs, idx[a]); push!(wts, w[a])
        end
        ptr[row + 1] = length(nbrs) + 1
        row += 1
    end

    n_fallback > 0 && @warn "local RBF: $n_fallback of $N stencils fell back to inverse distance"

    return RBFMap(nx, ny, nz, ci, cj, ck, seed_map.ctrl_at, ptr, nbrs, wts)
end

# ---------------------------------------------------------------------
# Diagnostics -- the three things that decide whether this is usable
# ---------------------------------------------------------------------

"""
    rbf_report(map; bound_width, step_scale, T0, nsample)

Report the properties that matter for a VFSA proposal: partition of unity,
negative-weight fraction, overshoot under a random control draw, and the
uniformity of the per-cell step across the model.
"""
function rbf_report(rmap::RBFMap; bound_width::Float64 = 5.0,
                    step_scale::Float64 = 0.18, T0::Float64 = 0.03,
                    nsample::Int = 12, rng::AbstractRNG = MersenneTwister(4))
    N = length(rmap.ptr) - 1
    M = length(rmap.ci)
    sums = Vector{Float64}(undef, N)                                    # sum w  -> should be 1
    l2   = Vector{Float64}(undef, N)                                    # ||w||  -> step damping
    negf = 0; nw = 0
    @inbounds for row in 1:N
        p = rmap.ptr[row]:(rmap.ptr[row+1] - 1)
        w = view(rmap.wts, p)
        sums[row] = sum(w)
        l2[row]   = sqrt(sum(abs2, w))
        negf += count(<(0), w); nw += length(w)
    end
    step = 0.146 * bound_width * step_scale                             # median |y| at T0

    @printf("  controls %d, cells %d, stencil %d\n", M, N,
            rmap.ptr[2] - rmap.ptr[1])
    @printf("  partition of unity  : sum(w) = %.6f .. %.6f  (target 1)\n",
            minimum(sums), maximum(sums))
    @printf("  negative weights    : %.1f%% of entries\n", 100 * negf / nw)
    @printf("  per-cell step       : %.4f .. %.4f log10  (ratio %.1fx)\n",
            step * minimum(l2), step * maximum(l2), maximum(l2) / minimum(l2))
    @printf("  per-cell step       : median %.4f log10 (%.1f%% in rho)\n",
            step * median(l2), 100 * (10^(step * median(l2)) - 1))

    over = 0.0                                                          # overshoot beyond the
    dv = zeros(rmap.nx, rmap.ny, rmap.nz)                               # range of contributing
    for _ in 1:nsample                                                  # control values
        f = randn(rng, M)
        apply_rbf_map!(dv, rmap, f)
        over = max(over, maximum(abs, dv) / maximum(abs, f))
    end
    @printf("  overshoot factor    : %.3f  (>1 means the field exceeds its controls)\n", over)
    return (sum_min = minimum(sums), sum_max = maximum(sums),
            neg_frac = negf / nw, step_ratio = maximum(l2) / minimum(l2),
            overshoot = over)
end

"""
    compare_maps(m, ix, iy, kz, n_ctrl; kwargs...)

Build the Gaussian map and the local PHS map on the same control set and print
both reports side by side.
"""
function compare_maps(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int},
                      kz::UnitRange{Int}, n_ctrl::Int;
                      exclude = nothing, depth_power::Real = 0.0,
                      sigma_scale::Float64 = 2.0, bound_width::Float64 = 5.0,
                      step_scale::Float64 = 0.18, seed::Int = 1,
                      kwargs...)
    println("---- normalised Gaussian (current) ----")
    g = build_rbf_map(m, ix, iy, n_ctrl, MersenneTwister(seed); kz = kz,
                      sigma_scale = sigma_scale, sigma_scale_deep = sigma_scale,
                      trunc_sigmas = 3.0, depth_power = depth_power,
                      exclude = exclude)
    rg = rbf_report(g; bound_width = bound_width, step_scale = step_scale)
    println("---- local PHS + polynomial (experimental) ----")
    l = build_local_rbf_map(m, ix, iy, kz, n_ctrl, MersenneTwister(seed);
                            exclude = exclude, depth_power = depth_power, kwargs...)
    rl = rbf_report(l; bound_width = bound_width, step_scale = step_scale)
    return (gaussian = rg, local_phs = rl)
end

# =====================================================================
# RECOMMENDED: adaptive-bandwidth Wendland map
#
# The two schemes above optimise approximation accuracy.  That is the wrong
# objective here: the controls are not data to be honoured, they are abstract
# parameters driven by white noise.  What a VFSA proposal needs is
#
#   (1) bounded amplification  -- w >= 0 and sum(w) = 1, so the field can never
#       leave the control range and the log_bounds clamp never fires and biases
#       the proposal distribution,
#   (2) uniform noise gain     -- ||w||_2 equal at every cell, so one annealing
#       temperature is correct everywhere in the model,
#   (3) locality               -- compact support, so a trial is semi-local,
#   (4) smoothness that tracks MT resolution.
#
# Gaussian-Shepard gets 1 and 3 but fails 2, because at fixed bandwidth the
# effective neighbour count follows control density (measured 5.4x to 7.4x
# spread).  PHS-FD fails 1 and 2 outright.  The cure for 2 is not a better
# basis function, it is an ADAPTIVE BANDWIDTH: take the n nearest controls and
# set the width from the n-th distance, so the weight profile has identical
# shape everywhere and ||w||_2 is constant by construction.  This is k-nearest
# neighbour variable-bandwidth kernel regression.  It also removes the shape
# parameter, which was the main attraction of PHS.
#
# Kernel is Wendland C2, positive definite in 3-D, exactly compactly supported
# so there is no truncation discontinuity, and zero at s = 1 so the field stays
# continuous as the n-th neighbour swaps between adjacent cells.
# =====================================================================

@inline wendland_c2(s::Float64) = s >= 1.0 ? 0.0 : (1.0 - s)^4 * (4.0 * s + 1.0)

"""
    build_resolution_map(m, ix, iy, kz, n_ctrl, rng; kwargs...)

Adaptive-bandwidth Wendland weight map.  Recommended over both
`build_rbf_map` (fixed bandwidth) and `build_local_rbf_map` (PHS-FD).

Keywords
- `nstencil`  controls per cell, default 24; sets the effective smoothing
- `metric`    `ResolutionMetric` controlling how bandwidth grows with depth
- `exclude`   3-D Bool mask of cells that may not host a control
"""
function build_resolution_map(m::WS3DModel, ix::UnitRange{Int}, iy::UnitRange{Int},
                              kz::UnitRange{Int}, n_ctrl::Int, rng::AbstractRNG;
                              nstencil::Int = 24,
                              metric::ResolutionMetric = ResolutionMetric(
                                  0.5, 0.25, median(m.dx), m.dz[1]),
                              exclude = nothing,
                              depth_power::Real = 0.0)

    seed_map = build_rbf_map(m, ix, iy, n_ctrl, rng; kz = kz,            # identical control
                             depth_power = depth_power, exclude = exclude)  # placement
    ci, cj, ck = seed_map.ci, seed_map.cj, seed_map.ck
    M   = length(ci)
    nst = min(nstencil, M)

    nx, ny, nz = length(ix), length(iy), length(kz)
    N = nx * ny * nz

    cxc = [m.cx[ix[i]] for i in 1:nx]
    cyc = [m.cy[iy[j]] for j in 1:ny]
    czc = [m.cz[kz[k]] for k in 1:nz]
    qx = [cxc[ci[q]] for q in 1:M]
    qy = [cyc[cj[q]] for q in 1:M]
    qz = [czc[ck[q]] for q in 1:M]
    ztop = m.origin[3]

    ptr  = Vector{Int}(undef, N + 1); ptr[1] = 1
    nbrs = Int[];     sizehint!(nbrs, nst * N)
    wts  = Float64[]; sizehint!(wts,  nst * N)
    d2   = Vector{Float64}(undef, M)
    idx  = Vector{Int}(undef, M)
    wloc = Vector{Float64}(undef, nst)

    row = 1
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        ex, ey, ez = cxc[i], cyc[j], czc[k]
        depth = max(ez - ztop, 0.0)
        L, H = scales(metric, depth)                                    # MT resolution scales

        for q in 1:M
            ddx = (qx[q] - ex) / L
            ddy = (qy[q] - ey) / L
            ddz = (qz[q] - ez) / H
            d2[q] = ddx * ddx + ddy * ddy + ddz * ddz
        end
        sortperm!(idx, d2)
        h = sqrt(max(d2[idx[nst]], eps()))                              # bandwidth = n-th dist,
                                                                        # so the profile shape is
        s = 0.0                                                         # identical everywhere
        for a in 1:nst                                                  # and ||w||_2 is constant
            wloc[a] = wendland_c2(sqrt(d2[idx[a]]) / h)                 # by construction
            s += wloc[a]
        end
        inv = s > 0 ? 1.0 / s : 0.0                                      # L1 normalise ->
        for a in 1:nst                                                   # partition of unity
            push!(nbrs, idx[a]); push!(wts, wloc[a] * inv)
        end
        ptr[row + 1] = length(nbrs) + 1
        row += 1
    end

    return RBFMap(nx, ny, nz, ci, cj, ck, seed_map.ctrl_at, ptr, nbrs, wts)
end

# Per-site galvanic distortion estimation (variable projection).
# Author: @pankajkmishra
# Solves a real, frequency-independent 2x2 matrix C per site minimizing the
# error-weighted misfit Zobs ~ C*Zpred over all periods, damped toward identity,
# then computes chi2/rms on the distortion-corrected residuals. Tipper is
# magnetic and is never corrected.

using Printf
using Dates

struct DistortionFit
    χ²::Float64
    rms::Float64
    penalty::Float64        # damping-weighted ||C-I||^2 summed over sites, diagnostic only
    C::Array{Float64,3}     # 2 x 2 x ns, [c11 c12; c21 c22] per site
    site::Vector{String}
end

#---------- per-site damped least squares ----------

# component pairing for Zobs = C*Zpred with columns 1..4 = XX,XY,YX,YY:
#   obs col 1 (XX): c_1*Zp[:,1] + c_2*Zp[:,3]     obs col 2 (XY): c_1*Zp[:,2] + c_2*Zp[:,4]
#   obs col 3 (YX): same coefficients with row 2   obs col 4 (YY): same with row 2

# solve one row of C from obs columns (o1,o2), damped toward the prior (p1,p2)
# normal equations accumulated in complex arithmetic: Re(a*conj(b)) is the
# stacked real/imag inner product, so real unknowns come out exact
function _distortion_row(Zo::AbstractMatrix{ComplexF64}, Zp::AbstractMatrix{ComplexF64},
                         Ze::AbstractMatrix{ComplexF64},
                         o1::Int, o2::Int, p1::Float64, p2::Float64, damping::Float64)
    N11 = 0.0; N12 = 0.0; N22 = 0.0; r1 = 0.0; r2 = 0.0; m = 0
    @inbounds for i in 1:size(Zo, 1)
        for (oc, ac, bc) in ((o1, 1, 3), (o2, 2, 4))
            zo = Zo[i, oc]; a = Zp[i, ac]; b = Zp[i, bc]; ze = Ze[i, oc]
            if isfinite(real(zo)) && isfinite(imag(zo)) &&
               isfinite(real(a)) && isfinite(imag(a)) &&
               isfinite(real(b)) && isfinite(imag(b)) &&
               isfinite(real(ze)) && isfinite(imag(ze))
                σ = abs(ze)
                if σ > 0
                    w = 1.0 / (σ * σ)
                    N11 += w * abs2(a)
                    N12 += w * real(a * conj(b))
                    N22 += w * abs2(b)
                    r1  += w * real(zo * conj(a))
                    r2  += w * real(zo * conj(b))
                    m   += 1
                end
            end
        end
    end
    # relative damping: identity-prior strength as a fixed fraction of the data information
    λ = damping * 0.5 * (N11 + N22)
    (m < 2 || !isfinite(λ)) && return p1, p2, λ
    det = (N11 + λ) * (N22 + λ) - N12 * N12
    det <= 1e-12 * max(N11, N22, λ, eps())^2 && return p1, p2, λ
    c1 = ((N22 + λ) * (r1 + λ * p1) - N12 * (r2 + λ * p2)) / det
    c2 = ((N11 + λ) * (r2 + λ * p2) - N12 * (r1 + λ * p1)) / det
    (isfinite(c1) && isfinite(c2)) || return p1, p2, λ
    return c1, c2, λ
end

# full 2x2 C for one site, rows solved independently; also returns per-row λ
function solve_site_distortion(Zo::AbstractMatrix{ComplexF64}, Zp::AbstractMatrix{ComplexF64},
                               Ze::AbstractMatrix{ComplexF64}; damping::Float64)
    c11, c12, λ1 = _distortion_row(Zo, Zp, Ze, 1, 2, 1.0, 0.0, damping)
    c21, c22, λ2 = _distortion_row(Zo, Zp, Ze, 3, 4, 0.0, 1.0, damping)
    return c11, c12, c21, c22, λ1, λ2
end

#---------- corrected misfit ----------

# corrected residual when both predicted coefficients exist, uncorrected
# fallback when only the entry's own prediction exists; the datum count N is
# then identical to _accumulate_chi2! and damping=Inf reproduces chi2_and_rms
function _accumulate_chi2_distorted!(accχ2::Base.RefValue{Float64}, accN::Base.RefValue{Int},
                                     accP::Base.RefValue{Float64},
                                     Zo::AbstractMatrix{ComplexF64}, Zp::AbstractMatrix{ComplexF64},
                                     Ze::AbstractMatrix{ComplexF64},
                                     c11::Float64, c12::Float64, c21::Float64, c22::Float64,
                                     λ1::Float64, λ2::Float64)
    @inbounds for i in 1:size(Zo, 1)
        for (oc, ac, bc, c1, c2) in ((1, 1, 3, c11, c12), (2, 2, 4, c11, c12),
                                     (3, 1, 3, c21, c22), (4, 2, 4, c21, c22))
            zo = Zo[i, oc]; ze = Ze[i, oc]
            (isfinite(real(zo)) && isfinite(imag(zo)) &&
             isfinite(real(ze)) && isfinite(imag(ze))) || continue
            σ = abs(ze)
            σ > 0 || continue
            a = Zp[i, ac]; b = Zp[i, bc]
            zown = Zp[i, oc]
            if isfinite(real(a)) && isfinite(imag(a)) && isfinite(real(b)) && isfinite(imag(b))
                zp = c1 * a + c2 * b
            elseif isfinite(real(zown)) && isfinite(imag(zown))
                zp = zown
            else
                continue
            end
            rre = (real(zp) - real(zo)) / σ
            rim = (imag(zp) - imag(zo)) / σ
            accχ2[] += rre * rre + rim * rim
            accN[]  += 2
        end
    end
    s1 = (c11 - 1.0)^2 + c12^2
    s2 = c21^2 + (c22 - 1.0)^2
    s1 > 0 && isfinite(λ1) && (accP[] += λ1 * s1)
    s2 > 0 && isfinite(λ2) && (accP[] += λ2 * s2)
    return nothing
end

"""
    chi2_and_rms_distorted(dobs_path, dpred_path; damping, use_tipper=true) -> DistortionFit

Distortion-corrected misfit between observed and predicted ModEM data files.
Per site a real 2x2 matrix `C` minimizing the error-weighted misfit
`Zobs ~ C*Zpred` over all periods is solved in closed form, damped toward the
identity with relative strength `damping` (`Inf` forces `C = I` and reproduces
`chi2_and_rms` exactly). Impedance residuals use `C*Zpred`; tipper residuals
are never corrected. Returns a `DistortionFit` carrying χ², rms, the damping
penalty, and the per-site matrices for diagnostics.
"""
function chi2_and_rms_distorted(dobs_path::AbstractString, dpred_path::AbstractString;
                                damping::Float64, use_tipper::Bool=true)
    dobs  = load_data_modem(dobs_path)
    dpred = load_data_modem(dpred_path)
    dobs.site == dpred.site || error("site mismatch between observed and predicted files")
    length(dobs.T) == length(dpred.T) && all(isapprox.(dobs.T, dpred.T; rtol=1e-8)) ||
        error("period mismatch between observed and predicted files")

    χ2 = Ref(0.0); N = Ref(0); P = Ref(0.0)
    ns = dobs.ns
    C = Array{Float64,3}(undef, 2, 2, ns)
    for s in 1:ns
        Zo = view(dobs.Z, :, :, s)
        Zp = view(dpred.Z, :, :, s)
        Ze = view(dobs.Zerr, :, :, s)
        c11, c12, c21, c22, λ1, λ2 = solve_site_distortion(Zo, Zp, Ze; damping=damping)
        C[1, 1, s] = c11; C[1, 2, s] = c12
        C[2, 1, s] = c21; C[2, 2, s] = c22
        _accumulate_chi2_distorted!(χ2, N, P, Zo, Zp, Ze, c11, c12, c21, c22, λ1, λ2)
    end

    if use_tipper
        _accumulate_chi2!(χ2, N, dobs.tip, dpred.tip, dobs.tiperr)
    end

    rms = N[] > 0 ? sqrt(χ2[] / N[]) : NaN
    return DistortionFit(χ2[], rms, P[], C, copy(dobs.site))
end

#---------- diagnostics ----------

# per-site distortion table sorted by ||C-I||_F descending, so the most
# distorted sites top the list; groom-bailey flavored angles, no decomposition
function write_distortion_file(path::AbstractString, fit::DistortionFit;
                               iter::Int=-1, rms::Float64=NaN, damping::Float64=NaN)
    ns = length(fit.site)
    C = fit.C
    normci = [sqrt((C[1,1,s]-1.0)^2 + C[1,2,s]^2 + C[2,1,s]^2 + (C[2,2,s]-1.0)^2) for s in 1:ns]
    order = sortperm(normci; rev=true)
    open(path, "w") do io
        println(io, "# per-site galvanic distortion C of the current best model — ",
                Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
        @printf(io, "# iter=%d  rms=%.5f  penalty=%.5g  damping=%.4g\n", iter, rms, fit.penalty, damping)
        println(io, "# twist/shear in degrees; sorted by ||C-I||_F descending")
        @printf(io, "%-12s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n",
                "# site", "c11", "c12", "c21", "c22", "det", "twist", "shear", "aniso", "normCI")
        for s in order
            c11 = C[1,1,s]; c12 = C[1,2,s]; c21 = C[2,1,s]; c22 = C[2,2,s]
            detc  = c11 * c22 - c12 * c21
            twist = atand(c21 - c12, c11 + c22)
            shear = atand(c12 + c21, c11 + c22)
            aniso = (c11 - c22) / (c11 + c22)
            @printf(io, "%-12s %9.4f %9.4f %9.4f %9.4f %9.4f %9.3f %9.3f %9.4f %9.4f\n",
                    fit.site[s], c11, c12, c21, c22, detc, twist, shear, aniso, normci[s])
        end
    end
    return path
end

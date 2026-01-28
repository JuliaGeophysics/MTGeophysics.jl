struct _Chi2RMS
    chi2::Float64
    rms::Float64
end

# Internal: accumulate chi^2 over complex arrays (Re & Im counted separately)
function _accumulate_chi2!(accχ2::Base.RefValue{Float64}, accN::Base.RefValue{Int},
                           Zo::AbstractArray{ComplexF64}, Zp::AbstractArray{ComplexF64}, Ze::AbstractArray{ComplexF64})
    @assert size(Zo) == size(Zp) == size(Ze)
    @inbounds for I in eachindex(Zo)
        zo = Zo[I]; zp = Zp[I]; ze = Ze[I]
        if isfinite(real(zo)) && isfinite(imag(zo)) &&
           isfinite(real(zp)) && isfinite(imag(zp)) &&
           isfinite(real(ze)) && isfinite(imag(ze))
            σ = abs(ze)
            if σ > 0
                rre = (real(zp) - real(zo)) / σ
                rim = (imag(zp) - imag(zo)) / σ
                accχ2[] += rre*rre + rim*rim
                accN[]  += 2
            end
        end
    end
    return nothing
end

"""
    chi2_and_rms(dobs_path::AbstractString, dpred_path::AbstractString;
                 use_impedance::Bool=true, use_tipper::Bool=true,
                 components::Vector{String}=String[])

Load observed & predicted ModEM ASCII files and return **(chi2_total, rms_total)**
using σ from the observed file as weights. `components` (e.g. ["ZXY","ZYX"]) restricts
impedance terms; tipper is included/excluded via `use_tipper`.
"""
function chi2_and_rms(dobs_path::AbstractString, dpred_path::AbstractString;
                      use_impedance::Bool=true, use_tipper::Bool=true,
                      components::Vector{String}=String[])
    dobs  = load_data_modem(dobs_path)
    dpred = load_data_modem(dpred_path)

    χ2 = Ref(0.0); N = Ref(0)

    if use_impedance
        compmap = Dict("ZXX"=>1, "ZXY"=>2, "ZYX"=>3, "ZYY"=>4)
        sel = isempty(components) ? (1:4) : sort([compmap[c] for c in components if haskey(compmap,c)])
        if !isempty(sel)
            _accumulate_chi2!(χ2, N, view(dobs.Z, :, sel, :), view(dpred.Z, :, sel, :), view(dobs.Zerr, :, sel, :))
        end
    end

    if use_tipper
        _accumulate_chi2!(χ2, N, dobs.tip, dpred.tip, dobs.tiperr)
    end

    rms = N[] > 0 ? sqrt(χ2[] / N[]) : NaN
    return _Chi2RMS(χ2[], rms)
end

# Example:
# s = chi2_and_rms("dobs.dat", "dpred.dat";
#                  use_impedance=true, use_tipper=false, components=["ZXY","ZYX"])
# @show s.chi2 s.rms
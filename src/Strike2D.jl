# 2D strike and data rotation
# Author: @pankajkmishra
# Geoelectric strike from the phase tensor (Caldwell et al., 2004) and rotation of the impedances and station
# coordinates into the strike frame of the 2D codes (x along strike, y along the profile). Angles are degrees
# clockwise from geographic north, the rotation of the ModEM data header; fwd.ctrl `Strike (deg)` picks it

using Printf
using Statistics

# angle into [-p/2, p/2]
_wrap_angle(a, p) = a - p * round(a / p)

# major axis of the station positions, degrees clockwise from the frame's x; 90 (the frame's y) without one
function _mt2d_profile_azimuth(data::DataFile2D)
    x, y = data.x_positions, data.receivers
    length(y) >= 2 || return 90.0
    x, y = x .- mean(x), y .- mean(y)
    cxx, cyy, cxy = sum(abs2, x), sum(abs2, y), sum(x .* y)
    cxx + cyy > 0 || return 90.0
    rad2deg(atan(2cxy, cxx - cyy)) / 2
end

"""
    EstimateStrike2D(data::DataFile2D) -> NamedTuple or nothing

Geoelectric strike from the phase tensor strike α − β of every site and period with a
full impedance tensor, averaged over 4θ (strike is defined modulo 90°) with the phase
tensor ellipticity as weight, so 1D-like tensors count little. Of the two strike axes,
the one closest to perpendicular to the station line is taken, then turned by 180° if
needed so the rotation from the file frame stays within ±90°.

Returns `strike` (absolute, degrees clockwise from north), `rotation` (from the file
frame), `consistency` (0 = random, 1 = one strike everywhere), `skew` (median |β|,
degrees; small in 2D), `count` (site-periods used) and `profile` (station line azimuth).
Returns `nothing` without full tensors or when every tensor is isotropic (1D).
"""
function EstimateStrike2D(data::DataFile2D)
    s, w, β = 0.0im, 0.0, Float64[]
    for i in eachindex(data.z_xy)
        pt = phase_tensor(data.z_xx[i], data.z_xy[i], data.z_yx[i], data.z_yy[i])
        pt === nothing && continue
        e = abs(pt.ellipticity)
        isfinite(e) && isfinite(pt.azimuth) || continue
        s += e * cis(deg2rad(4pt.azimuth))
        w += e
        push!(β, abs(pt.beta))
    end
    w > 0 || return nothing
    θ = rad2deg(angle(s)) / 4
    φ = _mt2d_profile_azimuth(data)
    # strike perpendicular to the station line, then the turn within ±90° of the file frame
    θ = argmin(t -> abs(_wrap_angle(t - (φ - 90), 180)), (θ, θ + 90))
    θ = _wrap_angle(θ, 180)
    (strike = data.rotation + θ, rotation = θ, consistency = abs(s) / w, skew = median(β), count = length(β),
     profile = data.rotation + φ)
end

"""
    RotateData2D(data::DataFile2D, strike) -> DataFile2D

Data in the frame whose x points along `strike` (degrees clockwise from north): the
impedance tensor Z' = R Z Rᵀ, R = [cos θ sin θ; −sin θ cos θ], θ = `strike − data.rotation`,
the errors propagated as independent variances, and the local station x, y rotated
(y becomes the profile position; lat/lon and z are unchanged). Needs the full tensor:
periods of a site without ZXX and ZYY are dropped (NaN) with a warning. The turn is
recorded as a `:strike` step of the rotation history, as `rotate_data` does in 3D.
"""
function RotateData2D(data::DataFile2D, strike::Real)
    θ = _wrap_angle(Float64(strike) - data.rotation, 360)
    abs(θ) < 1e-9 && return data
    c, s = cosd(θ), sind(θ)
    out = deepcopy(data)
    dropped = 0
    for i in eachindex(data.z_xy)
        r = _rotate_impedance([data.z_xx[i] data.z_xy[i]; data.z_yx[i] data.z_yy[i]],
                              [data.z_xx_error[i] data.z_xy_error[i]; data.z_yx_error[i] data.z_yy_error[i]], θ)
        if r === nothing
            dropped += any(isfinite, (data.z_xy[i], data.z_yx[i]))
            for k in (:z_xx, :z_xy, :z_yx, :z_yy)
                getproperty(out, k)[i] = complex(NaN, NaN)
                getproperty(out, Symbol(k, :_error))[i] = NaN
            end
            continue
        end
        out.z_xx[i], out.z_xy[i], out.z_yx[i], out.z_yy[i] = r.Z[1, 1], r.Z[1, 2], r.Z[2, 1], r.Z[2, 2]
        out.z_xx_error[i], out.z_xy_error[i], out.z_yx_error[i], out.z_yy_error[i] = r.σ[1, 1], r.σ[1, 2], r.σ[2, 1], r.σ[2, 2]
    end
    any(isfinite, out.z_xy) || any(isfinite, out.z_yx) ||
        error("rotating by $(θ)° needs ZXX and ZYY, which $(isempty(data.path) ? "the data" : data.path) lacks")
    dropped > 0 && @warn "$dropped site-period(s) without ZXX and ZYY dropped by the rotation"
    for j in axes(out.z_xy, 2), k in axes(out.z_xy, 1)
        xy = _impedance_to_rho_phase(out.z_xy[k, j], out.frequencies[k])
        yx = _impedance_to_rho_phase(out.z_yx[k, j], out.frequencies[k])
        out.rho_xy[k, j], out.phase_xy[k, j], out.rho_yx[k, j], out.phase_yx[k, j] = xy.rho, xy.phase, yx.rho, yx.phase
    end
    x, y = data.x_positions, data.receivers
    out.x_positions .= c .* x .+ s .* y
    out.receivers .= -s .* x .+ c .* y
    fields = NamedTuple{fieldnames(DataFile2D)}(Tuple(getfield(out, k) for k in fieldnames(DataFile2D)))
    DataFile2D(; merge(fields, (; rotation = _wrap_angle(Float64(strike), 360),
                                 rotations = vcat(data.rotations, RotationStep(:strike, θ))))...)
end

"""
    StrikeData2D(data::DataFile2D, strike=nothing) -> (; data, strike, rotation, note)

The data in the 2D strike frame. `strike = nothing` (fwd.ctrl `Strike (deg) : auto`)
takes `EstimateStrike2D`; without full tensors or with 1D data the file frame is kept.
A number is the strike in degrees clockwise from north, taken as given. `rotation` is
the turn applied from the file frame, `note` a one-line account for logs and summaries.
"""
function StrikeData2D(data::DataFile2D, strike::Union{Nothing, Real} = nothing)
    if strike === nothing
        e = EstimateStrike2D(data)
        if e === nothing
            note = @sprintf("auto, no full anisotropic impedance tensors to estimate it; file frame kept (%.2f°)", data.rotation)
            return (; data, strike = data.rotation, rotation = 0.0, note)
        end
        note = @sprintf("%.2f° (auto: phase tensor, %d site-periods, consistency %.2f, median |β| %.1f°, station line %.1f°)",
                        e.strike, e.count, e.consistency, e.skew, e.profile)
        θ = round(e.strike; digits = 2)          # far finer than the estimate, and a clean header line
    else
        θ = Float64(strike)
        note = @sprintf("%.2f° (fwd.ctrl)", θ)
    end
    rotated = RotateData2D(data, θ)
    turn = _wrap_angle(θ - data.rotation, 360)
    abs(turn) < 1e-9 && (turn = 0.0)
    (; data = rotated, strike = θ, rotation = turn, note = note * @sprintf(", data rotated by %.2f°", turn))
end

"""
    RotateToStrike2D(data_path; strike=nothing, output_path=<stem>-r<ext>) -> (; path, strike, rotation, note)

Write the data rotated into the 2D strike frame (see `StrikeData2D`) as a full-tensor
ModEM file whose header rotation holds the strike; `data.obs` gives `data-r.obs`.
Nothing is written when no rotation is needed (`path` is then the input).
"""
function RotateToStrike2D(data_path::AbstractString; strike::Union{Nothing, Real} = nothing,
                          output_path::AbstractString = _rotated_path(data_path))
    r = StrikeData2D(load_data2d(data_path), strike)
    println("Strike: ", r.note)
    path = r.rotation == 0 ? String(data_path) : write_data2d(output_path, r.data; full_tensor = true)
    (; path, strike = r.strike, rotation = r.rotation, note = r.note)
end

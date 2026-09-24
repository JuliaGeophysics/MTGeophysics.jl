# Data rotation
# Author: @pankajkmishra
# Turn the impedance tensor and tipper of a ModEM data file by an angle, for magnetic declination, a rotated mesh or
# a 2D strike, and keep the history in the file; the same algebra serves 1D, 2D and 3D, the 2D strike code included

using Printf

const ROTATION_KINDS = (:declination, :mesh, :strike)

# rotation of the axes by θ degrees clockwise, x' = R x
_rotation_matrix(θ) = (c = cosd(θ); s = sind(θ); [c s; -s c])

# Z' = R Z Rᵀ with independent errors propagated as variances; nothing when any entry is missing
function _rotate_impedance(Z::AbstractMatrix, σ::AbstractMatrix, θ::Real)
    all(isfinite, Z) && all(isfinite, σ) || return nothing
    R = _rotation_matrix(θ)
    R², σ² = R .^ 2, σ .^ 2
    (Z = R * Z * R', σ = sqrt.(R² * σ² * R²'))
end

# Hz = T H, so T' = T Rᵀ for the row vector T = [Tzx Tzy]
function _rotate_tipper(T::AbstractVector, σ::AbstractVector, θ::Real)
    all(isfinite, T) && all(isfinite, σ) || return nothing
    R = _rotation_matrix(θ)
    (T = vec(transpose(T) * R'), σ = sqrt.(vec(transpose(σ .^ 2) * (R .^ 2)')))
end

# data-r.obs next to data.obs; a file already named <stem>-r keeps its name
function _rotated_path(path::AbstractString, dir::AbstractString = dirname(abspath(path)))
    stem, ext = splitext(basename(path))
    joinpath(dir, (endswith(stem, "-r") ? stem : stem * "-r") * ext)
end

_wrap_rotation(a) = (a = mod(a + 180, 360) - 180; a == -180 ? 180.0 : a)

"""
    rotate_data(d::Data, angle; kind=:mesh) -> Data
    rotate_data(path, angle; kind=:mesh, output_path=<stem>-r<ext>) -> path

Rotate a ModEM data file's impedance tensor (Z' = R Z Rᵀ) and tipper (T' = T Rᵀ),
R = [cos θ sin θ; −sin θ cos θ], with the errors propagated as independent variances.
The header angle is the azimuth of the data's x axis, degrees clockwise from geographic
north. Each call appends a step to the rotation history on the file's first `#` line.
- `kind = :mesh` or `:strike`: turn the axes by `angle` (one number). Station x, y
  turn with them, lat/lon stay, and the header angle grows by `angle`. Use it to match
  a rotated mesh or a 2D strike.
- `kind = :declination`: data recorded with x on magnetic north but labelled
  geographic. The fields turn by −`angle` (the declination, east positive), so x
  becomes geographic north. The station positions, already geographic, and the header
  angle stay. `angle` may be one number or one per site.

A site-period without the full tensor (or both tipper components) is dropped with a
warning; nothing rotatable is an error. The file method writes `data-r.obs` for
`data.obs` (a file already named `-r` is rewritten), in the file's own units and sign.
"""
function rotate_data(d::Data, angle; kind::Symbol = :mesh)
    kind in ROTATION_KINDS || throw(ArgumentError("kind must be one of $(ROTATION_KINDS)"))
    angles = angle isa Real ? fill(Float64(angle), d.ns) : Float64.(collect(angle))
    length(angles) == d.ns || throw(ArgumentError("$(length(angles)) angles for $(d.ns) sites"))
    angle isa Real || kind == :declination ||
        throw(ArgumentError("one angle per site is for declination only; a $kind rotation turns the whole frame"))
    kind == :declination && any(r -> r.kind == :declination, d.rotations) &&
        @warn "the data were already corrected for declination; correcting again"
    turn = kind == :declination ? -angles : angles
    out = deepcopy(d)
    dropped, rotated = 0, 0
    for is in 1:d.ns, ip in 1:d.nf
        θ = turn[is]
        z = d.Z[ip, :, is]
        r = _rotate_impedance(permutedims(reshape(z, 2, 2)), permutedims(reshape(abs.(d.Zerr[ip, :, is]), 2, 2)), θ)
        if r === nothing
            dropped += any(isfinite, z)
            out.Z[ip, :, is] .= complex(NaN, NaN)
            out.Zerr[ip, :, is] .= complex(NaN, NaN)
        else
            out.Z[ip, :, is] .= vec(permutedims(r.Z))
            out.Zerr[ip, :, is] .= vec(permutedims(r.σ))
            rotated += 1
        end
        t = d.tip[ip, :, is]
        r = _rotate_tipper(t, abs.(d.tiperr[ip, :, is]), θ)
        if r === nothing
            dropped += any(isfinite, t)
            out.tip[ip, :, is] .= complex(NaN, NaN)
            out.tiperr[ip, :, is] .= complex(NaN, NaN)
        else
            out.tip[ip, :, is] .= r.T
            out.tiperr[ip, :, is] .= r.σ
            rotated += 1
        end
    end
    rotated > 0 || error("nothing to rotate: no site-period has the full impedance tensor or both tipper components")
    dropped > 0 && @warn "$dropped site-period block(s) without the full tensor or both tipper components dropped by the rotation"
    if kind != :declination
        c, s = cosd(angles[1]), sind(angles[1])
        out.x, out.y = c .* d.x .+ s .* d.y, -s .* d.x .+ c .* d.y
        out.zrot = _wrap_rotation.(d.zrot .+ angles[1])
        out.trot = copy(out.zrot)
    end
    out.ρ, out.φ, out.ρerr, out.φerr = calc_rho_pha(out.Z, out.Zerr, out.T)
    step = angle isa Real ? RotationStep(kind, angle) : RotationStep(kind, angles, copy(d.site))
    out.rotations = vcat(d.rotations, step)
    out
end

function rotate_data(path::AbstractString, angle; kind::Symbol = :mesh,
                     output_path::AbstractString = _rotated_path(path))
    d = redirect_stdout(() -> load_data_modem(path; warn_rotation = false), devnull)
    out = rotate_data(d, angle; kind)
    conv = _read_data_convention(path)
    tipper = any(isfinite, out.tip)
    redirect_stdout(devnull) do
        write_data_modem(output_path, out; sign = conv.sign, units = conv.units, include_tipper = tipper,
                         description = "Rotated by MTGeophysics.jl")
    end
    @printf("Rotated %s -> %s: %s, header %.2f°\n", basename(path), basename(output_path),
            _rotation_text(out.rotations[end:end]), isempty(out.zrot) ? 0.0 : out.zrot[1])
    String(output_path)
end

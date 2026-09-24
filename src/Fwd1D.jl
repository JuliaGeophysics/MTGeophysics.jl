# 1D MT forward modelling
# Author: @pankajkmishra
# A 1D model is a one-column MT2DMesh (dimension = 1) in the ModEM-style model layout, solved exactly as a
# layered earth; the 2D data file, inversion driver, algorithms and plots are reused as they are
# Fréchet derivatives by forward-mode differentiation of the layer recursion

#***********************************************************************
# Description: 1D MT forward problem
#
#   convention   exp(+iωt), quasi-static, μ = μ₀, z positive down from the surface
#   layers       ρ₁ … ρₙ with thicknesses h₁ … hₙ₋₁, the last layer a halfspace
#   intrinsic    kⱼ = √(iωμ₀/ρⱼ),   ζⱼ = iωμ₀ / kⱼ = √(iωμ₀ρⱼ)
#   recursion    Zₙ = ζₙ
#                Zⱼ = ζⱼ (Zⱼ₊₁ + ζⱼ tanh(kⱼhⱼ)) / (ζⱼ + Zⱼ₊₁ tanh(kⱼhⱼ))
#   data         Zxy = Z₁,  Zyx = -Z₁,  Zxx = Zyy = 0
#                ρa = |Z|² / (ωμ₀),   φ = atan(Im Z, Re Z), 45° over a halfspace
#***********************************************************************

using ForwardDiff
using Printf
using Statistics

#---------- impedance mode ----------

# XY (Zxy), YX (Zyx = -Z), XYYX (both) or DET (the determinant impedance √det Z, fitted as Zxy)
_ctrl_mode1d(s) = (m = Symbol(uppercase(strip(String(s)))); m in (:XY, :YX, :XYYX, :DET) ? m :
                   throw(ArgumentError("1D mode must be XY, YX, XYYX or DET")))

# the 2D driver's mode for a 1D mode: XY fits z_xy (TE), YX z_yx (TM), DET its determinant as z_xy
_mt1d_driver_mode(m::Symbol) = m in (:XY, :DET) ? :TE : m == :YX ? :TM : :TETM

# curve labels of a 1D mode for the data plots: 1D has no TE or TM, only the impedance fitted
_mt1d_plot_names(m::Symbol) = m == :XY ? (TE = "XY", TM = nothing) : m == :YX ? (TE = nothing, TM = "YX") :
                              m == :DET ? (TE = "DET", TM = nothing) : (TE = "XY", TM = "YX")

#---------- mesh ----------

"""
    Mesh1D(z_cell_sizes, frequencies) -> MT2DMesh

One-column mesh of a layered earth: `z_cell_sizes` from the surface down, the last layer
a halfspace; one receiver on the surface. Solved exactly by the layer recursion.
"""
Mesh1D(z_cell_sizes::AbstractVector{<:Real}, frequencies::AbstractVector{<:Real}) =
    MT2DMesh(y_nodes = [-0.5, 0.5], z_nodes = vcat(0.0, cumsum(Float64.(z_cell_sizes))), y_cell_sizes = [1.0],
             z_cell_sizes = Float64.(z_cell_sizes), receiver_positions = [0.0], frequencies = Float64.(frequencies),
             n_air_cells = 0, dimension = 1)

"""
    Mesh1DFromInputs(model::ModelFile2D, data) -> (mesh, ρ)

1D mesh and resistivity `(nz, 1)` of a one-column ModEM-layout model for the
frequencies of `data`.
"""
function Mesh1DFromInputs(model::ModelFile2D, data)
    size(model.resistivity, 2) == 1 || throw(ArgumentError("a 1D model has one column, this one has $(size(model.resistivity, 2))"))
    any(>(MT2D_AIR_THRESHOLD), model.resistivity) && throw(ArgumentError("1D models hold no air"))
    Mesh1D(model.z_cell_sizes, data.frequencies), Matrix{Float64}(model.resistivity)
end

"""
    MakeMesh1D(data; mode=:XYYX, first_layer_div=5.0, vertical_factor=1.1, depth_mult=4.0) -> meshes

Layering of every site of `data` (a `DataFile2D` or a data file path) from its own skin
depths: the background is the site's median apparent resistivity over the periods (of
the impedance `mode` fits), the first layer δ(f_max)/`first_layer_div`, layers growing by
`vertical_factor` down to `depth_mult` δ(f_min). Returns one `(; site, thicknesses,
background)` per site; `Invert1D` starts each site from its background halfspace.
"""
function MakeMesh1D(data::DataFile2D; mode = :XYYX, first_layer_div::Real = 5.0, vertical_factor::Real = 1.1,
                    depth_mult::Real = 4.0)
    mode = _ctrl_mode1d(string(mode))
    map(eachindex(data.site_names)) do i
        site = mt1d_site_data(data, i; mode)
        ρa = mode == :YX ? site.rho_yx : mode == :XYYX ? vcat(site.rho_xy, site.rho_yx) : site.rho_xy
        v = filter(x -> isfinite(x) && x > 0, vec(ρa))
        background = isempty(v) ? 100.0 : 10^median(log10.(v))
        thicknesses = mt2d_geometric_layers(data.frequencies; background_resistivity = background, first_layer_div,
                                            vertical_factor, depth_mult)
        (; site = data.site_names[i], thicknesses, background)
    end
end

MakeMesh1D(data_path::AbstractString; kwargs...) = MakeMesh1D(load_data2d(data_path); kwargs...)

#---------- forward ----------

"""
    mt1d_impedance(frequencies, ρ, h) -> Vector

Surface impedance (ohm, exp(+iωt)) of layers `ρ` (ohm m, the last a halfspace) with
thicknesses `h` (one fewer), by the layer recursion; generic in the number type.
"""
function mt1d_impedance(frequencies::AbstractVector{<:Real}, ρ::AbstractVector, h::AbstractVector{<:Real})
    length(h) >= length(ρ) - 1 || throw(ArgumentError("need a thickness for every layer above the halfspace"))
    map(frequencies) do f
        iωμ = 1im * 2π * f * μ₀_2D
        Z = sqrt(iωμ * ρ[end])
        for j in length(ρ)-1:-1:1
            ζ = sqrt(iωμ * ρ[j])
            t = tanh(sqrt(iωμ / ρ[j]) * h[j])
            Z = ζ * (Z + ζ * t) / (ζ + Z * t)
        end
        Z
    end
end

# response and cache of a 1D mesh; J = ∂Z/∂ρ (nf × nz) when cache_fields
struct MT1DCache
    mesh::MT2DMesh
    response::MT2DResponse
    J::Matrix{ComplexF64}
    mode::Symbol
    air::BitMatrix
end

function _mt1d_forward_cache(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real}; mode = :TETM, cache_fields = true)
    mode in (:TE, :TM, :TETM) || throw(ArgumentError("mode must be :TE, :TM, or :TETM"))
    size(resistivity) == (length(mesh.z_cell_sizes), 1) || throw(DimensionMismatch("1D resistivity must be (nz, 1)"))
    all(x -> isfinite(x) && x > 0, resistivity) || throw(ArgumentError("earth resistivities must be finite and positive"))
    f, h = mesh.frequencies, mesh.z_cell_sizes
    ρ = Float64.(resistivity[:, 1])
    Z = mt1d_impedance(f, ρ, h)
    J = cache_fields ? _mt1d_jacobian(f, ρ, h) : zeros(ComplexF64, 0, 0)
    nf = length(f)
    col(v) = reshape(v, nf, 1)
    on(p) = mode in (p, :TETM)
    z_xy = on(:TE) ? col(Z) : zeros(ComplexF64, nf, 1)
    z_yx = on(:TM) ? col(-Z) : zeros(ComplexF64, nf, 1)
    ρa(z) = abs2.(z) ./ (2π .* f .* μ₀_2D)
    response = MT2DResponse(frequencies = f, periods = 1 ./ f, receivers = mesh.receiver_positions,
                            rho_xy = on(:TE) ? ρa(z_xy) : zeros(nf, 1), phase_xy = on(:TE) ? rad2deg.(angle.(z_xy)) : zeros(nf, 1),
                            z_xy = z_xy, rho_yx = on(:TM) ? ρa(z_yx) : zeros(nf, 1),
                            phase_yx = on(:TM) ? rad2deg.(angle.(z_yx)) : zeros(nf, 1), z_yx = z_yx)
    response, MT1DCache(mesh, response, J, mode, falses(length(h), 1))
end

function _mt1d_jacobian(f, ρ, h)
    D = ForwardDiff.jacobian(x -> (Z = mt1d_impedance(f, x, h); vcat(real.(Z), imag.(Z))), ρ)
    nf = length(f)
    complex.(D[1:nf, :], D[nf+1:end, :])
end

#---------- fréchet derivatives on the 2D interface ----------

# δd = G δρ, the six response fields as in 2D
function _mt2d_frechet(c::MT1DCache, δρ)
    δZ = c.J * δρ[:, 1]
    ω = 2π .* c.mesh.frequencies
    fields(z, δz) = (rho = reshape(2 .* real.(conj.(z) .* δz) ./ (ω .* μ₀_2D), :, 1),
                     phase = reshape(rad2deg.(imag.(δz ./ z)), :, 1), z = reshape(δz, :, 1))
    none = (rho = zeros(length(ω), 1), phase = zeros(length(ω), 1), z = zeros(ComplexF64, length(ω), 1))
    xy = c.mode in (:TE, :TETM) ? fields(c.response.z_xy[:, 1], δZ) : none
    yx = c.mode in (:TM, :TETM) ? fields(c.response.z_yx[:, 1], -δZ) : none
    (rho_xy = xy.rho, phase_xy = xy.phase, z_xy = xy.z, rho_yx = yx.rho, phase_yx = yx.phase, z_yx = yx.z)
end

# δρ̂ = Gᵗ δd̂ with the pairing of 2D: real(dot(δd̂, δd))
function _mt2d_frechet_transpose(c::MT1DCache, δd̂)
    ω = 2π .* c.mesh.frequencies
    ρ̂ = zeros(size(c.J, 2))
    dual(k) = hasproperty(δd̂, k) ? getproperty(δd̂, k)[:, 1] : zeros(length(ω))
    for (on, zk, rk, pk, s) in ((c.mode in (:TE, :TETM), :z_xy, :rho_xy, :phase_xy, 1),
                                (c.mode in (:TM, :TETM), :z_yx, :rho_yx, :phase_yx, -1))
        on || continue
        z = getproperty(c.response, zk)[:, 1]
        ẑ = dual(zk) .+ (2 ./ (ω .* μ₀_2D)) .* dual(rk) .* z .+ (180 / π) .* dual(pk) .* (1im ./ conj.(z))
        ρ̂ .+= s .* real.(c.J' * ẑ)
    end
    reshape(ρ̂, :, 1)
end

# rows of G: ∂ real(conj(w) z[key][i]) / ∂ρ
function _mt2d_frechet_rows(c::MT1DCache, rows)
    out = zeros(length(rows), size(c.J, 2))
    for (k, row) in enumerate(rows)
        f = CartesianIndices(size(c.response.z_xy))[row.index][1]
        s = row.key == :z_xy ? 1 : -1
        out[k, :] = s .* real.(conj(row.weight) .* c.J[f, :])
    end
    out
end

#---------- data ----------

"""
    mt1d_site_data(data::DataFile2D, i; mode=:XYYX) -> DataFile2D

The data of site `i` as a one-site 1D survey at y = 0. With `mode = :DET` the ZXY
column holds the determinant impedance √det Z (Zxx, Zyy taken as 0 where missing) with
the mean of the ZXY and ZYX errors.
"""
function mt1d_site_data(data::DataFile2D, i::Integer; mode::Symbol = :XYYX)
    pick(a) = a[:, i:i]
    z_xy, e_xy = copy(pick(data.z_xy)), copy(pick(data.z_xy_error))
    if mode == :DET
        zero_nan(z) = isfinite(z) ? z : zero(z)
        det = zero_nan.(pick(data.z_xx)) .* zero_nan.(pick(data.z_yy)) .- pick(data.z_xy) .* pick(data.z_yx)
        z_xy = sqrt.(det)
        e_xy = (pick(data.z_xy_error) .+ pick(data.z_yx_error)) ./ 2
    end
    rp = [_impedance_to_rho_phase(z_xy[k], data.frequencies[k]) for k in axes(z_xy, 1), _ in 1:1]
    DataFile2D(title = data.title, periods = data.periods, frequencies = data.frequencies,
               site_names = data.site_names[i:i], receivers = [0.0], x_positions = data.x_positions[i:i],
               z_positions = data.z_positions[i:i], z_xy = z_xy, z_xy_error = e_xy,
               z_yx = copy(pick(data.z_yx)), z_yx_error = copy(pick(data.z_yx_error)),
               z_xx = copy(pick(data.z_xx)), z_xx_error = copy(pick(data.z_xx_error)),
               z_yy = copy(pick(data.z_yy)), z_yy_error = copy(pick(data.z_yy_error)),
               rho_xy = getfield.(rp, :rho), phase_xy = getfield.(rp, :phase),
               rho_yx = copy(pick(data.rho_yx)), phase_yx = copy(pick(data.phase_yx)), path = data.path,
               latitudes = isempty(data.latitudes) ? Float64[] : data.latitudes[i:i],
               longitudes = isempty(data.longitudes) ? Float64[] : data.longitudes[i:i], origin = data.origin)
end

"""
    ForwardSolve1D(model_path, data_path; mode=:XYYX, write_frechet=false, output_path,
                   add_noise=false, rng_seed=20260308) -> path

1D forward run: a one-column ModEM-layout model and a data file giving the sites, periods
and errors. Every site gets the response of the model, ZXY = Z and ZYX = -Z; `mode`
(`:XY`, `:YX`, `:XYYX`, `:DET`) picks what is written. Writes `data.pred` next to the data
by default, with `write_frechet` also G as `<output>.frechet`. A data file whose
impedances are all zero is a template: its errors are fractions of |Z|.
"""
function ForwardSolve1D(model_path::AbstractString, data_path::AbstractString; mode = :XYYX, write_frechet::Bool = false,
                        output_path::AbstractString = joinpath(dirname(abspath(data_path)), "data.pred"),
                        add_noise::Bool = false, rng_seed::Integer = 20260308)
    mode = _ctrl_mode1d(string(mode))
    template = load_data2d(data_path)
    mesh, ρ = Mesh1DFromInputs(ReadModel2D(model_path), template)
    one, _ = _mt1d_forward_cache(mesh, ρ; mode = :TETM, cache_fields = false)
    ns = length(template.receivers)
    rep(a) = repeat(a, 1, ns)
    response = MT2DResponse(frequencies = one.frequencies, periods = one.periods, receivers = template.receivers,
                            rho_xy = rep(one.rho_xy), phase_xy = rep(one.phase_xy), z_xy = rep(one.z_xy),
                            rho_yx = rep(one.rho_yx), phase_yx = rep(one.phase_yx), z_yx = rep(one.z_yx))
    errors = _resolve_forwardsolve2d_errors(template, response)
    predicted = data_from_response2d(response; errors...,
        site_names = template.site_names, x_positions = template.x_positions, z_positions = template.z_positions,
        latitudes = template.latitudes, longitudes = template.longitudes, origin = template.origin)
    mode == :XY && (predicted.z_yx .= NaN)
    mode == :YX && (predicted.z_xy .= NaN)
    written = write_data2d(output_path, add_noise ? _apply_mt2d_noise(predicted; rng_seed) : predicted)
    if write_frechet
        site = mt1d_site_data(predicted, 1; mode)
        WriteFrechet2D(splitext(written)[1] * ".frechet", mesh, ρ, site; mode = _mt1d_driver_mode(mode))
    end
    written
end

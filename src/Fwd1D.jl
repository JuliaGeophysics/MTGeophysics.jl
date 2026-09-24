# 1D MT forward modelling
# Author: @pankajkmishra
# A 1D model is a list of layer thicknesses and resistivities, stored as a one-column ModEM-style model file,
# solved exactly by the layer recursion. Only the data and model file formats are shared with 2D; the
# layering, forward, Fréchet derivatives and inversion (Inv1D.jl) are 1D's own, so 2D changes cannot break it
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
using Random
using Statistics

const μ₀_1D = 4π * 1e-7

#---------- impedance mode ----------

# XY (Zxy), YX (Zyx = -Z), XYYX (both) or DET (the determinant impedance √det Z, fitted as Zxy)
_ctrl_mode1d(s) = (m = Symbol(uppercase(strip(String(s)))); m in (:XY, :YX, :XYYX, :DET) ? m :
                   throw(ArgumentError("1D mode must be XY, YX, XYYX or DET")))

# data columns a mode fits: (key, error key, sign of Z)
_mt1d_components(m::Symbol) = m in (:XY, :DET) ? ((:z_xy, :z_xy_error, 1),) : m == :YX ? ((:z_yx, :z_yx_error, -1),) :
                              ((:z_xy, :z_xy_error, 1), (:z_yx, :z_yx_error, -1))

# curve labels of a 1D mode for the data plots: 1D has no TE or TM, only the impedance fitted
_mt1d_plot_names(m::Symbol) = m == :XY ? (TE = "XY", TM = nothing) : m == :YX ? (TE = nothing, TM = "YX") :
                              m == :DET ? (TE = "DET", TM = nothing) : (TE = "XY", TM = "YX")

#---------- layering ----------

mt1d_skin_depth(ρ::Real, f::Real) = sqrt(2 * ρ / (2π * f * μ₀_1D))

"""
    mt1d_layers(frequencies; background_resistivity=100.0, first_layer_div=5.0,
                vertical_factor=1.1, depth_mult=4.0) -> thicknesses

Layer thicknesses from the skin depths in `background_resistivity`: the first layer
δ(f_max)/`first_layer_div`, each next one `vertical_factor` thicker, down to
`depth_mult` δ(f_min). The last layer is the halfspace.
"""
function mt1d_layers(frequencies::AbstractVector{<:Real}; background_resistivity::Real = 100.0,
                     first_layer_div::Real = 5.0, vertical_factor::Real = 1.1, depth_mult::Real = 4.0)
    first_layer_div > 0 && vertical_factor >= 1 && depth_mult > 0 ||
        throw(ArgumentError("need first_layer_div > 0, vertical_factor >= 1, depth_mult > 0"))
    h = [mt1d_skin_depth(background_resistivity, maximum(frequencies)) / first_layer_div]
    bottom = depth_mult * mt1d_skin_depth(background_resistivity, minimum(frequencies))
    while sum(h) < bottom
        push!(h, h[end] * vertical_factor)
    end
    h
end

"""
    MakeMesh1D(data; mode=:XYYX, first_layer_div=5.0, vertical_factor=1.1, depth_mult=4.0) -> meshes

Layering of every site of `data` (a `DataFile2D` or a data file path) from its own skin
depths: the background is the site's median apparent resistivity over the periods (of
the impedance `mode` fits), then `mt1d_layers`. Returns one `(; site, thicknesses,
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
        thicknesses = mt1d_layers(data.frequencies; background_resistivity = background, first_layer_div,
                                  vertical_factor, depth_mult)
        (; site = data.site_names[i], thicknesses, background)
    end
end

MakeMesh1D(data_path::AbstractString; kwargs...) = MakeMesh1D(load_data2d(data_path); kwargs...)

#---------- forward ----------

"""
    mt1d_impedance(frequencies, ρ, h) -> Vector

Surface impedance (ohm, exp(+iωt)) of layers `ρ` (ohm m, the last a halfspace) with
thicknesses `h` (one fewer, extra ones ignored), by the layer recursion; generic in the
number type.
"""
function mt1d_impedance(frequencies::AbstractVector{<:Real}, ρ::AbstractVector, h::AbstractVector{<:Real})
    length(h) >= length(ρ) - 1 || throw(ArgumentError("need a thickness for every layer above the halfspace"))
    map(frequencies) do f
        iωμ = 1im * 2π * f * μ₀_1D
        Z = sqrt(iωμ * ρ[end])
        for j in length(ρ)-1:-1:1
            ζ = sqrt(iωμ * ρ[j])
            t = tanh(sqrt(iωμ / ρ[j]) * h[j])
            Z = ζ * (Z + ζ * t) / (ζ + Z * t)
        end
        Z
    end
end

"""
    mt1d_frechet(frequencies, ρ, h) -> G

Fréchet derivative G = ∂Z/∂log10 ρ (`nf × nlayers`, complex) of the surface impedance.
"""
function mt1d_frechet(frequencies::AbstractVector{<:Real}, ρ::AbstractVector{<:Real}, h::AbstractVector{<:Real})
    nf = length(frequencies)
    D = ForwardDiff.jacobian(m -> (Z = mt1d_impedance(frequencies, 10 .^ m, h); vcat(real.(Z), imag.(Z))), log10.(ρ))
    complex.(D[1:nf, :], D[nf+1:end, :])
end

# apparent resistivity and phase of an impedance
_mt1d_rho_phase(z, f) = isfinite(z) ? (abs2(z) / (2π * f * μ₀_1D), rad2deg(angle(z))) : (NaN, NaN)

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
    rp = [_mt1d_rho_phase(z_xy[k], data.frequencies[k]) for k in axes(z_xy, 1), _ in 1:1]
    DataFile2D(title = data.title, periods = data.periods, frequencies = data.frequencies,
               site_names = data.site_names[i:i], receivers = [0.0], x_positions = data.x_positions[i:i],
               z_positions = data.z_positions[i:i], z_xy = z_xy, z_xy_error = e_xy,
               z_yx = copy(pick(data.z_yx)), z_yx_error = copy(pick(data.z_yx_error)),
               z_xx = copy(pick(data.z_xx)), z_xx_error = copy(pick(data.z_xx_error)),
               z_yy = copy(pick(data.z_yy)), z_yy_error = copy(pick(data.z_yy_error)),
               rho_xy = first.(rp), phase_xy = last.(rp),
               rho_yx = copy(pick(data.rho_yx)), phase_yx = copy(pick(data.phase_yx)), path = data.path,
               latitudes = isempty(data.latitudes) ? Float64[] : data.latitudes[i:i],
               longitudes = isempty(data.longitudes) ? Float64[] : data.longitudes[i:i], origin = data.origin, rotation = data.rotation,
               rotations = data.rotations)
end

# the survey of `data` with Z at every site from one impedance column per site (nf × ns) and the errors of
# `data`; a template's (all impedances zero) are fractions of |Z|
function _mt1d_predicted(data::DataFile2D, Z::AbstractMatrix;
                         fractional::Bool = all(iszero, data.z_xy) && all(iszero, data.z_yx))
    err(e, z) = fractional ? max.(e .* abs.(z), 1e-6) : copy(e)
    z_xy, z_yx = ComplexF64.(Z), ComplexF64.(-Z)
    rp(z) = [_mt1d_rho_phase(z[k, j], data.frequencies[k]) for k in axes(z, 1), j in axes(z, 2)]
    xy, yx = rp(z_xy), rp(z_yx)
    DataFile2D(title = data.title, periods = data.periods, frequencies = data.frequencies, site_names = data.site_names,
               receivers = data.receivers, x_positions = data.x_positions, z_positions = data.z_positions,
               z_xy = z_xy, z_xy_error = err(data.z_xy_error, z_xy), z_yx = z_yx, z_yx_error = err(data.z_yx_error, z_yx),
               z_xx = fill(complex(NaN), size(Z)), z_xx_error = fill(NaN, size(Z)),
               z_yy = fill(complex(NaN), size(Z)), z_yy_error = fill(NaN, size(Z)),
               rho_xy = first.(xy), phase_xy = last.(xy), rho_yx = first.(yx), phase_yx = last.(yx),
               latitudes = data.latitudes, longitudes = data.longitudes, origin = data.origin, rotation = data.rotation,
               rotations = data.rotations)
end

# complex gaussian noise of each datum's error
function _mt1d_add_noise(data::DataFile2D, rng_seed::Integer)
    rng = MersenneTwister(rng_seed)
    out = deepcopy(data)
    for k in eachindex(out.z_xy)
        out.z_xy[k] += out.z_xy_error[k] / sqrt(2) * (randn(rng) + im * randn(rng))
        out.z_yx[k] += out.z_yx_error[k] / sqrt(2) * (randn(rng) + im * randn(rng))
        out.rho_xy[k], out.phase_xy[k] = _mt1d_rho_phase(out.z_xy[k], out.frequencies[CartesianIndices(out.z_xy)[k][1]])
        out.rho_yx[k], out.phase_yx[k] = _mt1d_rho_phase(out.z_yx[k], out.frequencies[CartesianIndices(out.z_yx)[k][1]])
    end
    out
end

# layer thicknesses and resistivities of a one-column model file
function _mt1d_read_model(path::AbstractString)
    model = ReadModel2D(path)
    size(model.resistivity, 2) == 1 || throw(ArgumentError("a 1D model has one column, this one has $(size(model.resistivity, 2))"))
    any(>(1e15), model.resistivity) && throw(ArgumentError("1D models hold no air"))
    Float64.(model.z_cell_sizes), Float64.(model.resistivity[:, 1])
end

"""
    WriteFrechet1D(path, h, ρ, site::DataFile2D; mode=:XYYX) -> path

G = ∂d/∂log10 ρ of a one-site survey's impedances: two rows (Re, Im) per period and fitted
component, one column per layer, in [mV/km]/[nT].
"""
function WriteFrechet1D(path::AbstractString, h::AbstractVector{<:Real}, ρ::AbstractVector{<:Real}, site::DataFile2D;
                        mode::Symbol = :XYYX)
    G = mt1d_frechet(site.frequencies, ρ, h) ./ (μ₀_1D * 1000)
    comps = [(k, s, k == :z_xy ? (mode == :DET ? "DET" : "ZXY") : "ZYX") for (k, _, s) in _mt1d_components(mode)]
    rows = [(k, s, c, i) for i in eachindex(site.frequencies) for (k, s, c) in comps if isfinite(getproperty(site, k)[i, 1])]
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# MTGeophysics 1D Fréchet derivative G = ∂d/∂m")
        println(io, "# d: impedance in [mV/km]/[nT], two rows (Re, Im) per data line, exp(+iωt)")
        println(io, "# m: log10 resistivity of the layers, from the surface down")
        println(io, "# > rows layers, then: period site component part g_1 ... g_layers")
        @printf(io, "> %d %d\n", 2 * length(rows), length(ρ))
        for (k, s, c, i) in rows, (part, v) in (("Re", real.(s .* G[i, :])), ("Im", imag.(s .* G[i, :])))
            print(io, @sprintf("%.8e %s %s %s", site.periods[i], site.site_names[1], c, part))
            foreach(x -> print(io, @sprintf(" %.6e", x)), v)
            println(io)
        end
    end
    String(path)
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
    h, ρ = _mt1d_read_model(model_path)
    Z = mt1d_impedance(template.frequencies, ρ, h)
    predicted = _mt1d_predicted(template, repeat(Z, 1, length(template.site_names)))
    add_noise && (predicted = _mt1d_add_noise(predicted, rng_seed))
    mode == :XY && (predicted.z_yx .= NaN)
    mode == :YX && (predicted.z_xy .= NaN)
    written = write_data2d(output_path, predicted)
    write_frechet && WriteFrechet1D(splitext(written)[1] * ".frechet", h, ρ, mt1d_site_data(predicted, 1; mode); mode)
    written
end

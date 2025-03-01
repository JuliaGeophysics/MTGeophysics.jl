#=
This function computes the magnetotelluric (MT) response given resistivity and thickness vectors.
The implementation is based on the equations from the book:
Ward, S. H., & Hohmann, G. W. (1988). Electromagnetic Theory for Geophysical Applications. In M. N. Nabighian (Ed.), 
Electromagnetic Methods in Applied Geophysics, Volume 1: Theory. Society of Exploration Geophysicists.

Author: pankajkmishra 
Last Edit: 2025-02-21
=#

function MT1D_response(frequencies, mesh, resistivities)
    μ0 = 4π * 1e-7  # Free-space permeability (H/m) (Eq. 3.96)
    N = length(mesh)
    
    if length(resistivities) != N - 1
        error("Length of resistivities must be equal to length(mesh) - 1")
    end

    thicknesses = diff(mesh)  # Compute thickness from mesh
    σ = 1.0 ./ resistivities  # Conductivity (Eq. 3.2)
    ω = 2π .* frequencies  # Angular frequency (Eq. 3.96)

    k = sqrt.(-1im .* ω .* μ0 .* σ[end])  # Wavenumber for bottom half-space (Eq. 3.96)
    Z = (ω .* μ0) ./ k  # Intrinsic impedance for last layer (Eq. 3.96)

    for i in (N-2):-1:1
        k_i = sqrt.(-1im .* ω .* μ0 .* σ[i])  # Wavenumber in layer (Eq. 3.96)
        Z_i = (ω .* μ0) ./ k_i  # Intrinsic impedance of layer (Eq. 3.96)
        tanh_term = tanh.(1im .* k_i .* thicknesses[i])  # Tangent hyperbolic term (Eq. 3.121)
        Z = Z_i .* (Z .+ Z_i .* tanh_term) ./ (Z_i .+ Z .* tanh_term)  # Recursive impedance (Eq. 3.121)
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)  # Apparent resistivity (Eq. 3.149)
    φ = atan.(imag.(Z), real.(Z)) .* (180 / π)  # Impedance phase (Eq. 3.150)
    
    return ρa, φ
end 


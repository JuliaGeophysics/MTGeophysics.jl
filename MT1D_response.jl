#=
This function computes the magnetotelluric (MT) response given resistivity and thickness vectors.
The implementation is based on the equations from the book:
Ward, S. H., & Hohmann, G. W. (1988). Electromagnetic Theory for Geophysical Applications. In M. N. Nabighian (Ed.), 
Electromagnetic Methods in Applied Geophysics, Volume 1: Theory. Society of Exploration Geophysicists.

Author: pankajkmishra 
Last Edit: 2025-02-21


# Arguments
- `frequencies`: Vector of frequencies (Hz)
- `resistivities`: Vector of resistivities (Ω·m)
- `thicknesses`: Vector of thicknesses (m)

# Returns
- `ρa`: Apparent resistivity (Ω·m)
- `φ`: Phase (degrees)
=#

function MT1D_response(frequencies, resistivities, thicknesses)
    μ0 = 4π * 1e-7  # Free-space permeability (H/m) (Ward and Hohmann, 1988, Eq. 3.96)
    m = length(resistivities)          
    σ = 1.0 ./ resistivities           
    ω = 2π .* frequencies              

    k = sqrt.(-1im .* ω .* μ0 .* σ[end])  
    Z = (ω .* μ0) ./ k  # impedance (Eq. 3.96)

    for i in (m-1):-1:1
        k_i = sqrt.(-1im .* ω .* μ0 .* σ[i])  # Wave number (Eq. 3.96)
        Z_i = (ω .* μ0) ./ k_i  # Intrinsic impedance (Eq. 3.96)
        tanh_term = tanh.(1im .* k_i .* thicknesses[i])  
        Z = Z_i .* (Z .+ Z_i .* tanh_term) ./ (Z_i .+ Z .* tanh_term)  # Recursive impedance (Eq. 3.121)
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)  # Apparent resistivity (Eq. 3.149)
    φ = atan.(imag.(Z), real.(Z)) .* (180 / π)  
    return ρa, φ
end 


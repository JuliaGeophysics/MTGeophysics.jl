#=
This function computes the magnetotelluric (MT) response given resistivity and thickness vectors.
The implementation is based on Finite-difference method 

Author: Pankaj K Mishra (pankaj.mishra@gtk.fi)
Date: 2025-02-14


# Arguments
- `frequencies`: Vector of frequencies (Hz)
- `resistivities`: Vector of resistivities (Ω·m)
- `thicknesses`: Vector of thicknesses (m)
- `dz`: Grid spacing (m)
- `z_max`: Maximum depth (m) 

# Returns
- `ρa`: Apparent resistivity (Ω·m)
- `φ`: Phase (degrees) 

=# 


function MT1D_response_FD(frequencies, resistivities, thicknesses; dz=5.0, z_max=6000.0)
    μ0 = 4π * 1e-7
    ω = 2π .* frequencies
    Z = Vector{ComplexF64}(undef, length(frequencies))
    
    for (j, f_val) in enumerate(frequencies)
        ω_val = 2π * f_val
        z = 0:dz:z_max
        N = length(z)
        
       
        sigma = zeros(ComplexF64, N)
        z_boundary1 = thicknesses[1]
        z_boundary2 = thicknesses[1] + thicknesses[2]
        for (i, z_val) in enumerate(z)
            if z_val < z_boundary1
                sigma[i] = 1 / resistivities[1]
            elseif z_val < z_boundary2
                sigma[i] = 1 / resistivities[2]
            else
                sigma[i] = 1 / resistivities[3]
            end
        end
        
       
        k2 = -1im * ω_val * μ0 .* sigma
        
        # FD system for the unknowns at grid points 2..N
        n_unknowns = N - 1
        A = zeros(ComplexF64, n_unknowns, n_unknowns)
        b = zeros(ComplexF64, n_unknowns)
        
        # Equation at grid point 2 (using E[1]=1 as the Dirichlet BC)
        A[1, 1] = -2/dz^2 + k2[2]
        if n_unknowns > 1
            A[1, 2] = 1/dz^2
        end
        b[1] = -1/dz^2 * 1.0  
        
        # Interior grid points 3..N-1.
        for i in 2:(n_unknowns-1)
            A[i, i-1] = 1/dz^2
            A[i, i]   = -2/dz^2 + k2[i+1]
            A[i, i+1] = 1/dz^2
        end
        
        # Bottom BC using a one-sided derivative with a radiation BC
        sigma_bottom = 1 / resistivities[end]
        k_bottom = sqrt(-1im * ω_val * μ0 * sigma_bottom)
        A[n_unknowns, n_unknowns-1] = -1/dz
        A[n_unknowns, n_unknowns]   = 1/dz + 1im*k_bottom
        b[n_unknowns] = 0.0
        
        # Solve the linear system for the electric field (excluding the known surface value)
        U = A \ b
        
        # Construct the E field vector with E(1)=1.
        E = ComplexF64[1.0]
        append!(E, U)
        
        # FD derivative at the surface
        dE_dz0 = (E[2] - E[1]) / dz
        
        # Z from the surface BC.
        Z[j] = -1im * ω_val * μ0 * E[1] / dE_dz0
    end
    
   
    ρa = abs.(Z).^2 ./ (ω .* μ0)
    φ = atan.(imag.(Z), real.(Z)) .* (180 / π)
    return ρa, φ
end

#=
Computes the magnetotelluric (MT) response using a finite-difference method with a non-uniform grid.

This function calculates the MT response based on provided resistivity and thickness vectors. Unlike the earlier MT1D_response_FD.jl implementation that used a uniform grid, this version employs a non-uniform grid—where the spacing increases with depth—and uses LinearSolve.jl to solve the resulting linear system.

Author: Pankaj K Mishra
Date: 2025-02-21

# What goes in 
- frequencies: Vector of frequencies (Hz)
- resistivities: Vector of resistivities (Ω·m)
- thicknesses: Vector of layer thicknesses (m)
- dz: (Optional) Initial grid spacing (m)
- z_max: (Optional) Maximum depth (m)
- z_scale: (Optional) Scaling factor for increasing grid spacing default is 1.05 

# What comes out 
- ρa: Apparent resistivity (Ω·m)
- φ: Phase (degrees)

=# 

using LinearSolve

function MT1D_response_FD_II(frequencies, resistivities, thicknesses; dz=nothing, z_max=nothing, z_scale=1.05)
    μ0 = 4π * 1e-7
   
    if dz === nothing
        f_max = maximum(frequencies)
        skin_depth_top = sqrt(2 * resistivities[1] / (2π * f_max * μ0))
        resolution_factor = 20.0
        dz = skin_depth_top / resolution_factor
    end
    if z_max === nothing
        z_max = 100000.0
    end

  
    z = [0.0]
    current_dz = dz
    while z[end] < z_max
        push!(z, z[end] + current_dz)
        current_dz *= z_scale
    end

    # Precompute conductivity based on the depth grid
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

    ω = 2π .* frequencies
    Z = Vector{ComplexF64}(undef, length(frequencies))

    # Frequency loop (Does it help to parallelize this loop?, maybe not for 1D)
    for (j, f_val) in enumerate(frequencies)
        ω_val = 2π * f_val
        k2 = -1im * ω_val * μ0 .* sigma
        n_unknowns = N - 1
        A = zeros(ComplexF64, n_unknowns, n_unknowns)
        b = zeros(ComplexF64, n_unknowns)
        if N >= 3
            h1 = z[2] - z[1]
            h2 = z[3] - z[2]
            A[1,1] = -2/(h1*h2) + k2[2]
            if n_unknowns > 1
                A[1,2] = 2/(h2*(h1+h2))
            end
            b[1] = -2/(h1*(h1+h2))
            for i in 3:(N-1)
                j_idx = i - 1
                h_previous = z[i] - z[i-1]
                h_next = z[i+1] - z[i]
                A[j_idx, j_idx-1] = 2/(h_previous*(h_previous+h_next))
                A[j_idx, j_idx]   = -2/(h_previous*h_next) + k2[i]
                A[j_idx, j_idx+1] = 2/(h_next*(h_previous+h_next))
            end
        else
            h = z[2] - z[1]
            A[1,1] = 1/h + 1im*sqrt(-1im*ω_val*μ0*(1/resistivities[end]))
            b[1] = -1/h
        end
        h_last = z[end] - z[end-1]
        sigma_bottom = 1 / resistivities[end]
        k_bottom = sqrt(-1im * ω_val * μ0 * sigma_bottom)
        A[end, end-1] = -1/h_last
        A[end, end]   = 1/h_last + 1im*k_bottom
        b[end] = 0.0
        prob = LinearProblem(A, b)
        # LinearSolve.jl can do many cool things here. diferent solvers!
        sol = solve(prob)
        U = sol.u
        E = ComplexF64[1.0]
        append!(E, U)
        h0 = z[2] - z[1]
        dE_dz0 = (E[2] - E[1]) / h0
        Z[j] = -1im * ω_val * μ0 * E[1] / dE_dz0
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)
    φ = atan.(imag.(Z), real.(Z)) .* (180/π)
    return ρa, φ
end

#=
This is a FD solver for the 1D magnetotelluric (MT) response. This reads a UBC 1D
mesh and model file and computes the apparent resistivity and phase response. 
- @pankajkmishra  
=#
using LinearSolve
function MT1D_response_FD(frequencies, mesh, model)
    μ0 = 4π * 1e-7  
    N = length(mesh)
    
  
    if mesh[1] != 0.0
        error("The depth mesh must start at 0.")
    end
    if length(model) != N - 1
        error("Length of model (resistivities) must equal length(mesh) - 1")
    end

    
    interfaces = mesh[2:end-1]  
    
    # Assign conductivity (σ = 1/ρ) at each mesh node
    sigma = zeros(ComplexF64, N)
    for (i, z_val) in enumerate(mesh)
        # For depths shallower than the first interface, assign cell  
        if isempty(interfaces) || z_val < interfaces[1]
            cell_index = 1
        else
            cell_index = searchsortedfirst(interfaces, z_val)
            if cell_index > length(model)
                cell_index = length(model)
            end
        end
        sigma[i] = 1 / model[cell_index]
    end

    ω = 2π .* frequencies
    Z = Vector{ComplexF64}(undef, length(frequencies))

    # Loop over frequencies
    for (j, f_val) in enumerate(frequencies)
        ω_val = 2π * f_val
        k2 = -1im * ω_val * μ0 .* sigma
        n_unknowns = N - 1
        A = zeros(ComplexF64, n_unknowns, n_unknowns)
        b = zeros(ComplexF64, n_unknowns)
        
        # Assemble the FD matrix 
        if N >= 3
            h1 = mesh[2] - mesh[1]
            h2 = mesh[3] - mesh[2]
            A[1,1] = -2/(h1*h2) + k2[2]
            if n_unknowns > 1
                A[1,2] = 2/(h2*(h1+h2))
            end
            b[1] = -2/(h1*(h1+h2))
            for i in 3:(N-1)
                j_idx = i - 1
                h_prev = mesh[i] - mesh[i-1]
                h_next = mesh[i+1] - mesh[i]
                A[j_idx, j_idx-1] = 2/(h_prev*(h_prev+h_next))
                A[j_idx, j_idx]   = -2/(h_prev*h_next) + k2[i]
                A[j_idx, j_idx+1] = 2/(h_next*(h_prev+h_next))
            end
        else
            h = mesh[2] - mesh[1]
            A[1,1] = 1/h + 1im*sqrt(-1im*ω_val*μ0*sigma[end])
            b[1] = -1/h
        end
        
        # Bottom BC 
        h_last = mesh[end] - mesh[end-1]
        sigma_bottom = sigma[end]
        k_bottom = sqrt(-1im * ω_val * μ0 * sigma_bottom)
        A[end, end-1] = -1/h_last
        A[end, end]   = 1/h_last + 1im*k_bottom
        b[end] = 0.0
        
        prob = LinearProblem(A, b)
        sol = solve(prob)
        U = sol.u
        
        # Impose the surface boundary condition: E(0) = 1
        E = ComplexF64[1.0]
        append!(E, U)
        h0 = mesh[2] - mesh[1]
        dE_dz0 = (E[2] - E[1]) / h0
        Z[j] = -1im * ω_val * μ0 * E[1] / dE_dz0
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)
    φ = atan.(imag.(Z), real.(Z)) .* (180/π)
    return ρa, φ
end

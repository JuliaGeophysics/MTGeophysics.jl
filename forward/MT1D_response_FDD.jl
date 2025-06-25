# FD computations in pure Julia so Enzyme can differentiate them 
# Once this is tested well, we can delete the previous FD one 
#  --- @Pankajkmishra 
function MT1D_response_FDD(frequencies::Vector{Float64}, mesh::Vector{Float64}, model::Vector{Float64})
    μ0 = 4π * 1e-7  
    N = length(mesh)

    if mesh[1] != 0.0
        error("The depth mesh must start at 0.")
    end
    if length(model) != N - 1
        error("Length of model (resistivities) must equal length(mesh) - 1")
    end

    interfaces = mesh[2:end-1]  
    sigma = zeros(ComplexF64, N)
    for (i, z_val) in enumerate(mesh)
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

    for (j, f_val) in enumerate(frequencies)
        ω_val = 2π * f_val
        k2 = -1im * ω_val * μ0 .* sigma
        n_unknowns = N - 1
        A = zeros(ComplexF64, n_unknowns, n_unknowns)
        b = zeros(ComplexF64, n_unknowns)

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

        h_last = mesh[end] - mesh[end-1]
        sigma_bottom = sigma[end]
        k_bottom = sqrt(-1im * ω_val * μ0 * sigma_bottom)
        A[end, end-1] = -1/h_last
        A[end, end]   = 1/h_last + 1im*k_bottom
        b[end] = 0.0

     
        U = A \ b

        E = ComplexF64[1.0]
        append!(E, U)
        h0 = mesh[2] - mesh[1]
        dE_dz0 = (E[2] - E[1]) / h0
        Z[j] = -1im * ω_val * μ0 * E[1] / dE_dz0
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)
    φ  = atan.(imag.(Z), real.(Z)) .* (180/π)
    return ρa, φ
end
#=
using Enzyme
function objective(model::Vector{Float64})
    frequencies = [1.0, 10.0, 100.0]
    mesh = [0.0, 50.0, 100.0, 200.0]
    ρa, φ = MT1D_response_FDD(frequencies, mesh, model)
    return sum(ρa)  # A dummy objective function 
end

model = [100.0, 200.0, 300.0]

# Compute gradient using Enzyme's reverse-mode AD
grad_model = Enzyme.gradient(set_runtime_activity(Reverse), objective, model)
println("Gradient of objective with respect to model: ", grad_model)
=#
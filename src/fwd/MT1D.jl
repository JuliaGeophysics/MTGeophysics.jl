# Author: Pankaj K Mishra (@pankajmishra), Last Updated: December 24, 2025 

abstract type MT1DMethod end
struct Analytical <: MT1DMethod end
struct FiniteDifference <: MT1DMethod end

function MT1D(frequencies::Vector{<:Real}, ρ::Vector{<:Real}, mesh::Vector{<:Real}, method::MT1DMethod=Analytical())
    freq = convert(Vector{Float64}, frequencies)
    resistivities = convert(Vector{Float64}, ρ)
    depths = convert(Vector{Float64}, mesh)
    return _mt1d_forward(freq, resistivities, depths, method)
end

function _mt1d_forward(frequencies::Vector{Float64}, ρ::Vector{Float64}, mesh::Vector{Float64}, ::Analytical)
    μ0 = 4π * 1e-7
    N = length(mesh)
    length(ρ) != N - 1 && error("Length of ρ must equal length(mesh) - 1")
    
    thicknesses = diff(mesh)
    σ = 1.0 ./ ρ
    ω = 2π .* frequencies
    k = sqrt.(-1im .* ω .* μ0 .* σ[end])
    Z = (ω .* μ0) ./ k

    for i in (N-2):-1:1
        k_i = sqrt.(-1im .* ω .* μ0 .* σ[i])
        Z_i = (ω .* μ0) ./ k_i
        tanh_term = tanh.(1im .* k_i .* thicknesses[i])
        Z = Z_i .* (Z .+ Z_i .* tanh_term) ./ (Z_i .+ Z .* tanh_term)
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)
    φ = atan.(imag.(Z), real.(Z)) .* (180/π)
    return ρa, φ
end

function _mt1d_forward(frequencies::Vector{Float64}, ρ::Vector{Float64}, mesh::Vector{Float64}, ::FiniteDifference)
    μ0 = 4π * 1e-7
    mesh[1] != 0.0 && error("The depth mesh must start at 0 for finite-difference method")
    length(ρ) != length(mesh) - 1 && error("Length of ρ must equal length(mesh) - 1")

    # Generate fine mesh for FD (refine between layer boundaries)
    fine_mesh, fine_ρ = _refine_mesh(mesh, ρ)
    N = length(fine_mesh)

    # Assign conductivity at each mesh node
    interfaces = fine_mesh[2:end-1]
    sigma = zeros(ComplexF64, N)
    for (i, z_val) in enumerate(fine_mesh)
        if isempty(interfaces) || z_val < interfaces[1]
            cell_index = 1
        else
            cell_index = searchsortedfirst(interfaces, z_val)
            cell_index > length(fine_ρ) && (cell_index = length(fine_ρ))
        end
        sigma[i] = 1 / fine_ρ[cell_index]
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
            h1, h2 = fine_mesh[2] - fine_mesh[1], fine_mesh[3] - fine_mesh[2]
            A[1,1] = -2/(h1*h2) + k2[2]
            n_unknowns > 1 && (A[1,2] = 2/(h2*(h1+h2)))
            b[1] = -2/(h1*(h1+h2))
            for i in 3:(N-1)
                j_idx = i - 1
                h_prev, h_next = fine_mesh[i] - fine_mesh[i-1], fine_mesh[i+1] - fine_mesh[i]
                A[j_idx, j_idx-1] = 2/(h_prev*(h_prev+h_next))
                A[j_idx, j_idx]   = -2/(h_prev*h_next) + k2[i]
                A[j_idx, j_idx+1] = 2/(h_next*(h_prev+h_next))
            end
        else
            h = fine_mesh[2] - fine_mesh[1]
            A[1,1] = 1/h + 1im*sqrt(-1im*ω_val*μ0*sigma[end])
            b[1] = -1/h
        end

        h_last = fine_mesh[end] - fine_mesh[end-1]
        k_bottom = sqrt(-1im * ω_val * μ0 * sigma[end])
        A[end, end-1] = -1/h_last
        A[end, end] = 1/h_last + 1im*k_bottom
        b[end] = 0.0

        U = A \ b
        E = ComplexF64[1.0]
        append!(E, U)
        h0 = fine_mesh[2] - fine_mesh[1]
        dE_dz0 = (E[2] - E[1]) / h0
        Z[j] = -1im * ω_val * μ0 * E[1] / dE_dz0
    end

    ρa = abs.(Z).^2 ./ (ω .* μ0)
    φ = atan.(imag.(Z), real.(Z)) .* (180/π)
    return ρa, φ
end

function _refine_mesh(mesh::Vector{Float64}, ρ::Vector{Float64}; d0=1.0, r=1.05, npad=20)
    fine_mesh = Float64[0.0]
    fine_ρ = Float64[]
    
    # Refine each layer with geometric spacing
    for i in 1:length(ρ)
        z_start, z_end = mesh[i], mesh[i+1]
        thickness = z_end - z_start
        
        # Generate geometrically spaced points within layer
        dz = d0
        z = z_start + dz
        while z < z_end
            push!(fine_mesh, z)
            push!(fine_ρ, ρ[i])
            dz *= r
            z += dz
        end
        push!(fine_mesh, z_end)
        push!(fine_ρ, ρ[i])
    end
    
    # Add padding below (extend half-space)
    dz = fine_mesh[end] - fine_mesh[end-1]
    for _ in 1:npad
        dz *= r
        push!(fine_mesh, fine_mesh[end] + dz)
        push!(fine_ρ, ρ[end])
    end
    
    return fine_mesh, fine_ρ
end

function MT1D(frequencies::Vector{<:Real}, ρ::Vector{<:Real}, mesh::Vector{<:Real}, method::Union{Symbol,String})
    m = lowercase(string(method))
    m == "analytical" && return MT1D(frequencies, ρ, mesh, Analytical())
    m == "fd" && return MT1D(frequencies, ρ, mesh, FiniteDifference())
    error("Unknown method: '$method'. Use :analytical or :fd")
end

#=
Creates UBC 1D format model and mesh filesfor a given model. adds padding as well 
- @pankajkmishra  
=#


using DelimitedFiles
using Plots

cd(@__DIR__)

model_boundaries = [0.0, 1000.0, 1500.0, 3000.0, 5000.0, 8000.0]
resistivities = [1000.0, 100.0, 1200.0, 200.0, 900.0]

function generate_fd_mesh(model_boundaries::Vector{Float64}; d0::Float64=1.0, r::Float64=1.05, pad::Int=15)
    fd_mesh = Float64[]
    current_z = model_boundaries[1]
    push!(fd_mesh, current_z)
    d = d0
    for i in 2:length(model_boundaries)
        target = model_boundaries[i]
        while current_z + d < target
            current_z += d
            push!(fd_mesh, current_z)
            d *= r
        end
        if abs(current_z - target) > 1e-8
            current_z = target
            push!(fd_mesh, current_z)
            d *= r
        end
    end
    for j in 1:pad
        current_z += d
        push!(fd_mesh, current_z)
        d *= r
    end
    return fd_mesh
end

fd_mesh = generate_fd_mesh(model_boundaries; d0=1.0, r=1.05, pad=20)

function get_resistivity(z::Float64, boundaries::Vector{Float64}, resistivities::Vector{Float64})
    for i in 1:length(resistivities)
        if z < boundaries[i+1]
            return resistivities[i]
        end
    end
    return resistivities[end]
end

fd_resistivities = Float64[]
for i in 1:length(fd_mesh)-1
    z_mid = (fd_mesh[i] + fd_mesh[i+1]) / 2.0
    push!(fd_resistivities, get_resistivity(z_mid, model_boundaries, resistivities))
end

open("model.mesh", "w") do io
    println(io, length(fd_mesh))
    for z in fd_mesh
        println(io, z)
    end
end

open("model.rho", "w") do io
    println(io, length(fd_resistivities))
    for r in fd_resistivities
        println(io, r)
    end
end

println("Model saved to model.mesh and model.rho.")



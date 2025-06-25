using DelimitedFiles

# Function to create UBC 1D format mesh and model files
function UBC_1D(model_boundaries::Vector{Float64}, resistivities::Vector{Float64};
                          d0::Float64=1.0, r::Float64=1.05, pad::Int=20)
    # Generate finite-difference mesh with padding
    function generate_fd_mesh(boundaries::Vector{Float64}; d0::Float64, r::Float64, pad::Int)
        fd_mesh = Float64[]
        current_z = boundaries[1]
        push!(fd_mesh, current_z)
        d = d0
        for i in 2:length(boundaries)
            target = boundaries[i]
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
        # Add padding below the last model boundary
        for j in 1:pad
            current_z += d
            push!(fd_mesh, current_z)
            d *= r
        end
        return fd_mesh
    end

    # Determine the resistivity for a given depth based on the model boundaries
    function get_resistivity(z::Float64, boundaries::Vector{Float64}, res::Vector{Float64})
        for i in 1:length(res)
            if z < boundaries[i+1]
                return res[i]
            end
        end
        return res[end]
    end

    # Generate the finite-difference mesh
    fd_mesh = generate_fd_mesh(model_boundaries; d0=d0, r=r, pad=pad)
    
    # Compute the resistivity for each cell (using the midpoint of each cell)
    fd_resistivities = Float64[]
    for i in 1:length(fd_mesh)-1
        z_mid = (fd_mesh[i] + fd_mesh[i+1]) / 2.0
        push!(fd_resistivities, get_resistivity(z_mid, model_boundaries, resistivities))
    end

    # Write the mesh file in UBC format
    open("model.mesh", "w") do io
        println(io, length(fd_mesh))
        for z in fd_mesh
            println(io, z)
        end
    end

    # Write the resistivity file in UBC format
    open("model.rho", "w") do io
        println(io, length(fd_resistivities))
        for r_val in fd_resistivities
            println(io, r_val)
        end
    end

    println("Model saved to model.mesh and model.rho.")
end

# Example usage:
#model_boundaries = [0.0, 1000.0, 1500.0, 3000.0, 5000.0, 8000.0]
#resistivities    = [1000.0, 100.0, 1200.0, 200.0, 900.0]
#create_ubc_model(model_boundaries, resistivities; d0=1.0, r=1.05, pad=20)

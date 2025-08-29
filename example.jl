# MTGeophysics.jl Usage Example

# This demonstrates basic usage of the MTGeophysics package

using Pkg
Pkg.activate(".")
using MTGeophysics

println("=== MTGeophysics.jl Example ===")

# Example 1: 1D MT forward modeling
println("\n1. MT 1D Forward Modeling:")

# Define a simple 3-layer earth model
frequencies = 10.0 .^ range(-2, 2, length=50)  # 0.01 to 100 Hz
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]  # depths in meters  
resistivities = [100.0, 10.0, 1000.0]  # resistivity in Ω⋅m

# Compute MT response using analytical solution
rho_a, phi = MT1D_response(frequencies, layer_boundaries, resistivities)

println("   Computed $(length(frequencies)) frequency points")
println("   Frequency range: $(minimum(frequencies)) - $(maximum(frequencies)) Hz")
println("   Apparent resistivity range: $(round(minimum(rho_a), digits=2)) - $(round(maximum(rho_a), digits=2)) Ω⋅m")
println("   Phase range: $(round(minimum(phi), digits=2))° - $(round(maximum(phi), digits=2))°")

# Example 2: Create UBC format files
println("\n2. UBC Format File Creation:")

# Create temporary directory for output
import Pkg: mkpath
tmpdir = mktempdir()
cd(tmpdir)

# Generate UBC mesh and model
UBC_1D(layer_boundaries, resistivities; d0=10.0, r=1.1, pad=15)

if isfile("model.mesh") && isfile("model.rho")
    println("   ✓ Created UBC format files: model.mesh and model.rho")
    mesh_lines = length(readlines("model.mesh"))
    model_lines = length(readlines("model.rho"))
    println("   ✓ Mesh has $mesh_lines nodes")
    println("   ✓ Model has $model_lines cells")
else
    println("   ✗ Error creating UBC files")
end

# Example 3: Data structures
println("\n3. Data Structures:")

# Create empty ModEM data structure
data = make_nan_data() 
println("   ✓ Created ModEMData structure")
println("   ✓ Initial sites: $(data.ns)")
println("   ✓ Initial frequencies: $(data.nf)")

println("\n=== Package successfully loaded and working! ===")
#=
Compares the analytical vs FD solver for 1D MT
@Pankajkmishra 
=#
using DelimitedFiles
using Plots
include("UBC_1D.jl") 
cd(@__DIR__)
include("MT1D_response.jl")
#include("MT1D_response_FD.jl")
include("MT1D_response_FDD.jl")


println("Reading the user input....")
model_boundaries = [0.0, 1000.0, 1500.0, 3000.0, 5000.0, 8000.0]
resistivities = [1000.0, 100.0, 1200.0, 200.0, 900.0] 


println("Making a UBC mesh and model from user inputs....")
UBC_1D(model_boundaries, resistivities; d0=1.0, r=1.05, pad=20)


println("Reading the UBC mesh and model we just created....")
function read_mesh(filename::String)
    lines = readlines(filename)
    [parse(Float64, line) for line in lines[2:end]]
end

function read_model(filename::String)
    lines = readlines(filename)
    [parse(Float64, line) for line in lines[2:end]]
end

fd_mesh_read = read_mesh("model.mesh")
fd_resistivities_read = read_model("model.rho")


println("Setting up frequency vector....")
frequencies = 10 .^ range(-3, 3, length=100)


println("Computing Analytical forward response....")
ρa_ana, φ_ana = MT1D_response(frequencies, model_boundaries, resistivities)
#ρa_FD, φ_FD = MT1D_response_FD(frequencies, fd_mesh_read, fd_resistivities_read)

println("Computing Finite-difference forward response ....")
ρa_FD, φ_FD = MT1D_response_FDD(frequencies, fd_mesh_read, fd_resistivities_read)






println("Plotting and saving ....")
p1 = plot(frequencies, ρa_ana, xscale=:log10, yscale=:log10, xflip=true, xlabel="Frequency (Hz)", ylabel="Apparent Resistivity (Ω·m)", lw=2, label="Analytical")
plot!(p1, frequencies, ρa_FD, lw=2, label="FD")
xlims!(p1, 1e-3, 1e3)
ylims!(p1, 1, 1e5)

p2 = plot(frequencies, φ_ana, xscale=:log10, xflip=true, xlabel="Frequency (Hz)", ylabel="Phase (°)", lw=2, legend=false)
plot!(p2, frequencies, φ_FD, lw=2)
xlims!(p2, 1e-3, 1e3)
ylims!(p2, 0, 90)

x_vals = Float64[]
y_vals = Float64[]
for i in 1:length(resistivities)
    push!(x_vals, resistivities[i])
    push!(y_vals, model_boundaries[i])
    push!(x_vals, resistivities[i])
    push!(y_vals, model_boundaries[i+1])
end
p3 = plot(x_vals, y_vals, seriestype=:steppre, yflip=true, xscale=:log10, xlabel="Resistivity (Ω·m)", ylabel="Depth (m)", lw=2, legend=false)
xlims!(p3, 1e0, 1e5)
ylims!(p3, 0, fd_mesh_read[end])
yticks!(p3, 0:500:fd_mesh_read[end], string.(0:500:fd_mesh_read[end]))

idx_pad = findfirst(x -> isapprox(x, model_boundaries[end], atol=1e-8), fd_mesh_read)
pad_x = Float64[]
pad_y = Float64[]
for i in idx_pad:length(fd_mesh_read)-1
    push!(pad_x, resistivities[end])
    push!(pad_y, fd_mesh_read[i])
    push!(pad_x, resistivities[end])
    push!(pad_y, fd_mesh_read[i+1])
end
plot!(p3, pad_x, pad_y, seriestype=:steppre, lw=2, linestyle=:dash, color=:black, label="Padding")

P = plot(p1, p2, p3, layout=@layout([grid(2,1) grid(1,1)]), size=(800,800))
display(P)
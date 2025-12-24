#= Compares the analytical vs FD solver for 1D MT - @Pankajkmishra =#
using Plots
cd(@__DIR__)
include("../src/fwd/MT1D.jl") 

# Define layer boundaries and resistivities
mesh = [0.0, 1000.0, 1500.0, 3000.0, 5000.0, 8000.0]
ρ = [1000.0, 100.0, 1200.0, 200.0, 900.0]
frequencies = 10 .^ range(-3, 3, length=40)

ρa_ana, φ_ana = MT1D(frequencies, ρ, mesh, :analytical)
ρa_FD, φ_FD = MT1D(frequencies, ρ, mesh, :fd)

# Plot
p1 = plot(frequencies, ρa_ana, xscale=:log10, yscale=:log10, xflip=true,  xlabel="Frequency (Hz)", ylabel="Apparent Resistivity (Ω·m)", lw=3, label="Analytical")
plot!(p1, frequencies, ρa_FD, lw=2, ls=:dash, label="FD")
xlims!(p1, 1e-3, 1e3); ylims!(p1, 1, 1e5)

p2 = plot(frequencies, φ_ana, xscale=:log10, xflip=true, xlabel="Frequency (Hz)", ylabel="Phase (°)", lw=3, label="Analytical")
plot!(p2, frequencies, φ_FD, lw=2, ls=:dash, label="FD")
xlims!(p2, 1e-3, 1e3); ylims!(p2, 0, 90)

x_vals, y_vals = Float64[], Float64[]
for i in 1:length(ρ)
    push!(x_vals, ρ[i]); push!(y_vals, mesh[i])
    push!(x_vals, ρ[i]); push!(y_vals, mesh[i+1])
end
p3 = plot(x_vals, y_vals, seriestype=:steppre, yflip=true, xscale=:log10, xlabel="Resistivity (Ω·m)", ylabel="Depth (m)", lw=2, legend=false)
xlims!(p3, 1e0, 1e5); ylims!(p3, 0, mesh[end])

P = plot(p1, p2, p3, layout=@layout([grid(2,1) grid(1,1)]), size=(800,800))
display(P)
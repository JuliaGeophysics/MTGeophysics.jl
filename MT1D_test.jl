using Plots
include("MT1D_forward_analytical.jl")


ρ = [100.0, 10.0, 100.0]  # Resistivities (Ω·m)
h = [2000.0, 1000.0]      # Thicknesses (m)
f = 10 .^ range(-3, 2, length=100)  


ρa, φ = MT1D_response(f, ρ, h)


p1 = plot(f, ρa, xscale=:log10, yscale=:log10, linewidth=2, xflip=true,
    xlabel="Frequency (Hz)", ylabel="Apparent Resistivity",
    title="Apparent Resistivity", legend=false)
xlims!(p1, 10^-3, 10^2)  
ylims!(p1, 1, 1000)  


p2 = plot(f, φ, xscale=:log10, linewidth=2, xflip=true,
    xlabel="Frequency (Hz)", ylabel="Phase (degrees)",
    title="Phase", legend=false)
xlims!(p2, 10^-3, 10^2)  
ylims!(p2, 0, 90)  


plot(p1, p2, layout=(2, 1), size=(800, 600))

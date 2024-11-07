using Plots
using LinearAlgebra

function MT1D(T, rho, h)
    mu0 = 4 * π * 1e-7
    m = length(rho)
    k = zeros(Complex{Float64}, m, length(T))
    
    for n in 1:m
        k[n, :] = sqrt.(-im * 2 * π * mu0 ./ (T .* rho[n]))
    end
    
    Z = -(im * mu0 * 2 * π) ./ (T .* k[m, :])
    
    for n in (m-1):-1:1
        A = -(im * mu0 * 2 * π) ./ (T .* k[n, :])
        B = exp.(-2 .* k[n, :] .* h[n])
        Z = A .* (A .* (1 .- B) .+ Z .* (1 .+ B)) ./ (A .* (1 .+ B) .+ Z .* (1 .- B))
    end
    
    rho_a = (T ./ (mu0 * 2 * π)) .* abs.(Z).^2
    phase = -atan.(imag.(Z) ./ real.(Z)) .* 180 / π
    
    return rho_a, phase
end

# Test parameters
T = 10 .^ range(-3, 3, length=100)
rho = [100.0, 10.0, 500.0, 50.0]
h = [500.0, 300.0, 1000.0]

rho_a, phase = MT1D(T, rho, h)

# Set up the plot with subplots
p1 = plot(T, rho_a, xscale=:log10, yscale=:log10, xlabel="Period (s)", ylabel="Apparent Resistivity (Ohm-m)", label="Apparent Resistivity", linewidth=2, title="Apparent Resistivity")
p2 = plot(T, phase, xscale=:log10, xlabel="Period (s)", ylabel="Phase (degrees)", label="Phase", linewidth=2, title="Phase")

# Combine and save the plots
plot(p1, p2, layout=(2, 1), size=(800, 600))
savefig("MT1D_results.png")
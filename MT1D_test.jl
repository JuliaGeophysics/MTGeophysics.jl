using Plots
include("MT1D_response.jl")  ;
include("MT1D_response_FD.jl") ;

ρ = [100.0, 10.0, 100.0] ;
h = [2000.0, 1000.0] ;     
f = 10 .^ range(-3, 2, length=100);


ρa_ana, φ_ana = MT1D_response(f, ρ, h);
ρa_fd, φ_fd   = MT1D_response_FD(f, ρ, h; dz=5.0, z_max=6000.0);


p1 = plot(f, ρa_ana, xscale=:log10, yscale=:log10, linewidth=2, xflip=true,
    xlabel="Frequency (Hz)", ylabel="Apparent Resistivity (Ω·m)",
    title="Apparent Resistivity", label="Analytical");
plot!(p1, f, ρa_fd, linestyle=:dash, linewidth=2, label="FD");
xlims!(p1, 10^-3, 10^2);
ylims!(p1, 1, 1000);


p2 = plot(f, φ_ana, xscale=:log10, linewidth=2, xflip=true,
    xlabel="Frequency (Hz)", ylabel="Phase (degrees)",
    title="Phase", label="Analytical");
plot!(p2, f, φ_fd, linestyle=:dash, linewidth=2, label="FD");
xlims!(p2, 10^-3, 10^2);
ylims!(p2, 0, 90);

plot(p1, p2, layout=(2, 1), size=(800,600)) 
# 1D forward solver
# Author: @pankajkmishra
# Ensures the layered impedance matches the half-space and closed-form two-layer answers and the 2D solver on a
# laterally uniform model, and that the Fréchet derivative matches finite differences

using Test, LinearAlgebra, Random

@testset "1D forward" begin
    f = collect(10 .^ range(-2, 2, length = 12))

    @testset "layered forward" begin
        Z = mt1d_impedance(f, [100.0], Float64[])
        @test abs2.(Z) ./ (2π .* f .* 4π * 1e-7) ≈ fill(100.0, 12)
        @test all(rad2deg.(angle.(Z)) .≈ 45)
        # two layers: the top at high frequency, the basement at low frequency
        Z2 = mt1d_impedance([1e4, 1e-4], [10.0, 1000.0], [500.0])
        ρa = abs2.(Z2) ./ (2π .* [1e4, 1e-4] .* 4π * 1e-7)
        @test ρa[1] ≈ 10 rtol = 1e-3
        # closed-form two-layer impedance; 50 S of cover still pulls the lowest frequency below 1000
        two(ω) = (k = ρ -> sqrt(1im * ω * 4π * 1e-7 / ρ); Zi = ρ -> 1im * ω * 4π * 1e-7 / k(ρ);
                  t = tanh(k(10.0) * 500.0); Zi(10.0) * (Zi(1000.0) + Zi(10.0) * t) / (Zi(10.0) + Zi(1000.0) * t))
        @test ρa ≈ abs2.(two.(2π .* [1e4, 1e-4])) ./ (2π .* [1e4, 1e-4] .* 4π * 1e-7) rtol = 1e-8
        @test 900 < ρa[2] < 1000

        # the 2D finite-difference solver on a laterally uniform model agrees
        mesh2 = BuildMesh2D(frequencies = f, y_core_range = (-2000.0, 2000.0), y_core_cell = 500.0, y_padding = 20_000.0,
                            air_top = -50_000.0, air_cells = 10, receiver_positions = [0.0],
                            ground_layers = mt2d_geometric_layers(f; first_layer_div = 10.0, vertical_factor = 1.05))
        layers = [(0.0, 150.0, 100.0), (150.0, 500.0, 20.0), (500.0, Inf, 500.0)]
        ρ2 = build_mt2d_halfspace_model(mesh2)
        zc = MTGeophysics.mt2d_z_centers(mesh2)
        for iz in mesh2.n_air_cells+1:length(zc), (top, bottom, ρ) in layers
            top <= zc[iz] < bottom && (ρ2[iz, :] .= ρ)
        end
        r2 = run_mt2d_forward(mesh2, ρ2)
        Z1 = mt1d_impedance(f, ρ2[mesh2.n_air_cells+1:end, 1], mesh2.z_cell_sizes[mesh2.n_air_cells+1:end])
        @test maximum(abs.(abs2.(Z1) ./ (2π .* f .* 4π * 1e-7) ./ r2.rho_xy[:, 1] .- 1)) < 0.03
        @test maximum(abs.(rad2deg.(angle.(Z1)) .- r2.phase_xy[:, 1])) < 1.0
    end

    @testset "Fréchet derivatives" begin
        rng = MersenneTwister(3)
        h = fill(100.0, 12)
        ρ = 100 .* exp.(0.5randn(rng, 12))
        δ = randn(rng, 12)
        G = mt1d_frechet(f, ρ, h)
        ε = 1e-6
        fd = (mt1d_impedance(f, ρ .* 10 .^ (ε .* δ), h) .- mt1d_impedance(f, ρ .* 10 .^ (-ε .* δ), h)) ./ 2ε
        @test G * δ ≈ fd rtol = 1e-6
        @test size(G) == (12, 12)
        @test sum(mt1d_layers(f; background_resistivity = 100.0)) >= 4 * mt1d_skin_depth(100.0, minimum(f))
    end
end

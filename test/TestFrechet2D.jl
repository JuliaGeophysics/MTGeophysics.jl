# 2D Fréchet derivatives
# Author: @pankajkmishra
# Ensures G, Gδm and Gᵗδd̂ match central finite differences and pass the dot-product test for TE, TM and TETM in
# every parameterization, flat and with topography, that air cells have zero columns and bad inputs are refused

using Test, LinearAlgebra, Random

@testset "2D Fréchet derivatives" begin
    mesh = BuildMesh2D(frequencies=[0.3, 3.0], y_core_range=(-500.,500.),
        y_core_cell=250., y_padding=500., air_cells=2, air_top=-1000.,
        ground_layers=[100.,300.,900.], receiver_positions=[-375.,125.,400.])
    rho = build_mt2d_halfspace_model(mesh)
    rng = MersenneTwister(127)
    rho[3:end,:] .*= exp.(0.3randn(rng, size(rho,1)-2, size(rho,2)))
    direction = randn(rng, size(rho))
    fields = (:rho_xy, :phase_xy, :z_xy, :rho_yx, :phase_yx, :z_yx)
    h = 1e-4

    @testset "Directional derivatives, all response fields and modes" begin
        for mode in (:TE, :TM, :TETM)
            tangent = ApplyFrechet2D(mesh, rho, direction; mode, parameterization=:log_resistivity)
            rp = run_mt2d_forward(mesh, rho .* exp.(h*direction); mode)
            rm = run_mt2d_forward(mesh, rho .* exp.(-h*direction); mode)
            for key in fields
                fd = (getproperty(rp,key) - getproperty(rm,key))/(2h)
                @test getproperty(tangent,key) ≈ fd rtol=2e-5 atol=1e-10
            end
            weights = (rho_xy=randn(rng,2,3), phase_xy=randn(rng,2,3),
                z_xy=randn(rng,ComplexF64,2,3), rho_yx=randn(rng,2,3),
                phase_yx=randn(rng,2,3), z_yx=randn(rng,ComplexF64,2,3))
            gradient = ApplyFrechetTranspose2D(mesh, rho, weights; mode, parameterization=:log_resistivity)
            @test dot(gradient,direction) ≈ sum(real(dot(getproperty(weights,k),getproperty(tangent,k))) for k in fields) rtol=1e-8
            @test all(iszero,gradient[1:mesh.n_air_cells,:])
        end
    end

    @testset "Fréchet ordering, boundaries, parameterizations and fixed air" begin
        cells = [CartesianIndex(3,1), CartesianIndex(5,size(rho,2)),
                 CartesianIndex(3,4), CartesianIndex(4,5), CartesianIndex(1,2)]
        for parameterization in (:resistivity,:log_resistivity,:log10_resistivity)
            sens = FrechetDerivative2D(mesh,rho;active_cells=cells,parameterization)
            @test sens.cells == cells
            @test size(sens.z_xy) == (6,5)
            for (j,cell) in enumerate(cells)
                plus, minus = copy(rho), copy(rho)
                if parameterization == :resistivity
                    plus[cell] += h; minus[cell] -= h
                elseif parameterization == :log_resistivity
                    plus[cell] *= exp(h); minus[cell] *= exp(-h)
                else
                    plus[cell] *= 10.0^h; minus[cell] *= 10.0^-h
                end
                rp,rm = run_mt2d_forward(mesh,plus),run_mt2d_forward(mesh,minus)
                for key in (:z_xy,:z_yx)
                    fd = vec(getproperty(rp,key)-getproperty(rm,key))/(2h)
                    @test getproperty(sens,key)[:,j] ≈ fd rtol=2e-4 atol=1e-10
                end
            end
            @test all(iszero,sens.z_xy[:,end])
            @test all(iszero,sens.z_yx[:,end])
        end
        air = zeros(size(rho)); air[1:2,:] .= 1
        @test all(k -> all(iszero,getproperty(ApplyFrechet2D(mesh,rho,air),k)),fields)
        @test all(iszero,ApplyFrechetTranspose2D(mesh,rho,(;)))
        @test_throws ArgumentError run_mt2d_forward(mesh,rho;mode=:bad)
        @test_throws ArgumentError ApplyFrechet2D(mesh,rho,direction;parameterization=:bad)
        @test_throws DimensionMismatch ApplyFrechet2D(mesh,rho,zeros(2,2))
        @test_throws DimensionMismatch FrechetDerivative2D(mesh,rho;active_cells=falses(2,2))
        bad = copy(rho); bad[3,2] = -1
        @test_throws ArgumentError run_mt2d_forward(mesh,bad)
    end
end

# TM surface fields next to air involve tiny currents times ρ_air, so the finite-difference
# check uses 1e6 ohm m air to keep its noise below the tolerance
@testset "2D Fréchet derivatives with topography" begin
    mesh = _remesh(_hill_mesh(); air_resistivity = 1e6)
    air = mt2d_air_mask(mesh)
    rng = MersenneTwister(4)
    ρ = build_mt2d_halfspace_model(mesh)
    ρ[3:end, :] .*= exp.(0.3randn(rng, size(ρ, 1) - 2, size(ρ, 2)))
    fields = (:rho_xy, :phase_xy, :z_xy, :rho_yx, :phase_yx, :z_yx)
    response = run_mt2d_forward(mesh, ρ)
    @test all(k -> all(isfinite, getproperty(response, k)), fields)
    direction = randn(rng, size(ρ))
    h = 1e-4
    for mode in (:TE, :TM, :TETM)
        tangent = ApplyFrechet2D(mesh, ρ, direction; mode, parameterization = :log_resistivity)
        rp = run_mt2d_forward(mesh, ρ .* exp.(h * direction); mode)
        rm = run_mt2d_forward(mesh, ρ .* exp.(-h * direction); mode)
        for key in fields
            fd = (getproperty(rp, key) - getproperty(rm, key)) / (2h)
            @test getproperty(tangent, key) ≈ fd rtol = 2e-5 atol = 1e-10
        end
        weights = NamedTuple{fields}((randn(rng, 2, 6), randn(rng, 2, 6), randn(rng, ComplexF64, 2, 6),
                                      randn(rng, 2, 6), randn(rng, 2, 6), randn(rng, ComplexF64, 2, 6)))
        gradient = ApplyFrechetTranspose2D(mesh, ρ, weights; mode, parameterization = :log_resistivity)
        @test dot(gradient, direction) ≈ sum(real(dot(getproperty(weights, k), getproperty(tangent, k))) for k in fields) rtol = 1e-8
        @test all(iszero, gradient[air])
    end
    hilltop_air = CartesianIndex(3, findfirst(>(0), mesh.topo_air))
    G = FrechetDerivative2D(mesh, ρ; active_cells = [hilltop_air, CartesianIndex(6, 4)])
    @test all(iszero, G.z_xy[:, 1]) && all(iszero, G.z_yx[:, 1]) && any(!iszero, G.z_xy[:, 2])
    @test FrechetDerivative2D(mesh, ρ; active_cells = findall(.!air)[1:2]).cells == findall(.!air)[1:2]
end

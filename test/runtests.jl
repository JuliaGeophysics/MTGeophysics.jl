using Test
using MTGeophysics

@testset "MTGeophysics.jl Tests" begin
    
    @testset "Data Structures" begin
        # Test ModEMData creation
        data = make_nan_data()
        @test isa(data, ModEMData)
        @test data.ns == 0  # Initially no sites
        @test data.nf == 0  # Initially no frequencies
    end

    @testset "MT 1D Forward Modeling" begin
        # Test parameters
        frequencies = [1.0, 10.0, 100.0]
        mesh = [0.0, 1000.0, 2000.0] 
        resistivities = [100.0, 1000.0]
        
        rho_a, phi = MT1D(frequencies, resistivities, mesh, Analytical())
        @test length(rho_a) == length(frequencies)
        @test length(phi) == length(frequencies)
        @test all(rho_a .> 0)
        @test all(-90 .< phi .< 90)

        rho_fd, phi_fd = MT1D(frequencies, resistivities, mesh, FiniteDifference())
        @test length(rho_fd) == length(frequencies)
        @test length(phi_fd) == length(frequencies)
        @test all(rho_fd .> 0)
        @test rho_fd ≈ rho_a rtol=0.2
        @test phi_fd ≈ phi rtol=0.2
    end
end
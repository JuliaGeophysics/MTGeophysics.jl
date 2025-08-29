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
        
        # Test analytical solution
        rho_a, phi = MT1D_response(frequencies, mesh, resistivities)
        @test length(rho_a) == length(frequencies)
        @test length(phi) == length(frequencies)
        @test all(rho_a .> 0)  # Apparent resistivity should be positive
        @test all(-90 .< phi .< 90)  # Phase should be reasonable
        
        # Test FDD solution
        rho_a_fdd, phi_fdd = MT1D_response_FDD(frequencies, mesh, resistivities)
        @test length(rho_a_fdd) == length(frequencies)
        @test length(phi_fdd) == length(frequencies) 
        @test all(rho_a_fdd .> 0)  # Apparent resistivity should be positive
    end
    
    @testset "UBC Format" begin
        # Test UBC file generation (in temporary directory)
        mktempdir() do tmpdir
            cd(tmpdir)
            boundaries = [0.0, 100.0, 200.0]
            resistivities = [50.0, 100.0]
            
            # Should not throw error
            UBC_1D(boundaries, resistivities)
            
            # Check files were created
            @test isfile("model.mesh")
            @test isfile("model.rho")
        end
    end
end
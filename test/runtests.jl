# Package smoke tests.
# Author: @pankajkmishra
# This file checks core constructors and basic forward-model outputs for expected behavior.
# It runs both the original 3D data-structure tests and the 1D/2D regression checks.

using Test
using MTGeophysics

include(joinpath(dirname(@__DIR__), "Helpers", "benchmarks_1d.jl"))
include(joinpath(dirname(@__DIR__), "Helpers", "benchmarks_2d.jl"))

@testset "MTGeophysics.jl Tests" begin

    @testset "Data Structures" begin
        data = make_nan_data()
        @test isa(data, Data)
        @test data.ns == 0
        @test data.nf == 0
    end

    include(joinpath(@__DIR__, "TestForward1D.jl"))
    include(joinpath(@__DIR__, "TestForward2D.jl"))
    include(joinpath(@__DIR__, "TestIO2D.jl"))
end

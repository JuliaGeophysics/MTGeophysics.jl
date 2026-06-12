# Package smoke tests.
# Author: @pankajkmishra
# This file checks core constructors and basic forward-model outputs for expected behavior.
# It runs both the original 3D data-structure tests and the 1D/2D regression checks.

using Test
using MTGeophysics

include(joinpath(dirname(@__DIR__), "helpers", "benchmarks_1d.jl"))
include(joinpath(dirname(@__DIR__), "helpers", "benchmarks_2d.jl"))

@testset "MTGeophysics.jl Tests" begin

    include(joinpath(@__DIR__, "TestIO3D.jl"))
    include(joinpath(@__DIR__, "TestCore3D.jl"))
    include(joinpath(@__DIR__, "TestShapefileOverlay.jl"))
    include(joinpath(@__DIR__, "TestForward1D.jl"))
    include(joinpath(@__DIR__, "TestForward2D.jl"))
    include(joinpath(@__DIR__, "TestIO2D.jl"))
end

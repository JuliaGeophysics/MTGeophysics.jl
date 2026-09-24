# Package tests
# Author: @pankajkmishra
# Runs every test file: 3D data, masks and VFSA helpers, phase tensors and overlays, then 1D and 2D from the forward
# solver up to the file workflows; the benchmark writers and the shared 2D fixtures are loaded first

using Test
using MTGeophysics

include(joinpath(dirname(@__DIR__), "helpers", "benchmarks_1D.jl"))
include(joinpath(dirname(@__DIR__), "helpers", "benchmarks_2D.jl"))
include(joinpath(@__DIR__, "Helpers2D.jl"))

@testset "MTGeophysics.jl Tests" begin

    #---------- 3D ----------
    include(joinpath(@__DIR__, "TestIO3D.jl"))
    include(joinpath(@__DIR__, "TestRotate.jl"))
    include(joinpath(@__DIR__, "TestDistortion3D.jl"))
    include(joinpath(@__DIR__, "TestCore3D.jl"))
    include(joinpath(@__DIR__, "TestMask3D.jl"))
    include(joinpath(@__DIR__, "TestRBFPadding3D.jl"))
    include(joinpath(@__DIR__, "TestPhaseTensor.jl"))
    include(joinpath(@__DIR__, "TestShapefileOverlay.jl"))

    #---------- 1D ----------
    include(joinpath(@__DIR__, "TestForward1D.jl"))
    include(joinpath(@__DIR__, "TestInversion1D.jl"))

    #---------- 2D ----------
    include(joinpath(@__DIR__, "TestForward2D.jl"))
    include(joinpath(@__DIR__, "TestFrechet2D.jl"))
    include(joinpath(@__DIR__, "TestMesh2D.jl"))
    include(joinpath(@__DIR__, "TestModelFile2D.jl"))
    include(joinpath(@__DIR__, "TestDataFile2D.jl"))
    include(joinpath(@__DIR__, "TestControl2D.jl"))
    include(joinpath(@__DIR__, "TestStrike2D.jl"))
    include(joinpath(@__DIR__, "TestTopography2D.jl"))
    include(joinpath(@__DIR__, "TestBenchmarks2D.jl"))
    include(joinpath(@__DIR__, "TestInversion2D.jl"))
    include(joinpath(@__DIR__, "TestVFSA2D.jl"))
    include(joinpath(@__DIR__, "TestMakeMesh2D.jl"))
    include(joinpath(@__DIR__, "TestWorkflows2D.jl"))
end

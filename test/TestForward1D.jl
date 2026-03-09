# This script checks that the 1D analytical and finite-difference solvers stay aligned on a simple layered-earth example.

using Test

frequencies = collect(10 .^ range(-2, 2, length = 12))
thicknesses = [150.0, 350.0, 900.0]
resistivities = [100.0, 20.0, 500.0, 1000.0]

result = Forward1D(
    frequencies,
    resistivities,
    thicknesses;
    first_cell = 15.0,
    growth_factor = 1.08,
    padding_cells = 20,
)

analytical = result.responses[:analytical]
fd = result.responses[:fd]

@test maximum(abs.((fd.apparent_resistivity .- analytical.apparent_resistivity) ./ analytical.apparent_resistivity)) < 0.30
@test maximum(abs.(fd.phase .- analytical.phase)) < 12.0

mktempdir() do temp_dir
    mesh_result = MakeMesh1D(output_dir = temp_dir)
    data_path = write_mt1d_data_template(joinpath(temp_dir, "Data.dat"), frequencies)
    observed_path = ForwardSolve1D(mesh_result.model_path, data_path)

    @test basename(mesh_result.model_path) == "Layered1D.true"
    @test basename(observed_path) == "Data.obs"
    @test isfile(mesh_result.model_path)
    @test isfile(observed_path)
    @test length(load_mt1d_data_spec(observed_path).frequencies) == length(frequencies)
    @test isfile(PlotData1D(observed_path, observed_path; output_path = joinpath(temp_dir, "DataCompare1D.png")))
    @test isfile(PlotModel1D(mesh_result.model_path, observed_path; output_path = joinpath(temp_dir, "ModelDepth1D.png"), maximum_depth_m = 500.0))

    benchmarks = SaveBenchmarks1D(output_root = temp_dir)
    @test isfile(benchmarks[1].model_path)
    @test isfile(benchmarks[1].reference_path)
    @test isfile(benchmarks[1].observed_path)
    @test isdir(benchmarks[1].case_dir)
end

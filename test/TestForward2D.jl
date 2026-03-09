# This script checks that the 2D forward solver returns stable TE and TM responses on a simple half-space profile.

using Test

mesh = BuildMesh2D(
    frequencies = collect(10 .^ range(-1, 1, length = 4)),
    y_core_range = (-2000.0, 2000.0),
    y_core_cell = 500.0,
    y_padding = 2500.0,
    pad_factor = 1.20,
    air_top = -3000.0,
    air_cells = 4,
    ground_layers = [200.0, 200.0, 400.0, 400.0, 800.0, 800.0],
    receiver_stride = 2,
)

resistivity = MTGeophysics.build_mt2d_halfspace_model(mesh; background_resistivity = 100.0)
response = run_mt2d_forward(mesh, resistivity)
station_index = MTGeophysics.mt2d_center_station(mesh)

benchmark_mesh = BuildMesh2D(
    frequencies = [1.0],
    y_core_range = (-10_000.0, 10_000.0),
    y_core_cell = 200.0,
    y_padding = 20_000.0,
    pad_factor = 1.30,
    air_top = -30_000.0,
    air_cells = 10,
    ground_layers = vcat(fill(100.0, 10), fill(200.0, 15), fill(500.0, 20), fill(1000.0, 25)),
    receiver_positions = collect(-10_000.0:400.0:10_000.0),
)
benchmark_models = MTGeophysics.build_mt2d_comemi_models(benchmark_mesh)
geometry_checks = MTGeophysics.validate_mt2d_comemi_models(benchmark_mesh, benchmark_models)

@test size(response.rho_xy) == (length(mesh.frequencies), length(mesh.receiver_positions))
@test size(response.rho_yx) == (length(mesh.frequencies), length(mesh.receiver_positions))
@test all(isfinite, response.rho_xy)
@test all(isfinite, response.rho_yx)
@test all(isfinite, response.phase_xy)
@test all(isfinite, response.phase_yx)
@test minimum(response.rho_xy) > 0
@test minimum(response.rho_yx) > 0
@test maximum(abs.((response.rho_xy[:, station_index] .- response.rho_yx[:, station_index]) ./ response.rho_xy[:, station_index])) < 0.40
@test maximum(abs.(response.phase_xy[:, station_index] .- MTGeophysics._phase_fold_to_0_90(response.phase_yx[:, station_index]))) < 35.0
@test length(benchmark_models) == 3
@test length(geometry_checks) == 12

mktempdir() do temp_dir
    mesh_result = MakeMesh2D(output_dir = temp_dir)
    @test basename(mesh_result.model_paths["comemi2d_case1_dyke"]) == "Comemi2D1.true"
    @test basename(mesh_result.model_paths["comemi2d_case2_resistive_blocks"]) == "Comemi2D2.true"
    @test basename(mesh_result.model_paths["comemi2d_case3_mixed"]) == "Comemi2D3.true"

    benchmarks = SaveBenchmarks2D(output_root = temp_dir)
    @test isfile(benchmarks[1].model_path)
    @test isfile(benchmarks[1].start_model_path)
    @test isfile(benchmarks[1].reference_path)
    @test isfile(benchmarks[1].observed_path)
    @test isdir(benchmarks[1].case_dir)
end

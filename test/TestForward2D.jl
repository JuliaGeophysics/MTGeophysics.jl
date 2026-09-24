# This script checks that the 2D forward solver returns stable TE and TM responses on a simple half-space profile.

using Test

# resistivity at a set of points inside each COMEMI anomaly and host
function check_comemi_models(mesh, models)
    lookup = Dict(m.name => m.resistivity for m in models)
    checks = [
        ("comemi2d_case1_dyke", 0.0, 1000.0, 5.0),
        ("comemi2d_case1_dyke", 4000.0, 1000.0, 100.0),
        ("comemi2d_case1_dyke", 4000.0, 3000.0, 500.0),
        ("comemi2d_case2_resistive_blocks", -2000.0, 1200.0, 800.0),
        ("comemi2d_case2_resistive_blocks", 4000.0, 3000.0, 400.0),
        ("comemi2d_case2_resistive_blocks", 0.0, 1000.0, 30.0),
        ("comemi2d_case2_resistive_blocks", 0.0, 3000.0, 100.0),
        ("comemi2d_case3_mixed", -4000.0, 1000.0, 3.0),
        ("comemi2d_case3_mixed", 3000.0, 5000.0, 1000.0),
        ("comemi2d_case3_mixed", 0.0, 500.0, 80.0),
        ("comemi2d_case3_mixed", 0.0, 2000.0, 20.0),
        ("comemi2d_case3_mixed", 0.0, 8000.0, 300.0),
    ]
    map(checks) do (name, y, z, expected)
        iy = searchsortedlast(mesh.y_nodes, y)
        iz = searchsortedlast(mesh.z_nodes, z)
        lookup[name][iz, iy] ≈ expected
    end
end

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
station_index = cld(length(mesh.receiver_positions), 2)

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
geometry_checks = check_comemi_models(benchmark_mesh, benchmark_models)

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
@test length(geometry_checks) == 12 && all(geometry_checks)

mktempdir() do temp_dir
    r = only(SaveBenchmarks2D(output_root = temp_dir, cases = ["2D-I"]))
    @test r.case_dir == joinpath(temp_dir, "2D-I")
    @test all(isfile, (r.true_model_path, r.data_path, r.start_model_path, r.prior_model_path, r.cov_path, r.mask_path))
    @test readdir(r.case_dir) == ["cov.ctrl", "data.dat", "mask.ctrl", "model.prior", "model.start", "model.true"]   # controls ship in examples/ctrl/2D
    @test !occursin('#', read(r.cov_path, String))
    observed = load_data2d(r.data_path)
    @test observed.site_names[1] == "JYV001" && issorted(observed.longitudes)
    @test all(isapprox.(observed.latitudes, 62.25; atol = 0.01)) && observed.origin == [62.25, 25.75]
    @test size(ReadModel2D(r.start_model_path).resistivity) == size(ReadCov2D(r.cov_path).mask)
    @test length(ReadModel2D(r.true_model_path).resistivity) > length(ReadModel2D(r.start_model_path).resistivity)
end

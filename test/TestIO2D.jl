# This script checks the 2D model/data round-trip and runs a one-iteration VFSA smoke test.

using Test

mesh = build_default_mt2d_mesh()
true_resistivity = only(filter(model -> model.name == "comemi2d_case1_dyke", MTGeophysics.build_mt2d_comemi_models(mesh))).resistivity
true_response = run_mt2d_forward(mesh, true_resistivity)

mktempdir() do temp_dir
    model_path = write_model2d(joinpath(temp_dir, "model.rho"), mesh, true_resistivity; title = "roundtrip model")
    data_path = write_data2d(
        joinpath(temp_dir, "data.dat"),
        data_from_response2d(true_response; impedance_error_fraction = 0.05, title = "roundtrip data"),
    )

    loaded_model = load_model2d(model_path)
    loaded_data = load_data2d(data_path)
    roundtrip_mesh = build_mesh_from_model2d(
        loaded_model;
        frequencies = loaded_data.frequencies,
        receiver_positions = loaded_data.receivers,
    )
    roundtrip_response = run_mt2d_forward(roundtrip_mesh, loaded_model.resistivity)
    roundtrip_predicted = data_from_response2d(
        roundtrip_response;
        z_xy_error = loaded_data.z_xy_error,
        z_yx_error = loaded_data.z_yx_error,
        z_xx_error = loaded_data.z_xx_error,
        z_yy_error = loaded_data.z_yy_error,
        site_names = loaded_data.site_names,
        x_positions = loaded_data.x_positions,
        z_positions = loaded_data.z_positions,
    )
    fit = chi2_rms2d(loaded_data, roundtrip_predicted)

    @test size(loaded_model.resistivity) == size(true_resistivity)
    @test loaded_data.receivers == true_response.receivers
    @test fit.rms < 1e-5
    @test isfile(PlotModel2D(model_path; output_path = joinpath(temp_dir, "ModelPlot2DPadding.png"), show_padding = true, maximum_depth_km = Inf))

    mesh_result = MakeMesh2D(output_dir = temp_dir)
    template_path = write_mt2d_data_template(joinpath(temp_dir, "Data.dat"), mesh_result.mesh)
    observed_path = ForwardSolve2D(mesh_result.model_paths["comemi2d_case1_dyke"], template_path)
    forward_data = load_data2d(observed_path)

    @test basename(mesh_result.model_paths["comemi2d_case1_dyke"]) == "Comemi2D1.true"
    @test basename(observed_path) == "Data.obs"
    @test isfile(observed_path)
    @test maximum(abs.(forward_data.z_xy_error .- 0.05 .* abs.(forward_data.z_xy))) < 1e-8
    @test maximum(abs.(forward_data.z_yx_error .- 0.05 .* abs.(forward_data.z_yx))) < 1e-8

    benchmarks = SaveBenchmarks2D(output_root = temp_dir)
    @test isfile(benchmarks[1].model_path)
    @test isfile(benchmarks[1].start_model_path)
    @test isfile(benchmarks[1].reference_path)
    @test isfile(benchmarks[1].observed_path)
    @test isdir(benchmarks[1].case_dir)

    smoke_config = VFSA2DMTConfig(
        n_ctrl = 8,
        n_chains = 2,
        max_iter = 1,
        n_trials = 1,
        keep_models = true,
        perturb_depth_m = 600.0,
        snapshot_interval = 1,
        output_root = temp_dir,
    )
    smoke_run = joinpath(temp_dir, "VFSA2DMT_Test")
    inversion = VFSA2DMT(
        benchmarks[1].start_model_path,
        benchmarks[1].observed_path;
        run_dir = smoke_run,
        config = smoke_config,
    )

    @test isfile(joinpath(smoke_run, "chain_01", "0vfsa2DMT.log"))
    @test isfile(joinpath(smoke_run, "chain_01", "0vfsa2DMT.details"))
    @test isfile(joinpath(smoke_run, "chain_02", "0vfsa2DMT.log"))
    @test isfile(joinpath(smoke_run, "chain_02", "0vfsa2DMT.details"))
    @test isdir(joinpath(smoke_run, "snapshots"))
    @test isfile(joinpath(smoke_run, "snapshots", "avg_iter_00001.rho"))
    @test isfile(joinpath(smoke_run, "chain_01", "best_iter_00001.rho"))
    @test isfile(joinpath(smoke_run, "chain_02", "best_iter_00001.rho"))
    @test isfile(joinpath(smoke_run, "data.obs"))
    @test isfile(joinpath(smoke_run, "model.true"))
    @test isfile(joinpath(smoke_run, "chain_01", "itr_00001_trial_01.rho"))
    @test isfile(joinpath(smoke_run, "model.c1best"))
    @test !isfile(joinpath(smoke_run, "chain_01_final.rho"))
    @test isfile(joinpath(smoke_run, "data.c1best"))
    @test !isfile(joinpath(smoke_run, "dpred_01_00001_01.pred"))
    @test isfile(joinpath(smoke_run, "model.mean"))
    @test isfile(joinpath(smoke_run, "model.median"))
    @test isfile(joinpath(smoke_run, "model.std"))
    @test isfile(joinpath(smoke_run, "Summary.md"))
    @test length(inversion.best_chain.iterations) == 1
    @test inversion.best_chain.best_chi2 > 0
    @test inversion.ensemble.chain_count == 2
    @test inversion.chains[1].final_model_path == ""
    @test maximum(inversion.chains[1].control_z_m) <= 600.0 + 1e-9

    stats = AnalyseEnsemble2D(smoke_run; copy_script = false)
    @test isfile(stats.summary_path)
    @test isfile(stats.mean_model_path)

    zero_iter_run = MTGeophysics.run_mt2d_vfsa(
        benchmarks[1].start_model_path,
        benchmarks[1].observed_path;
        run_dir = joinpath(temp_dir, "VFSA2DMT_ZeroIter"),
        config = VFSA2DMTConfig(
            n_ctrl = 8,
            n_chains = 1,
            max_iter = 0,
            n_trials = 1,
            perturb_depth_m = 600.0,
            keep_models = false,
            output_root = temp_dir,
        ),
    )
    zero_iter_start = load_model2d(benchmarks[1].start_model_path)
    zero_iter_data = load_data2d(zero_iter_run.observed_data_path)
    zero_iter_mesh = build_mesh_from_model2d(
        zero_iter_start;
        frequencies = zero_iter_data.frequencies,
        receiver_positions = zero_iter_data.receivers,
    )
    zero_iter_ground = zero_iter_run.chains[1].final_resistivity[(zero_iter_mesh.n_air_cells + 1):end, :]
    perturbable_z = MTGeophysics._perturbable_ground_indices(zero_iter_mesh, 600.0)
    @test all(isapprox.(zero_iter_ground, 100.0; atol = 1e-10, rtol = 0.0))
    @test last(perturbable_z) < length(zero_iter_mesh.z_cell_sizes)

    example_run = VFSA2DMT(
        VFSA2DMTParams(
            script_path = joinpath(temp_dir, \"run_vfsa2dmt.jl\"),
            data_path = benchmarks[1].observed_path,
            start_model_path = benchmarks[1].start_model_path,
            model_path = benchmarks[1].model_path,
            config = VFSA2DMTConfig(
                n_ctrl = 4,
                n_chains = 1,
                max_iter = 0,
                n_trials = 1,
                perturb_depth_m = 600.0,
                keep_models = false,
            ),
        ),
    )
    @test dirname(example_run.run_info.run_dir) == temp_dir
    @test startswith(basename(example_run.run_info.run_dir), "run_VFSA2DMT_")
end

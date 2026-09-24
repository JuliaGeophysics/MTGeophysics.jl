# This script checks the older 2D model layout round trip and runs an in-memory VFSA smoke test.

using Test

mesh = build_default_mt2d_mesh()
true_resistivity = only(filter(model -> model.name == "comemi2d_case1_dyke", MTGeophysics.build_mt2d_comemi_models(mesh))).resistivity
true_response = run_mt2d_forward(mesh, true_resistivity)

mktempdir() do temp_dir
    model_path = MTGeophysics._write_model2d_legacy(joinpath(temp_dir, "model.rho"), mesh, true_resistivity; title = "roundtrip model")
    data_path = write_data2d(
        joinpath(temp_dir, "data.dat"),
        data_from_response2d(true_response; impedance_error_fraction = 0.05, title = "roundtrip data"),
    )

    loaded_model = MTGeophysics._load_model2d_legacy(model_path)
    loaded_data = load_data2d(data_path)
    roundtrip_mesh = MTGeophysics._mesh_from_legacy_model2d(
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

    # vfsa on the noisy dyke data, from a halfspace
    observed = MTGeophysics._apply_mt2d_noise(data_from_response2d(true_response; impedance_error_fraction = 0.05);
                                              rng_seed = 20260308)
    start = build_mt2d_halfspace_model(mesh)
    smoke_config = VFSA2DConfig(n_ctrl = 8, n_chains = 2, max_iter = 2, perturb_depth_m = 600.0,
                                snapshot_interval = 1, verbose = false)
    smoke_run = joinpath(temp_dir, "VFSA2D_Test")
    inversion = VFSA2D(mesh, start, observed; config = smoke_config, run_dir = smoke_run)
    vdir = joinpath(smoke_run, "vfsa")
    for k in 1:2, f in ("best.rho", "History.csv", "best_iter_00001.rho")
        @test isfile(joinpath(vdir, "chain_0$k", f))
    end
    for f in ("model.mean.rho", "model.median.rho", "model.p05.rho", "model.p95.rho", "model.best.rho",
              "Uncertainty.csv", "Chains.csv")
        @test isfile(joinpath(vdir, f))
    end
    @test length(inversion.chains) == 2 && all(c -> length(c.history) <= 3, inversion.chains)
    @test inversion.ensemble.count == 2 && all(isfinite, inversion.ensemble.std)
    @test isfinite(inversion.rms) && inversion.best_rms <= minimum(c.history[1].rms for c in inversion.chains)
    @test all(c -> maximum(mesh.z_nodes[k[1]+1] for k in c.controls) <= 600.0 + 1e-6, inversion.chains)
    @test countlines(joinpath(vdir, "Uncertainty.csv")) == 1 + count(.!mt2d_air_mask(mesh))
    stats = AnalyseEnsemble2D(smoke_run)
    @test length(stats.chains) == 2 && all(isfile, stats.paths)

    # no iterations: the ensemble is the start model
    zero_run = VFSA2D(mesh, start, observed; config = VFSA2DConfig(n_ctrl = 8, n_chains = 1, max_iter = 0, verbose = false))
    na = mesh.n_air_cells
    @test all(isapprox.(zero_run.resistivity[na+1:end, :], start[na+1:end, :]; rtol = 1e-9))
end

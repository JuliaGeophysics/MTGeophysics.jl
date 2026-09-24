# 2D VFSA
# Author: @pankajkmishra
# Ensures the in-memory VFSA runs its chains, writes chains, snapshots and ensemble, keeps its controls in the core,
# carries the core bottom down a third per layer, rebuilds the ensemble from the chains, and returns the start
# model after zero iterations

using Test

@testset "2D VFSA" begin
    mesh = build_default_mt2d_mesh()
    true_resistivity = only(filter(model -> model.name == "comemi2d_case1_dyke", MTGeophysics.build_mt2d_comemi_models(mesh))).resistivity
    true_response = run_mt2d_forward(mesh, true_resistivity)
    mktempdir() do temp_dir
        # vfsa on the noisy dyke data, from a halfspace
        observed = MTGeophysics._apply_mt2d_noise(data_from_response2d(true_response; impedance_error_fraction = 0.05);
                                                  rng_seed = 20260308)
        start = build_mt2d_halfspace_model(mesh)
        smoke_config = VFSA2DConfig(n_ctrl = 8, n_chains = 2, max_iter = 2, z_core_cells = 3,
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
        na0 = mesh.n_air_cells
        core = MTGeophysics._core_range(mesh.y_cell_sizes)
        @test all(c -> all(k -> na0 < k[1] <= na0 + 3 && k[2] in core, c.controls), inversion.chains)
        # below the core each column carries the core bottom down a third per layer, towards the start
        b = inversion.chains[1].best.m
        j = first(core) + 2
        @test b[na0+4, j] - log10(start[na0+4, j]) ≈ (b[na0+3, j] - log10(start[na0+3, j])) / 3 atol = 1e-9
        @test countlines(joinpath(vdir, "Uncertainty.csv")) == 1 + count(.!mt2d_air_mask(mesh))
        stats = AnalyseEnsemble2D(smoke_run)
        @test length(stats.chains) == 2 && all(isfile, stats.paths)

        # no iterations: the ensemble is the start model
        zero_run = VFSA2D(mesh, start, observed; config = VFSA2DConfig(n_ctrl = 8, n_chains = 1, max_iter = 0, verbose = false))
        na = mesh.n_air_cells
        @test all(isapprox.(zero_run.resistivity[na+1:end, :], start[na+1:end, :]; rtol = 1e-9))
    end
end

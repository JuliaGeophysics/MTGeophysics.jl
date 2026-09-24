# 1D forward, Fréchet derivatives, MakeMesh1D and Invert1D, 1D's own code; the 2D solver only as a cross-check
# Author: @pankajkmishra

using Test, LinearAlgebra, Random

@testset "1D" begin
    f = collect(10 .^ range(-2, 2, length = 12))

    @testset "layered forward" begin
        Z = mt1d_impedance(f, [100.0], Float64[])
        @test abs2.(Z) ./ (2π .* f .* 4π * 1e-7) ≈ fill(100.0, 12)
        @test all(rad2deg.(angle.(Z)) .≈ 45)
        # two layers: the top at high frequency, the basement at low frequency
        Z2 = mt1d_impedance([1e4, 1e-4], [10.0, 1000.0], [500.0])
        ρa = abs2.(Z2) ./ (2π .* [1e4, 1e-4] .* 4π * 1e-7)
        @test ρa[1] ≈ 10 rtol = 1e-3
        # closed-form two-layer impedance; 50 S of cover still pulls the lowest frequency below 1000
        two(ω) = (k = ρ -> sqrt(1im * ω * 4π * 1e-7 / ρ); Zi = ρ -> 1im * ω * 4π * 1e-7 / k(ρ);
                  t = tanh(k(10.0) * 500.0); Zi(10.0) * (Zi(1000.0) + Zi(10.0) * t) / (Zi(10.0) + Zi(1000.0) * t))
        @test ρa ≈ abs2.(two.(2π .* [1e4, 1e-4])) ./ (2π .* [1e4, 1e-4] .* 4π * 1e-7) rtol = 1e-8
        @test 900 < ρa[2] < 1000

        # the 2D finite-difference solver on a laterally uniform model agrees
        mesh2 = BuildMesh2D(frequencies = f, y_core_range = (-2000.0, 2000.0), y_core_cell = 500.0, y_padding = 20_000.0,
                            air_top = -50_000.0, air_cells = 10, receiver_positions = [0.0],
                            ground_layers = mt2d_geometric_layers(f; first_layer_div = 10.0, vertical_factor = 1.05))
        layers = [(0.0, 150.0, 100.0), (150.0, 500.0, 20.0), (500.0, Inf, 500.0)]
        ρ2 = build_mt2d_halfspace_model(mesh2)
        zc = MTGeophysics.mt2d_z_centers(mesh2)
        for iz in mesh2.n_air_cells+1:length(zc), (top, bottom, ρ) in layers
            top <= zc[iz] < bottom && (ρ2[iz, :] .= ρ)
        end
        r2 = run_mt2d_forward(mesh2, ρ2)
        Z1 = mt1d_impedance(f, ρ2[mesh2.n_air_cells+1:end, 1], mesh2.z_cell_sizes[mesh2.n_air_cells+1:end])
        @test maximum(abs.(abs2.(Z1) ./ (2π .* f .* 4π * 1e-7) ./ r2.rho_xy[:, 1] .- 1)) < 0.03
        @test maximum(abs.(rad2deg.(angle.(Z1)) .- r2.phase_xy[:, 1])) < 1.0
    end

    @testset "Fréchet derivatives" begin
        rng = MersenneTwister(3)
        h = fill(100.0, 12)
        ρ = 100 .* exp.(0.5randn(rng, 12))
        δ = randn(rng, 12)
        G = mt1d_frechet(f, ρ, h)
        ε = 1e-6
        fd = (mt1d_impedance(f, ρ .* 10 .^ (ε .* δ), h) .- mt1d_impedance(f, ρ .* 10 .^ (-ε .* δ), h)) ./ 2ε
        @test G * δ ≈ fd rtol = 1e-6
        @test size(G) == (12, 12)
        @test sum(mt1d_layers(f; background_resistivity = 100.0)) >= 4 * mt1d_skin_depth(100.0, minimum(f))
    end

    @testset "files and inversion" begin
        mktempdir() do dir
            r = only(SaveBenchmarks1D(output_root = dir))
            @test sort(readdir(r.case_dir)) == ["data.dat", "model.true"]
            ctrl = joinpath(dirname(@__DIR__), "examples", "ctrl", "1D")
            observed = load_data2d(r.data_path)
            @test observed.site_names == ["Fin001"]
            @test observed.z_yx ≈ -observed.z_xy rtol = 0.5

            # forward from the true model reproduces the noise-free data, and writes G on request
            pred = ForwardSolve1D(r.true_model_path, r.data_path; mode = :DET, write_frechet = true,
                                  output_path = joinpath(dir, "true.pred"))
            @test isfile(joinpath(dir, "true.frechet"))
            p = load_data2d(pred)
            @test chi2_rms2d(observed, p).rms < 1.6
            @test isfile(PlotModel1D(r.true_model_path; output_path = joinpath(dir, "True.png")))

            # the mesh follows the site's own skin depths
            m = only(MakeMesh1D(observed))
            @test m.site == "Fin001" && 30 < m.background < 1000
            @test m.thicknesses[1] ≈ mt1d_skin_depth(m.background, maximum(observed.frequencies)) / 5
            @test sum(m.thicknesses) >= 4 * mt1d_skin_depth(m.background, minimum(observed.frequencies)) - m.thicknesses[end]

            # control files: one reader, each algorithm with its own keys
            gc = InvCtrl1D(algorithm = :gn, mode = :DET, target_rms = 1.0, max_iter = 20, lambda = 1.0)
            @test ReadInvCtrl1D(WriteInvCtrl1D(joinpath(dir, "gn.ctrl"), gc)) == gc
            vc = InvCtrl1D(algorithm = :vfsa, target_rms = 1.0, max_iter = 200, chains = 2)
            @test ReadInvCtrl1D(WriteInvCtrl1D(joinpath(dir, "vfsa.ctrl"), vc)) == vc
            @test !occursin("lambda", read(joinpath(dir, "vfsa.ctrl"), String))
            write(joinpath(dir, "bad.ctrl"), read(joinpath(ctrl, "InvCtrl.GN"), String) * "Number of chains : 2\n")
            @test_throws ErrorException ReadInvCtrl1D(joinpath(dir, "bad.ctrl"))
            write(joinpath(dir, "bad.ctrl"), read(joinpath(ctrl, "InvCtrl.GN"), String) * "Log10 resistivity bounds : 0 4\n")
            @test_throws ErrorException ReadInvCtrl1D(joinpath(dir, "bad.ctrl"))
            @test !occursin("bounds", read(joinpath(dir, "gn.ctrl"), String)) && occursin("bounds", read(joinpath(dir, "vfsa.ctrl"), String))
            write(joinpath(dir, "bad.ctrl"), replace(read(joinpath(ctrl, "InvCtrl.GN"), String), "GN" => "NLCG"))
            @test_throws ErrorException ReadInvCtrl1D(joinpath(dir, "bad.ctrl"))

            gn = Invert1D(r.data_path, joinpath(ctrl, "InvCtrl.GN"); run_dir = joinpath(dir, "gn"))
            @test gn.rms < 1.2 && gn.sites[1].converged
            @test all(isfile, joinpath.(dir, "gn", ["data.pred", "Summary.txt", "inputs/InvCtrl.GN", "Fin001/model.start",
                                                    "Fin001/model.rho", "Fin001/History.csv"]))
            @test all(isfile, PlotInversion1D(gn; true_model_path = r.true_model_path))

            v = Invert1D(r.data_path, joinpath(dir, "vfsa.ctrl"), MakeMesh1D(observed); run_dir = joinpath(dir, "vfsa"))
            @test isfinite(v.rms) && length(v.sites[1].vfsa.chains) == 2 && v.algorithm == :vfsa
            @test all(isfile, joinpath.(dir, "vfsa", "Fin001", "vfsa", ["model.mean.rho", "model.p05.rho", "History_chain_01.csv"]))
            @test v.sites[1].vfsa.best_rms < v.sites[1].vfsa.chains[1].history[1].rms
            @test isfile(joinpath(dir, "vfsa", "Fin001", "vfsa", "data.best.pred"))
            @test all(isfile, PlotInversion1D(v))

            det = Invert1D(r.data_path, joinpath(dir, "gn.ctrl"); run_dir = joinpath(dir, "det"))
            @test det.rms < 1.3
        end
    end
end

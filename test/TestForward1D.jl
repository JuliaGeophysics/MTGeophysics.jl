# 1D forward, Fréchet derivatives, MakeMesh1D and Invert1D, which reuse the 2D machinery
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
        mesh1 = Mesh1D(mesh2.z_cell_sizes[mesh2.n_air_cells+1:end], f)
        r1 = run_mt2d_forward(mesh1, ρ2[mesh2.n_air_cells+1:end, 1:1])
        @test maximum(abs.(r1.rho_xy ./ r2.rho_xy .- 1)) < 0.03
        @test maximum(abs.(r1.phase_xy .- r2.phase_xy)) < 1.0
        @test r1.z_yx ≈ -r1.z_xy
    end

    @testset "Fréchet derivatives" begin
        rng = MersenneTwister(3)
        mesh = Mesh1D(fill(100.0, 12), f)
        ρ = 100 .* exp.(0.5randn(rng, 12, 1))
        δ = randn(rng, 12, 1)
        fields = (:rho_xy, :phase_xy, :z_xy, :rho_yx, :phase_yx, :z_yx)
        h = 1e-5
        for mode in (:TE, :TM, :TETM)
            t = ApplyFrechet2D(mesh, ρ, δ; mode, parameterization = :log_resistivity)
            rp, rm = run_mt2d_forward(mesh, ρ .* exp.(h * δ); mode), run_mt2d_forward(mesh, ρ .* exp.(-h * δ); mode)
            for k in fields
                @test getproperty(t, k) ≈ (getproperty(rp, k) - getproperty(rm, k)) / 2h rtol = 1e-6 atol = 1e-10
            end
            w = NamedTuple{fields}((randn(rng, 12, 1), randn(rng, 12, 1), randn(rng, ComplexF64, 12, 1),
                                    randn(rng, 12, 1), randn(rng, 12, 1), randn(rng, ComplexF64, 12, 1)))
            g = ApplyFrechetTranspose2D(mesh, ρ, w; mode, parameterization = :log_resistivity)
            @test dot(g, δ) ≈ sum(real(dot(getproperty(w, k), getproperty(t, k))) for k in fields) rtol = 1e-10
        end
        G = FrechetDerivative2D(mesh, ρ; parameterization = :log10_resistivity)
        @test size(G.z_xy) == (12, 12)
    end

    @testset "files and inversion" begin
        mktempdir() do dir
            r = only(SaveBenchmarks1D(output_root = dir))
            @test sort(readdir(r.case_dir)) == ["data.dat", "model.true"]
            ctrl = joinpath(dirname(@__DIR__), "examples", "ctrl", "1D")
            observed = load_data2d(r.data_path)
            @test observed.site_names == ["JYV1D"]
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
            @test m.site == "JYV1D" && 30 < m.background < 1000
            @test m.thicknesses[1] ≈ mt2d_skin_depth(m.background, maximum(observed.frequencies)) / 5
            @test sum(m.thicknesses) >= 4 * mt2d_skin_depth(m.background, minimum(observed.frequencies)) - m.thicknesses[end]

            # control files: one reader, each algorithm with its own keys
            gc = InvCtrl1D(algorithm = :gn, mode = :DET, target_rms = 1.0, max_iter = 20, lambda = 1.0)
            @test ReadInvCtrl1D(WriteInvCtrl1D(joinpath(dir, "gn.ctrl"), gc)) == gc
            vc = InvCtrl1D(algorithm = :vfsa, target_rms = 1.0, max_iter = 20, chains = 2, control_points = 12)
            @test ReadInvCtrl1D(WriteInvCtrl1D(joinpath(dir, "vfsa.ctrl"), vc)) == vc
            @test !occursin("lambda", read(joinpath(dir, "vfsa.ctrl"), String))
            write(joinpath(dir, "bad.ctrl"), read(joinpath(ctrl, "InvCtrl.GN"), String) * "Number of chains : 2\n")
            @test_throws ErrorException ReadInvCtrl1D(joinpath(dir, "bad.ctrl"))
            write(joinpath(dir, "bad.ctrl"), replace(read(joinpath(ctrl, "InvCtrl.GN"), String), "GN" => "NLCG"))
            @test_throws ErrorException ReadInvCtrl1D(joinpath(dir, "bad.ctrl"))

            gn = Invert1D(r.data_path, joinpath(ctrl, "InvCtrl.GN"); run_dir = joinpath(dir, "gn"))
            @test gn.rms < 1.2 && gn.sites[1].converged
            @test all(isfile, joinpath.(dir, "gn", ["data.pred", "Summary.txt", "inputs/InvCtrl.GN", "JYV1D/model.start",
                                                    "JYV1D/model.rho", "JYV1D/History.csv"]))
            @test all(isfile, PlotInversion1D(gn; true_model_path = r.true_model_path))

            v = Invert1D(r.data_path, joinpath(dir, "vfsa.ctrl"), MakeMesh1D(observed); run_dir = joinpath(dir, "vfsa"))
            @test isfinite(v.rms) && length(v.sites[1].vfsa.chains) == 2 && v.algorithm == :vfsa
            @test isdir(joinpath(dir, "vfsa", "JYV1D", "vfsa"))
            @test all(isfile, PlotInversion1D(v))

            det = Invert1D(r.data_path, joinpath(dir, "gn.ctrl"); run_dir = joinpath(dir, "det"))
            @test det.rms < 1.3
        end
    end
end

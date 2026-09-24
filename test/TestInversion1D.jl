# 1D files, mesh and inversion
# Author: @pankajkmishra
# Ensures the 1D benchmark case, ForwardSolve1D (with G), MakeMesh1D and the InvCtrl1D files behave, and that
# Invert1D fits the benchmark with GN, VFSA and the determinant, writing its run folder and plots

using Test

@testset "1D inversion" begin
    mktempdir() do dir
        r = only(SaveBenchmarks1D(output_root = dir))
        @test sort(readdir(r.case_dir)) == ["data.dat", "model.true"]
        ctrl = joinpath(dirname(@__DIR__), "examples", "ctrl", "1D")
        observed = load_data2d(r.data_path)
        @test observed.site_names == ["UA01"]
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
        @test m.site == "UA01" && 30 < m.background < 1000
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
        @test all(isfile, joinpath.(dir, "gn", ["data.pred", "Summary.txt", "inputs/InvCtrl.GN", "UA01/model.start",
                                                "UA01/model.rho", "UA01/History.csv"]))
        @test all(isfile, PlotInversion1D(gn; true_model_path = r.true_model_path))

        v = Invert1D(r.data_path, joinpath(dir, "vfsa.ctrl"), MakeMesh1D(observed); run_dir = joinpath(dir, "vfsa"))
        @test isfinite(v.rms) && length(v.sites[1].vfsa.chains) == 2 && v.algorithm == :vfsa
        @test all(isfile, joinpath.(dir, "vfsa", "UA01", "vfsa", ["model.mean.rho", "model.p05.rho", "History_chain_01.csv"]))
        @test v.sites[1].vfsa.best_rms < v.sites[1].vfsa.chains[1].history[1].rms
        @test isfile(joinpath(dir, "vfsa", "UA01", "vfsa", "data.best.pred"))
        @test all(isfile, PlotInversion1D(v))

        det = Invert1D(r.data_path, joinpath(dir, "gn.ctrl"); run_dir = joinpath(dir, "det"))
        @test det.rms < 1.3
    end
end

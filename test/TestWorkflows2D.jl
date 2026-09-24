# 2D file workflows
# Author: @pankajkmishra
# Ensures the ModEM-style file runs work end to end: ForwardSolve2D fills a template and writes G, the six-file
# Invert2D and five-file VFSA2D fit, fix masked cells and write their run folders and plots, and with topography
# air and water stay fixed while bad masks are refused

using Test

@testset "2D file workflows" begin
    @testset "flat survey" begin
        mktempdir() do dir
            # a small survey with the shipped controls
            path(f) = joinpath(dir, f)
            ctrl_dir = joinpath(dirname(@__DIR__), "examples", "ctrl", "2D")
            cp(joinpath(ctrl_dir, "FwdCtrl"), path("fwd.ctrl"))
            cp(joinpath(ctrl_dir, "InvCtrl.GN"), path("inv.ctrl"))
            f, y = [1.0, 10.0, 100.0], [-750.0, -250.0, 250.0, 750.0]
            mesh = BuildMesh2D(frequencies = f, receiver_positions = y, y_core_range = (-1000.0, 1000.0), y_core_cell = 250.0,
                               y_padding = 4000.0, pad_factor = 1.5, air_top = -20_000.0, air_cells = 8, max_core_layers = 6)
            WriteModel2D(path("model.rho"), mesh, build_mt2d_halfspace_model(mesh; background_resistivity = 100.0))
            nf, ns = length(f), length(y)
            nan = fill(NaN, nf, ns)
            write_data2d(path("data.dat"), DataFile2D(title = "survey", periods = 1 ./ f, frequencies = f,
                site_names = ["TK" * lpad(i, 2, '0') for i in 1:ns], receivers = y, x_positions = zeros(ns), z_positions = zeros(ns),
                z_xy = zeros(ComplexF64, nf, ns), z_xy_error = fill(0.05, nf, ns),
                z_yx = zeros(ComplexF64, nf, ns), z_yx_error = fill(0.05, nf, ns),
                z_xx = complex.(nan), z_xx_error = copy(nan), z_yy = complex.(nan), z_yy_error = copy(nan),
                rho_xy = copy(nan), phase_xy = copy(nan), rho_yx = copy(nan), phase_yx = copy(nan),
                latitudes = fill(62.25, ns), longitudes = 25.75 .+ y ./ 51_900, origin = [62.25, 25.75]))
            WriteCov2D(path("cov.ctrl"), Cov2D(length(mesh.z_cell_sizes) - mesh.n_air_cells, length(mesh.y_cell_sizes)))
            fwd = ReadFwdCtrl2D(path("fwd.ctrl"))
            WriteFwdCtrl2D(path("fwd.ctrl"), FwdCtrl2D(fwd.mode, fwd.air_layers, fwd.air_thickness, fwd.air_growth,
                                                        fwd.air_resistivity, true, fwd.dipole_length, fwd.strike))
            pred = ForwardSolve2D(path("model.rho"), path("data.dat"), path("fwd.ctrl"))
            @test pred == path("data.pred")
            p = load_data2d(pred)
            @test all(isfinite, p.z_xy) && all(isfinite, p.z_yx)
            @test p.z_xy_error ≈ 0.05 .* abs.(p.z_xy) rtol = 1e-6                 # template errors are fractions
            lines = readlines(path("data.frechet"))
            dims = parse.(Int, split(lines[findfirst(startswith(">"), lines)])[2:end])
            model = ReadModel2D(path("model.rho"))
            @test dims[2] == length(model.resistivity) && dims[3:4] == collect(size(model.resistivity))
            @test dims[1] == 4 * length(p.z_xy) && count(l -> !startswith(l, "#") && !startswith(l, ">"), lines) == dims[1]

            # observed data from a model with a conductor, inverted from a masked halfspace
            truth = copy(model.resistivity)
            truth[2:4, 5:8] .= 10.0
            WriteModel2D(path("model.true"), model.y_cell_sizes, model.z_cell_sizes, truth)
            ForwardSolve2D(path("model.true"), path("data.dat"), path("fwd.ctrl"); output_path = path("obs.dat"))
            cov = ReadCov2D(path("cov.ctrl"))
            cov.mask[1, :] .= 0
            WriteCov2D(path("cov.ctrl"), cov)
            WriteModel2D(path("model.prior"), model.y_cell_sizes, model.z_cell_sizes, model.resistivity)
            write(path("inv.ctrl"), replace(read(path("inv.ctrl"), String),
                                            r"Maximum number of iterations\s*:\s*\d+" => "Maximum number of iterations : 1"))
            run = Invert2D(path("model.rho"), path("obs.dat"), path("fwd.ctrl"), path("inv.ctrl"), path("cov.ctrl"),
                           path("model.prior"))
            @test dirname(run.run_dir) == dir && startswith(basename(run.run_dir), "run_")
            @test all(f -> isfile(joinpath(run.run_dir, f)), ("model.rho", "data.pred", "History.csv", "Summary.txt"))
            @test isfile(joinpath(run.run_dir, "inputs", "inv.ctrl"))
            @test length(run.history) == 2 && run.rms < run.history[1].rms
            final = ReadModel2D(joinpath(run.run_dir, "model.rho")).resistivity
            @test final[1, :] ≈ model.resistivity[1, :] rtol = 1e-7                # masked layer stays fixed
            @test !(final[2:end, :] ≈ model.resistivity[2:end, :])
            plots = PlotInversion2D(run; true_model_path = path("model.true"), maximum_depth_km = 5.0)
            @test all(isfile, plots) && all(p -> startswith(p, joinpath(run.run_dir, "plots")), plots)
            @test isfile(PlotData2D(path("obs.dat"); predicted_path = pred, output_path = path("DataFit.png")))

            # vfsa from five files: the mask alone, no covariance and no prior
            WriteVFSACtrl2D(path("vfsa.ctrl"), VFSACtrl2D(target_rms = 1.0, max_iter = 1, chains = 2, control_points = 8,
                                                          log_bounds = (0.0, 4.0)))
            mask = ones(Int, size(model.resistivity))
            mask[1, :] .= 0
            WriteMask2D(path("mask.ctrl"), mask)
            @test_throws ErrorException Invert2D(path("model.rho"), path("obs.dat"), path("fwd.ctrl"), path("vfsa.ctrl"),
                                                 path("cov.ctrl"), path("model.prior"))
            vrun = VFSA2D(path("model.rho"), path("obs.dat"), path("fwd.ctrl"), path("vfsa.ctrl"), path("mask.ctrl"))
            @test vrun.algorithm == :vfsa && isfinite(vrun.rms) && vrun.history === nothing
            @test isfile(joinpath(vrun.run_dir, "inputs", "mask.ctrl")) && !isfile(joinpath(vrun.run_dir, "inputs", "cov.ctrl"))
            vfinal = ReadModel2D(joinpath(vrun.run_dir, "model.rho")).resistivity
            @test vfinal[1, :] ≈ model.resistivity[1, :] rtol = 1e-5               # masked layer stays fixed
            @test vrun.converged == (vrun.rms <= 1.0) && vrun.reason in (:target_rms, :max_iter)
            @test all(f -> isfile(joinpath(vrun.run_dir, f)), ("model.rho", "data.pred", "Summary.txt"))
            @test isdir(joinpath(vrun.run_dir, "vfsa"))
        end
    end

    @testset "topography and water" begin
        (; topo, data, model, lake) = _topo_case()
        t = Topography2D(model, data, topo; water = [lake])
        fwd = FwdCtrl2D(mode = :TETM, air_layers = 4, air_thickness = 20_000.0, air_growth = 3.0, air_resistivity = 1e9)
        m, r = Mesh2DFromInputs(t.model, t.data, fwd)
        mktempdir() do dir
            # six-file run: water and air fixed, bad masks refused
            observed = data_from_response2d(run_mt2d_forward(m, r); site_names = t.data.site_names,
                                            z_positions = t.data.z_positions, latitudes = data.latitudes,
                                            longitudes = data.longitudes, origin = data.origin)
            paths = (start = WriteModel2D(joinpath(dir, "model.start"), model.y_cell_sizes, model.z_cell_sizes, t.model.resistivity),
                     data = write_data2d(joinpath(dir, "data.dat"), observed),
                     fwd = WriteFwdCtrl2D(joinpath(dir, "fwd.ctrl"), fwd),
                     inv = WriteInvCtrl2D(joinpath(dir, "inv.ctrl"), InvCtrl2D(algorithm = :gn, lambda = 1.0, target_rms = 0.0, max_iter = 1)),
                     cov = WriteCov2D(joinpath(dir, "cov.ctrl"), Cov2D(sy = fill(0.3, 15), sz = 0.3, n_smooth = 1, mask = t.mask)))
            run = Invert2D(paths.start, paths.data, paths.fwd, paths.inv, paths.cov, paths.start; run_dir = joinpath(dir, "run"))
            @test count(run.water) == 4 && !any(run.active .& (run.water .| mt2d_air_mask(run.mesh)))
            @test run.final[run.water] == run.start[run.water]
            summary = read(joinpath(dir, "run", "Summary.txt"), String)
            @test occursin("Water cells: 4", summary) && occursin("Station snapping", summary)
            @test count(>(1e15), ReadModel2D(joinpath(dir, "run", "model.rho")).resistivity) == count(==(0), t.mask)

            # vfsa from five files with mask.ctrl: air, water and fixed cells stay, the ensemble is written
            vctrl = WriteVFSACtrl2D(joinpath(dir, "vfsa.ctrl"), VFSACtrl2D(target_rms = 0.0, max_iter = 3, chains = 2,
                                                                          control_points = 10))
            mpath = WriteMask2D(joinpath(dir, "mask.ctrl"), t.mask)
            vrun = VFSA2D(paths.start, paths.data, paths.fwd, vctrl, mpath; run_dir = joinpath(dir, "vrun"))
            fixed = .!vrun.active
            @test vrun.final[fixed] ≈ vrun.start[fixed]
            @test count(vrun.water) == 4 && occursin("Water cells: 4", read(joinpath(dir, "vrun", "Summary.txt"), String))
            @test length(vrun.vfsa.chains) == 2 && all(isfinite, vrun.vfsa.ensemble.std)
            @test all(isfile, joinpath.(dir, "vrun", "vfsa", ["model.mean.rho", "Uncertainty.csv", "chain_01/best.rho", "data.best.pred"]))
            @test count(>(1e15), ReadModel2D(joinpath(dir, "vrun", "model.rho")).resistivity) == count(==(0), t.mask)
            WriteMask2D(mpath, ones(Int, 15, 20))
            @test_throws ErrorException VFSA2D(paths.start, paths.data, paths.fwd, vctrl, mpath; run_dir = joinpath(dir, "bad0"))

            WriteCov2D(paths.cov, Cov2D(sy = fill(0.3, 15), sz = 0.3, n_smooth = 1, mask = ones(Int, 15, 20)))
            @test_throws ErrorException Invert2D(paths.start, paths.data, paths.fwd, paths.inv, paths.cov, paths.start;
                                                 run_dir = joinpath(dir, "bad1"))
            wet = copy(t.mask); wet[end, 11] = 9       # under the station at y = 0
            WriteCov2D(paths.cov, Cov2D(sy = fill(0.3, 15), sz = 0.3, n_smooth = 1, mask = wet))
            @test_throws ErrorException Invert2D(paths.start, paths.data, paths.fwd, paths.inv, paths.cov, paths.start;
                                                 run_dir = joinpath(dir, "bad2"))
        end
    end
end

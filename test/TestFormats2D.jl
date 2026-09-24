# 2D file formats and the ModEM-style forward and inversion workflows
# Author: @pankajkmishra

using Test

@testset "2D formats and file workflows" begin

    #---------- control files ----------

    @testset "control files" begin
        mktempdir() do dir
            fwd = FwdCtrl2D(mode = :TE, air_layers = 7, air_thickness = 30_000.0, air_growth = 1.5,
                            air_resistivity = 1e8, write_frechet = true)
            @test ReadFwdCtrl2D(WriteFwdCtrl2D(joinpath(dir, "fwd.ctrl"), fwd)) == fwd
            inv = InvCtrl2D(algorithm = :nlcg, lambda = 3.0, target_rms = 1.2, max_iter = 40,
                            smooth_y = 2.0, nlcg_precondition = false)
            @test ReadInvCtrl2D(WriteInvCtrl2D(joinpath(dir, "inv.ctrl"), inv)) == inv
            # bounds are a VFSA setting: GN and NLCG are unbounded and reject the key
            @test !occursin("bounds", read(joinpath(dir, "inv.ctrl"), String))
            write(joinpath(dir, "bounded.ctrl"), read(joinpath(dir, "inv.ctrl"), String) * "Log10 resistivity bounds : 0 4\n")
            @test_throws ErrorException ReadInvCtrl2D(joinpath(dir, "bounded.ctrl"))
            vfsa = VFSACtrl2D(target_rms = 1.1, max_iter = 9, seed = 7, chains = 3, rbf_top = 1.5, core_expansion = 2, log_bounds = (0.5, 3.5))
            @test ReadVFSACtrl2D(WriteVFSACtrl2D(joinpath(dir, "vfsa.ctrl"), vfsa)) == vfsa
            mask = [1 1 0 1; 1 9 9 1; 1 1 1 0]
            @test ReadMask2D(WriteMask2D(joinpath(dir, "mask.ctrl"), mask)) == mask
            @test readlines(joinpath(dir, "mask.ctrl"))[1] == "4 3"
            # each file holds only its own algorithm's keys
            @test !occursin("VFSA", read(joinpath(dir, "inv.ctrl"), String)) && !occursin("GN damping", read(joinpath(dir, "inv.ctrl"), String))
            @test !any(k -> occursin(k, read(joinpath(dir, "vfsa.ctrl"), String)), ("Smallness", "Algorithm", "lambda"))
            cov = Cov2D(sy = [0.1, 0.2, 0.3], sz = 0.4, n_smooth = 2, exceptions = [(2, 1, 0.0)], mask = [1 1 0 1; 1 2 2 1; 1 1 1 1])
            back = ReadCov2D(WriteCov2D(joinpath(dir, "cov.ctrl"), cov))
            @test back.sy == cov.sy && back.sz == cov.sz && back.n_smooth == 2
            @test back.exceptions == cov.exceptions && back.mask == cov.mask
            @test !any(f -> occursin('#', read(joinpath(dir, f), String)), ("fwd.ctrl", "inv.ctrl", "vfsa.ctrl", "cov.ctrl", "mask.ctrl"))

            # only the required keys, the rest take defaults
            path = joinpath(dir, "min.ctrl")
            write(path, "Algorithm : GN\nInitial damping factor lambda : 1\nExit search when rms is less than : 1.05\n" *
                        "Maximum number of iterations : 5   # trailing comment\n")
            c = ReadInvCtrl2D(path)
            @test c.algorithm == :gn && c.max_iter == 5 && c.target_rms == 1.05 && c.mode == :TETM

            write(path, "Algorithm : VFSA\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\n" *
                        "Maximum number of iterations : 5\n")
            @test_throws ErrorException ReadInvCtrl2D(path)                       # vfsa has its own control
            write(path, "Exit search when rms is less than : 1\nMaximum number of iterations : 5\nVFSA chains : 2\n")
            @test_throws ErrorException ReadVFSACtrl2D(path)                      # old prefixed keys
            write(path, "1 2\n1\n")
            @test_throws ErrorException ReadMask2D(path)                          # one value short
            write(path, "Algorithm : GN\nMaximum number of iterations : 5\n")
            @test_throws ErrorException ReadInvCtrl2D(path)                       # missing required keys
            write(path, "Mode : TE\nAir layers : 5\nAir thickness (m) : 1e4\nAir growth factor : 2\n" *
                        "Air resistivity (ohm m) : 1e9\nAir colour : blue\n")
            @test_throws ErrorException ReadFwdCtrl2D(path)                       # unknown key
            write(path, "Mode : TE\nAir layers : 5\nAir thickness (m) : 1e4\nAir growth factor : 2\n")
            @test_throws ErrorException ReadFwdCtrl2D(path)                       # no default air resistivity
            write(path, "Mode : TE\nMode : TM\n")
            @test_throws ErrorException ReadFwdCtrl2D(path)                       # duplicate key
        end

        # the shipped controls read cleanly, one per algorithm
        ctrl_dir = joinpath(dirname(@__DIR__), "examples", "ctrl", "2D")
        @test ReadFwdCtrl2D(joinpath(ctrl_dir, "FwdCtrl")).mode == :TETM
        @test [ReadInvCtrl2D(joinpath(ctrl_dir, "InvCtrl.$a")).algorithm for a in ("GN", "NLCG")] == [:gn, :nlcg]
        @test ReadVFSACtrl2D(joinpath(ctrl_dir, "InvCtrl.VFSA")).chains >= 1
    end

    #---------- model files and mesh assembly ----------

    @testset "model files and mesh" begin
        mesh = BuildMesh2D(frequencies = [1.0, 10.0], receiver_positions = [-500.0, 500.0], y_core_range = (-1000.0, 1000.0),
                           y_core_cell = 250.0, y_padding = 3000.0, pad_factor = 1.5, air_top = -10_000.0, air_cells = 5,
                           max_core_layers = 6)
        ρ = build_mt2d_halfspace_model(mesh; background_resistivity = 100.0)
        ρ[mesh.n_air_cells+2, 3:5] .= 7.5
        mktempdir() do dir
            path = WriteModel2D(joinpath(dir, "model.rho"), mesh, ρ)
            lines = readlines(path)
            @test startswith(lines[1], "# 2D MT model")
            @test split(lines[2]) == [string(length(mesh.y_cell_sizes)), string(length(mesh.z_cell_sizes) - mesh.n_air_cells), "LOGE"]
            # fixed-width columns: every full row of sizes or values has the same length
            @test length(unique(length.(filter(l -> length(split(l)) == 10, lines[3:end])))) <= 2
            m = ReadModel2D(path)
            @test m.n_air_cells == 0
            @test m.resistivity ≈ ρ[mesh.n_air_cells+1:end, :] rtol = 1e-5
            @test m.origin[2] ≈ mesh.y_nodes[1]
            # older layout with air rows reads as the same earth model
            legacy = ReadModel2D(MTGeophysics._write_model2d_legacy(joinpath(dir, "legacy.rho"), mesh, ρ))
            @test legacy.resistivity ≈ m.resistivity rtol = 1e-5
            @test all(abs.(legacy.z_cell_sizes .- m.z_cell_sizes) .<= 5e-4)
            shifted = MT2DMesh(y_nodes = mesh.y_nodes .+ 10, z_nodes = mesh.z_nodes, y_cell_sizes = mesh.y_cell_sizes,
                               z_cell_sizes = mesh.z_cell_sizes, receiver_positions = Float64[], frequencies = [1.0],
                               n_air_cells = mesh.n_air_cells)
            @test_throws ArgumentError WriteModel2D(joinpath(dir, "off.rho"), shifted, ρ)
        end

        air = mt2d_air_layers(6, 20_000.0, 2.0)
        @test sum(air) ≈ 20_000 && all(air[1:end-1] ./ air[2:end] .≈ 2.0)
        @test mt2d_air_layers(4, 100.0, 1.0) ≈ fill(25.0, 4)

        model = ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = fill(100.0, 10), z_cell_sizes = fill(50.0, 4),
                            resistivity = fill(30.0, 4, 10), n_air_cells = 0, origin = [0.0, -500.0, 0.0], rotation = 0.0,
                            format = "LOGE")
        survey = (receivers = [-250.0, 250.0], frequencies = [1.0, 10.0], z_positions = zeros(2), site_names = ["A", "B"])
        fwd = FwdCtrl2D(mode = :TETM, air_layers = 4, air_thickness = 8000.0, air_growth = 2.0, air_resistivity = 1e7)
        m, r = Mesh2DFromInputs(model, survey, fwd)
        @test m.n_air_cells == 4 && m.air_resistivity == 1e7 && size(r) == (8, 10)
        @test all(r[1:4, :] .== 1e7) && all(r[5:end, :] .== 30.0)
        @test m.z_nodes[5] ≈ 0 && m.z_nodes[1] ≈ -8000 && m.y_nodes[1] ≈ -500 && m.y_nodes[end] ≈ 500
        @test_throws ArgumentError Mesh2DFromInputs(model, merge(survey, (receivers = [-250.0, 900.0],)), fwd)
        # stations snap to the ground, warned beyond half a surface cell (25 m here)
        @test_logs Mesh2DFromInputs(model, merge(survey, (z_positions = [0.0, 12.0],)), fwd)
        @test_logs (:warn, r"moved to the ground") Mesh2DFromInputs(model, merge(survey, (z_positions = [0.0, 40.0],)), fwd)

        # geometric air rounds in the older model layout; the surface row must still be found
        mktempdir() do dir
            reloaded = MTGeophysics._load_model2d_legacy(MTGeophysics._write_model2d_legacy(joinpath(dir, "air.rho"), m, r))
            back = MTGeophysics._mesh_from_legacy_model2d(reloaded; frequencies = [1.0, 10.0], receiver_positions = [-250.0, 250.0])
            @test all(isfinite, run_mt2d_forward(back, reloaded.resistivity).z_xy)
        end

        layers = mt2d_geometric_layers([0.1, 1000.0]; background_resistivity = 100.0, first_layer_div = 5.0,
                                       vertical_factor = 1.1, depth_mult = 4.0)
        @test layers[1] ≈ mt2d_skin_depth(100.0, 1000.0) / 5
        @test all(layers[2:end] ./ layers[1:end-1] .≈ 1.1)
        @test sum(layers[1:end-1]) < 4 * mt2d_skin_depth(100.0, 0.1) <= sum(layers)
    end

    #---------- data format and sign convention ----------

    @testset "data format and TE/TM sign" begin
        mesh = BuildMesh2D(frequencies = [0.1, 1.0, 10.0], receiver_positions = [-400.0, 0.0, 400.0],
                           y_core_range = (-1000.0, 1000.0), y_core_cell = 200.0, y_padding = 20_000.0, pad_factor = 1.3,
                           air_top = -40_000.0, air_cells = 10, max_core_layers = 30)
        ρ = build_mt2d_halfspace_model(mesh; background_resistivity = 100.0)
        response = run_mt2d_forward(mesh, ρ)
        # exp(+iωt): Zxy in the first quadrant, Zyx in the third, and Zyx = -Zxy over a halfspace
        @test all(z -> real(z) > 0 && imag(z) > 0, response.z_xy)
        @test all(z -> real(z) < 0 && imag(z) < 0, response.z_yx)
        @test response.z_yx ≈ -response.z_xy rtol = 0.02
        @test all(isapprox.(response.rho_xy, 100.0; rtol = 0.05))

        data = data_from_response2d(response; site_names = ["A1", "A2", "A3"], latitudes = [62.2, 62.25, 62.3],
                                    longitudes = [25.7, 25.75, 25.8], origin = [62.25, 25.75])
        mktempdir() do dir
            path = write_data2d(joinpath(dir, "data.dat"), data)
            text = read(path, String)
            @test occursin("> Full_Impedance", text) && occursin("> exp(+i\\omega t)", text)
            @test occursin("> [mV/km]/[nT]", text) && !occursin("ZXX", text) && !occursin("TX", text)
            @test startswith(text, "# 2D MT data")
            # fixed-width columns, ModEM style
            @test length(unique(length.(filter(l -> !startswith(l, r"[#>]"), split(chomp(text), '\n'))))) == 1
            back = load_data2d(path)
            @test back.z_xy ≈ data.z_xy rtol = 1e-6
            @test back.z_yx ≈ data.z_yx rtol = 1e-6
            @test back.z_xy_error ≈ data.z_xy_error rtol = 1e-6
            @test back.site_names == data.site_names && back.receivers == data.receivers
            @test back.latitudes == data.latitudes && back.longitudes == data.longitudes && back.origin == [62.25, 25.75]
            @test issorted(back.frequencies)
            # the 3D reader sees the same impedances, so one file serves 1D, 2D and 3D
            d3 = load_data_modem(path)
            @test d3.Z[:, 2, :] ≈ data.z_xy[end:-1:1, :] rtol = 1e-6
            @test d3.Z[:, 3, :] ≈ data.z_yx[end:-1:1, :] rtol = 1e-6
        end
    end

    #---------- forward and inversion from files ----------

    @testset "file workflows" begin
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
                site_names = ["SITE$i" for i in 1:ns], receivers = y, x_positions = zeros(ns), z_positions = zeros(ns),
                z_xy = zeros(ComplexF64, nf, ns), z_xy_error = fill(0.05, nf, ns),
                z_yx = zeros(ComplexF64, nf, ns), z_yx_error = fill(0.05, nf, ns),
                z_xx = complex.(nan), z_xx_error = copy(nan), z_yy = complex.(nan), z_yy_error = copy(nan),
                rho_xy = copy(nan), phase_xy = copy(nan), rho_yx = copy(nan), phase_yx = copy(nan),
                latitudes = fill(62.25, ns), longitudes = 25.75 .+ y ./ 51_900, origin = [62.25, 25.75]))
            WriteCov2D(path("cov.ctrl"), Cov2D(length(mesh.z_cell_sizes) - mesh.n_air_cells, length(mesh.y_cell_sizes)))
            fwd = ReadFwdCtrl2D(path("fwd.ctrl"))
            WriteFwdCtrl2D(path("fwd.ctrl"), FwdCtrl2D(fwd.mode, fwd.air_layers, fwd.air_thickness, fwd.air_growth,
                                                        fwd.air_resistivity, true, fwd.dipole_length))
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
end

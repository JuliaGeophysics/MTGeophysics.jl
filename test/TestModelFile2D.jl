# 2D model files
# Author: @pankajkmishra
# Ensures the ModEM-layout model file round-trips with fixed-width columns and its origin, that the older layout with
# air rows reads as the same earth model and reproduces the forward response, and that misplaced grids are refused

using Test

@testset "2D model files" begin
    @testset "ModEM layout" begin
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
    end

    @testset "older layout round trip" begin
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
        end
    end
end

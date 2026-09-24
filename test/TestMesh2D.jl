# 2D mesh building
# Author: @pankajkmishra
# Ensures the skin-depth, geometric and air layers have the documented sizes, and that a model, a survey and fwd.ctrl
# assemble into the solver mesh with the air on top, stations inside it and snapping beyond half a cell warned about

using Test

@testset "2D mesh building" begin
    @testset "skin-depth layers" begin
        f = [0.1, 1.0, 10.0]
        δ = mt2d_skin_depth(100.0, 0.1)
        @test δ ≈ 503.29*sqrt(100/0.1) rtol=1e-4
        layers = mt2d_skin_depth_layers(f; background_resistivity=100.0, z_core_cell=250.0)
        core = findall(==(250.0), layers)
        @test core == 1:length(core)
        @test sum(layers[core]) >= δ && sum(layers[core]) < δ + 250
        @test sum(layers) >= 4δ
        @test all(diff(layers[length(core):end]) .> 0)
        mesh = BuildMesh2D(frequencies=f)
        @test mesh.z_cell_sizes[mesh.n_air_cells+1:end] ≈ mt2d_skin_depth_layers(f)
    end

    @testset "geometric and air layers" begin
        layers = mt2d_geometric_layers([0.1, 1000.0]; background_resistivity = 100.0, first_layer_div = 5.0,
                                       vertical_factor = 1.1, depth_mult = 4.0)
        @test layers[1] ≈ mt2d_skin_depth(100.0, 1000.0) / 5
        @test all(layers[2:end] ./ layers[1:end-1] .≈ 1.1)
        @test sum(layers[1:end-1]) < 4 * mt2d_skin_depth(100.0, 0.1) <= sum(layers)
        air = mt2d_air_layers(6, 20_000.0, 2.0)
        @test sum(air) ≈ 20_000 && all(air[1:end-1] ./ air[2:end] .≈ 2.0)
        @test mt2d_air_layers(4, 100.0, 1.0) ≈ fill(25.0, 4)
    end

    @testset "mesh from a model, a survey and fwd.ctrl" begin
        model = ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = fill(100.0, 10), z_cell_sizes = fill(50.0, 4),
                            resistivity = fill(30.0, 4, 10), n_air_cells = 0, origin = [0.0, -500.0, 0.0], rotation = 0.0,
                            format = "LOGE")
        survey = (receivers = [-250.0, 250.0], frequencies = [1.0, 10.0], z_positions = zeros(2), site_names = ["TK01", "TK02"])
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
    end
end

# 2D data files and sign convention
# Author: @pankajkmishra
# Ensures exp(+iωt) puts ZXY in the first quadrant and ZYX in the third, and that the ModEM Full_Impedance file
# round-trips in [mV/km]/[nT] with fixed-width columns, lat/lon and origin, and reads the same through the 3D reader

using Test

@testset "2D data files" begin
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

        data = data_from_response2d(response; site_names = ["TK01", "TK02", "TK03"], latitudes = [62.2, 62.25, 62.3],
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
end

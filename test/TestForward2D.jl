# 2D forward solver
# Author: @pankajkmishra
# Ensures the TE/TM solver returns finite, positive and physically consistent responses on a half-space profile

using Test

@testset "2D forward, half-space" begin
    mesh = BuildMesh2D(
        frequencies = collect(10 .^ range(-1, 1, length = 4)),
        y_core_range = (-2000.0, 2000.0),
        y_core_cell = 500.0,
        y_padding = 2500.0,
        pad_factor = 1.20,
        air_top = -3000.0,
        air_cells = 4,
        ground_layers = [200.0, 200.0, 400.0, 400.0, 800.0, 800.0],
        receiver_stride = 2,
    )

    resistivity = MTGeophysics.build_mt2d_halfspace_model(mesh; background_resistivity = 100.0)
    response = run_mt2d_forward(mesh, resistivity)
    station_index = cld(length(mesh.receiver_positions), 2)

    @test size(response.rho_xy) == (length(mesh.frequencies), length(mesh.receiver_positions))
    @test size(response.rho_yx) == (length(mesh.frequencies), length(mesh.receiver_positions))
    @test all(isfinite, response.rho_xy)
    @test all(isfinite, response.rho_yx)
    @test all(isfinite, response.phase_xy)
    @test all(isfinite, response.phase_yx)
    @test minimum(response.rho_xy) > 0
    @test minimum(response.rho_yx) > 0
    @test maximum(abs.((response.rho_xy[:, station_index] .- response.rho_yx[:, station_index]) ./ response.rho_xy[:, station_index])) < 0.40
    @test maximum(abs.(response.phase_xy[:, station_index] .- MTGeophysics._phase_fold_to_0_90(response.phase_yx[:, station_index]))) < 35.0
end

@testset "2D forward, TM receiver Ey at a lateral contact" begin
    Jy, dz = 1e-4, 10.0
    for (y_nodes, ys) in [([-1000.0, -500.0, 0.0, 500.0, 1000.0], [-375.0, -125.0, 0.0, 125.0, 375.0, 500.0]),
                          ([-1000.0, -400.0, 0.0, 150.0, 900.0], [-300.0, -100.0, 0.0, 50.0, 100.0, 150.0])]
        ρ = [10.0, 10.0, 100.0, 100.0]
        H = hcat(fill(1.0 + 0im, 5), fill(1.0 + Jy * dz + 0im, 5))
        surface = MTGeophysics._mt2d_surface_fields_tm(0.0, y_nodes, dz, 1 ./ ρ, H)
        plan = (stations = [(group = 1, window = nothing) for _ in ys],)
        locations = hcat(ys, zeros(length(ys)))
        Z = MTGeophysics._mt2d_station_impedance(:TM, plan, [surface], locations, y_nodes)
        Erec, Hrec = MTGeophysics.compute_fields_at_receivers_tm(0.0, locations, y_nodes, dz, 1 ./ ρ, H)
        exact = [ρ[searchsortedlast(y_nodes, y)] * Jy for y in ys]
        @test Z ≈ exact rtol = 1e-10
        @test Erec ./ Hrec ≈ exact rtol = 1e-10
    end
end

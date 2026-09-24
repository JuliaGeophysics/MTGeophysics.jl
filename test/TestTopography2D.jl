# 2D topography
# Author: @pankajkmishra
# Ensures stations sit on the surface rows of a staircase ground, uniform topography equals a flat surface one row down,
# models with topographic air round-trip and snap stations, Topography2D cuts in air, water and station columns (up or
# down), and the trapezoidal hill of Wannamaker, Stodt & Rijo (1986) is reproduced

using Test, LinearAlgebra, Random

@testset "2D topography" begin
    mesh = _hill_mesh()
    air = mt2d_air_mask(mesh)
    rng = MersenneTwister(4)
    ρ = build_mt2d_halfspace_model(mesh)
    ρ[3:end, :] .*= exp.(0.3randn(rng, size(ρ, 1) - 2, size(ρ, 2)))
    fields = (:rho_xy, :phase_xy, :z_xy, :rho_yx, :phase_yx, :z_yx)

    @testset "mesh helpers" begin
        @test mt2d_receiver_depths(mesh) ≈ [100.0, 50.0, 0.0, 0.0, 50.0, 100.0]
        @test count(air) == 2 * size(ρ, 2) + sum(mesh.topo_air)
        t, _, _, locations = MTGeophysics._assemble_mt2d_system(mesh, ρ)
        plan = MTGeophysics._mt2d_station_plan(mesh, t, locations)
        @test [g.z_index for g in plan.groups] == [3, 4, 5]
        @test all(p -> p.window !== nothing, plan.stations)          # every station is next to a step here
        flat = MTGeophysics._mt2d_station_plan(_remesh(mesh; topo_air = Int[]), t, hcat(locations[:, 1], fill(1000.0, 6)))
        @test all(p -> p.window === nothing, flat.stations)
    end

    @testset "uniform topography equals a flat surface one row down" begin
        shifted = _remesh(mesh; topo_air = fill(1, length(mesh.y_cell_sizes)))
        flat = _remesh(mesh; n_air_cells = 3, topo_air = Int[])
        a, b = run_mt2d_forward(shifted, ρ), run_mt2d_forward(flat, ρ)
        @test all(k -> getproperty(a, k) ≈ getproperty(b, k), fields)
        w = (z_xy = randn(rng, ComplexF64, 2, 6), z_yx = randn(rng, ComplexF64, 2, 6))
        @test ApplyFrechetTranspose2D(shifted, ρ, w) ≈ ApplyFrechetTranspose2D(flat, ρ, w)
    end

    @testset "model files and station snapping" begin
        mktempdir() do dir
            path = WriteModel2D(joinpath(dir, "topo.rho"), mesh, ρ)
            model = ReadModel2D(path)
            @test count(>(1e15), model.resistivity) == sum(mesh.topo_air)
            fwd = FwdCtrl2D(mode = :TETM, air_layers = 2, air_thickness = 1000.0, air_growth = 1.0, air_resistivity = 1e9)
            data = (receivers = mesh.receiver_positions, frequencies = mesh.frequencies,
                    z_positions = mt2d_receiver_depths(mesh), site_names = ["TK" * lpad(i, 2, '0') for i in 1:6])
            m2, r2 = @test_logs Mesh2DFromInputs(model, data, fwd)
            @test m2.topo_air == mesh.topo_air && all(r2[mt2d_air_mask(m2)] .== 1e9)
            @test run_mt2d_forward(m2, r2).z_xy ≈ run_mt2d_forward(mesh, ρ).z_xy rtol = 1e-4
            @test all(o -> o.offset == 0, mt2d_station_offsets(m2, data))
            @test_logs (:warn, r"moved to the ground") Mesh2DFromInputs(model, merge(data, (z_positions = data.z_positions .+ 40,)), fwd)
            stray = deepcopy(model)
            stray.resistivity[end, 1] = 1e17
            @test_throws ArgumentError Mesh2DFromInputs(stray, data, fwd)
        end
    end

    # published behaviour of the trapezoidal hill of Wannamaker, Stodt & Rijo (1986): 450 m high, 2 km base,
    # 450 m flat top, 100 ohm m, 2 Hz; TM low on the hilltop, high at its foot, TE mildly raised on top, both
    # back to 100 ohm m away from the hill. 50 m staircase cells, stations at cell centres, mirror symmetric
    @testset "trapezoidal hill" begin
        h = 50.0
        elevation(y) = 450 * clamp((1000 - abs(y)) / 775, 0, 1)
        core = collect(-3500.0:h:3500.0)
        pad(n0) = (nodes = Float64[]; Δ = h; y = n0;
                   while abs(y) < 60_000; Δ *= 1.2; y += sign(n0) * Δ; push!(nodes, y); end; nodes)
        y_nodes = vcat(reverse(pad(core[1])), core, pad(core[end]))
        dz = fill(h, ceil(Int, 950 / h))
        while sum(dz) < 40_000
            push!(dz, dz[end] * 1.15)
        end
        air = mt2d_air_layers(12, 60_000.0, 2.0)
        zc = cumsum(dz) .- dz ./ 2
        yc = (y_nodes[1:end-1] .+ y_nodes[2:end]) ./ 2
        hill = MT2DMesh(y_nodes = y_nodes, z_nodes = vcat(0.0, cumsum(vcat(air, dz))) .- sum(air), y_cell_sizes = diff(y_nodes),
                        z_cell_sizes = vcat(air, dz), receiver_positions = collect(-2975.0:50.0:2975.0), frequencies = [2.0],
                        n_air_cells = length(air), topo_air = [count(<(450 - elevation(y)), zc) for y in yc])
        r = run_mt2d_forward(hill, fill(100.0, length(hill.z_cell_sizes), length(yc)))
        y = hill.receiver_positions
        top, far = argmin(abs.(y)), [argmin(y), argmax(y)]
        @test r.rho_yx[1, top] < 40 && 105 < r.rho_xy[1, top] < 125
        @test 850 < abs(y[argmax(r.rho_yx[1, :])]) < 1150 && maximum(r.rho_yx) > 140
        @test all(abs.(r.rho_yx[1, far] ./ 100 .- 1) .< 0.05) && all(abs.(r.rho_xy[1, far] ./ 100 .- 1) .< 0.05)
        @test maximum(abs.(r.rho_yx[1, :] ./ reverse(r.rho_yx[1, :]) .- 1)) < 0.02
        @test maximum(abs.(r.rho_xy[1, :] ./ reverse(r.rho_xy[1, :]) .- 1)) < 0.02
    end

    @testset "topography builder" begin
        (; stations, elevation, ty, topo, data, model, lake) = _topo_case()
        mktempdir() do dir
            topopath = WriteTopo2D(joinpath(dir, "topo.dat"), topo)
            back = ReadTopo2D(topopath)
            @test maximum(abs.(back.elevations .- topo.elevations)) <= 0.005
            @test maximum(abs.(back.latitudes .- topo.latitudes)) <= 1e-6
            y, _ = mt2d_profile_topography(back, data)
            @test maximum(abs.(y .- ty)) < 0.2

            t = Topography2D(model, data, back; water = [lake])
            @test t.datum ≈ 160 atol = 0.01
            @test maximum(abs.(t.data.z_positions .- (160 .- elevation.(stations)))) < 0.05
            @test t.mask[1, 10] == 1 && all(t.mask[1:3, 20] .== 0) && t.mask[4, 20] == 1
            @test count(==(9), t.mask) == 4 && all(t.model.resistivity[t.mask .== 9] .== 100.0)
            @test all(t.model.resistivity[t.mask .== 0] .== 1e17)
            @test_throws ArgumentError Topography2D(model, data, back; water = [(y_range = (-500.0, -300.0), level = 200.0)])

            fwd = FwdCtrl2D(mode = :TETM, air_layers = 4, air_thickness = 20_000.0, air_growth = 3.0, air_resistivity = 1e9)
            m, r = @test_logs Mesh2DFromInputs(t.model, t.data, fwd)
            @test all(o -> abs(o.offset) <= o.tolerance, mt2d_station_offsets(m, t.data))

            # a station in a dip narrower than its column lowers that column's ground, not only raises it
            dip(y) = 100 - 60exp(-((y + 90) / 15)^2)
            fine = collect(-1500.0:5.0:1500.0)
            sy = [-400.0, -90.0, 300.0]
            ddata = (; receivers = sy, frequencies = [1.0, 10.0], z_positions = zeros(3), site_names = ["TK01", "TK02", "TK03"],
                       _wgs(sy)..., origin = [62.25, 25.75])
            dt = Topography2D(model, ddata, Topo2D(; _wgs(fine)..., elevations = dip.(fine)))
            @test all(dt.mask[1:3, 10] .== 0) && dt.mask[4, 10] == 1 && all(dt.mask[1, [9, 11]] .== 1)
            dm, _ = @test_logs Mesh2DFromInputs(dt.model, dt.data, fwd)
            @test all(o -> abs(o.offset) <= o.tolerance, mt2d_station_offsets(dm, dt.data))
        end
    end
end

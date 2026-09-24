# 2D topography: forward, Fréchet derivatives, model files, station snapping, regularization,
# the topography builder and the six-file inversion with air and water
# Author: @pankajkmishra

using Test, LinearAlgebra, Random, SparseArrays

# a small hill: no topographic air in the middle, one cell on the flanks, two outside;
# the receivers sit on three surface rows, which share node rows
function _hill_mesh()
    base = BuildMesh2D(frequencies = [0.3, 3.0], y_core_range = (-750.0, 750.0), y_core_cell = 250.0,
                       y_padding = 600.0, air_cells = 2, air_top = -1000.0, ground_layers = [50.0, 50.0, 100.0, 300.0, 900.0],
                       receiver_positions = [-625.0, -375.0, -125.0, 125.0, 400.0, 600.0])
    yc = (base.y_nodes[1:end-1] .+ base.y_nodes[2:end]) ./ 2
    topo = [abs(y) < 300 ? 0 : abs(y) < 600 ? 1 : 2 for y in yc]
    _remesh(base; topo_air = topo)
end

_remesh(m; kw...) = MT2DMesh(; (k => getfield(m, k) for k in fieldnames(MT2DMesh) if !haskey(kw, k))..., kw...)

# E-W stations near Jyväskylä at profile positions y
function _wgs(y)
    trans = MTGeophysics.Proj.Transformation(MTGeophysics._local_tm_proj_string(62.25, 25.75), "EPSG:4326"; always_xy = true)
    p = [trans(500_000.0 + yi, 0.0) for yi in y]
    (latitudes = [q[2] for q in p], longitudes = [q[1] for q in p])
end

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

    # TM surface fields next to air involve tiny currents times ρ_air, so the finite-difference
    # check uses 1e6 ohm m air to keep its noise below the tolerance
    @testset "Fréchet derivatives with topography" begin
        mesh = _remesh(mesh; air_resistivity = 1e6)
        response = run_mt2d_forward(mesh, ρ)
        @test all(k -> all(isfinite, getproperty(response, k)), fields)
        direction = randn(rng, size(ρ))
        h = 1e-4
        for mode in (:TE, :TM, :TETM)
            tangent = ApplyFrechet2D(mesh, ρ, direction; mode, parameterization = :log_resistivity)
            rp = run_mt2d_forward(mesh, ρ .* exp.(h * direction); mode)
            rm = run_mt2d_forward(mesh, ρ .* exp.(-h * direction); mode)
            for key in fields
                fd = (getproperty(rp, key) - getproperty(rm, key)) / (2h)
                @test getproperty(tangent, key) ≈ fd rtol = 2e-5 atol = 1e-10
            end
            weights = NamedTuple{fields}((randn(rng, 2, 6), randn(rng, 2, 6), randn(rng, ComplexF64, 2, 6),
                                          randn(rng, 2, 6), randn(rng, 2, 6), randn(rng, ComplexF64, 2, 6)))
            gradient = ApplyFrechetTranspose2D(mesh, ρ, weights; mode, parameterization = :log_resistivity)
            @test dot(gradient, direction) ≈ sum(real(dot(getproperty(weights, k), getproperty(tangent, k))) for k in fields) rtol = 1e-8
            @test all(iszero, gradient[air])
        end
        hilltop_air = CartesianIndex(3, findfirst(>(0), mesh.topo_air))
        G = FrechetDerivative2D(mesh, ρ; active_cells = [hilltop_air, CartesianIndex(6, 4)])
        @test all(iszero, G.z_xy[:, 1]) && all(iszero, G.z_yx[:, 1]) && any(!iszero, G.z_xy[:, 2])
        @test FrechetDerivative2D(mesh, ρ; active_cells = findall(.!air)[1:2]).cells == findall(.!air)[1:2]
    end

    @testset "model files and station snapping" begin
        mktempdir() do dir
            path = WriteModel2D(joinpath(dir, "topo.rho"), mesh, ρ)
            model = ReadModel2D(path)
            @test count(>(1e15), model.resistivity) == sum(mesh.topo_air)
            fwd = FwdCtrl2D(mode = :TETM, air_layers = 2, air_thickness = 1000.0, air_growth = 1.0, air_resistivity = 1e9)
            data = (receivers = mesh.receiver_positions, frequencies = mesh.frequencies,
                    z_positions = mt2d_receiver_depths(mesh), site_names = ["S$i" for i in 1:6])
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

    @testset "regularization skips air and water" begin
        water = falses(size(ρ)); water[5:6, 1] .= true
        excluded = air .| water
        R = MTGeophysics._inv2d_regularizer(mesh, ρ, Inv2DOptions(); excluded)
        @test nnz(R[:, findall(vec(excluded))]) == 0
        @test nnz(R[:, findall(vec(.!excluded))]) > 0

        observed = data_from_response2d(run_mt2d_forward(mesh, ρ))
        start = build_mt2d_halfspace_model(mesh)
        result = Invert2D(mesh, start, observed; algorithm = GaussNewton2DConfig(),
                          options = Inv2DOptions(max_iter = 2, target_rms = 0.0, verbose = false), water_cells = water)
        @test all(c -> !air[c] && !water[c], result.active_cells)
        @test result.resistivity[water] == start[water]
        @test result.history[end].rms < result.history[1].rms
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

    @testset "topography builder and six-file inversion" begin
        stations = collect(-400.0:200.0:400.0)
        elevation(y) = 100 + 60exp(-(y / 300)^2) - 40exp(-((y + 800) / 100)^2)
        ty = collect(-1500.0:50.0:1500.0)
        topo = Topo2D(; _wgs(ty)..., elevations = elevation.(ty))
        data = (; receivers = stations, frequencies = [1.0, 10.0], z_positions = zeros(5), site_names = ["JYV$i" for i in 1:5],
                  _wgs(stations)..., origin = [62.25, 25.75])
        model = ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = fill(100.0, 20),
                            z_cell_sizes = vcat(fill(20.0, 10), [50.0, 100.0, 200.0, 400.0, 800.0]),
                            resistivity = fill(100.0, 15, 20), n_air_cells = 0, origin = [0.0, -1000.0, 0.0],
                            rotation = 0.0, format = "LOGE")
        lake = (y_range = (-950.0, -700.0), level = 105.0)

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

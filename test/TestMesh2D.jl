# 2D mesh tool and covariance masks
# Author: @pankajkmishra

using Test

@testset "2D mesh tool" begin
    trans = MTGeophysics.Proj.Transformation(MTGeophysics._local_tm_proj_string(62.25, 25.75), "EPSG:4326"; always_xy = true)
    wgs(y) = (p = [trans(500_000.0 + v, 0.0) for v in y]; ([q[2] for q in p], [q[1] for q in p]))
    y = collect(-4000.0:1000.0:4000.0)
    f = collect(10 .^ range(-1, 3, length = 9))
    nf, ns = length(f), length(y)
    lat, lon = wgs(y)
    nan = fill(NaN, nf, ns)
    Z = repeat(mt1d_impedance(f, [100.0], Float64[]), 1, ns)          # a 100 ohm m halfspace
    data = DataFile2D(title = "", periods = 1 ./ f, frequencies = f, site_names = ["S$i" for i in 1:ns], receivers = y,
                      x_positions = zeros(ns), z_positions = zeros(ns), z_xy = Z,
                      z_xy_error = 0.05 .* abs.(Z), z_yx = -Z, z_yx_error = 0.05 .* abs.(Z),
                      z_xx = complex.(nan), z_xx_error = copy(nan), z_yy = complex.(nan), z_yy_error = copy(nan),
                      rho_xy = fill(100.0, nf, ns), phase_xy = fill(45.0, nf, ns), rho_yx = fill(100.0, nf, ns),
                      phase_yx = fill(-135.0, nf, ns), latitudes = lat, longitudes = lon, origin = [62.25, 25.75])
    ty = collect(-80_000.0:200.0:80_000.0)
    tlat, tlon = wgs(ty)
    relief(v) = 150 + 60sin(2π * v / 7000) - 70exp(-((v + 30_000) / 6000)^2)
    topo = Topo2D(latitudes = tlat, longitudes = tlon, elevations = relief.(ty))

    mktempdir() do dir
        dpath = write_data2d(joinpath(dir, "data.dat"), data)
        tpath = WriteTopo2D(joinpath(dir, "topo.dat"), topo)

        flat = MakeMesh2D(dpath; out_dir = joinpath(dir, "flat"))
        @test all(isfile, [flat.paths.start_model_path, flat.paths.cov_path, flat.paths.fwd_path, flat.paths.inv_path,
                           joinpath(dir, "flat", "Mesh.png")])
        m = ReadModel2D(flat.paths.start_model_path)
        @test all(==(1), ReadCov2D(flat.paths.cov_path).mask)
        @test ReadMask2D(flat.paths.mask_path) == ReadCov2D(flat.paths.cov_path).mask
        @test ReadVFSACtrl2D(flat.paths.vfsa_path).chains >= 1
        @test all(isapprox.(m.resistivity, 100.0; rtol = 1e-4))
        @test abs(sum(m.y_cell_sizes) / 2 + m.origin[2]) < 1e-6
        @test count(==(500.0), m.y_cell_sizes) >= 16                  # half the 1 km spacing over the core
        @test isempty(flat.notes)

        lake = [(y_range = (-35_000.0, -25_000.0), level = 110.0)]
        t = MakeMesh2D(dpath; out_dir = joinpath(dir, "topo"), topo_path = tpath, water = lake, fixed_below_m = 30_000.0)
        cov = ReadCov2D(t.paths.cov_path)
        @test ReadMask2D(t.paths.mask_path) == cov.mask
        mt = ReadModel2D(t.paths.start_model_path)
        @test any(==(0), cov.mask) && any(==(9), cov.mask) && size(cov.mask) == size(mt.resistivity)
        @test all(cov.mask[mt.resistivity .> 1e15] .== 0)
        zt = vcat(0.0, cumsum(mt.z_cell_sizes)[1:end-1])
        @test all(cov.mask[zt .>= 30_000.0, :] .== 0)
        observed = load_data2d(t.paths.data_path)
        @test minimum(observed.z_positions) ≈ 0 atol = 1e-6
        @test maximum(observed.z_positions) > 50
        fwd = ReadFwdCtrl2D(t.paths.fwd_path)
        @test fwd.air_layers == 10 && fwd.dipole_length == 100
        # the tool only raises the ground in station columns (as MakeMesh3D), so a station in a dip narrower
        # than a column sits below its column's ground and is snapped up, with a warning and a note
        mesh, _ = Mesh2DFromInputs(mt, observed, fwd; warn = false)
        buried = filter(o -> abs(o.offset) > o.tolerance, mt2d_station_offsets(mesh, observed))
        @test all(o -> o.z > o.surface, buried)
        @test isempty(buried) || any(n -> occursin("snap by more than half a cell", n), t.notes)
        @test any(>(0), mt2d_topo_air(mesh))
        yc, ground = mt2d_ground(mt)
        @test length(yc) == size(mt.resistivity, 2) && maximum(ground) > 0

        mask = Mask2D(mt; fixed = [(y_range = (-1000.0, 1000.0), z_range = (0.0, 500.0))])
        yc0 = mt.origin[2] .+ cumsum(mt.y_cell_sizes) .- mt.y_cell_sizes ./ 2
        @test all(mask[1:2, findall(abs.(yc0) .< 900)] .== 0)
        @test_throws ArgumentError MakeMesh2D(dpath; out_dir = joinpath(dir, "bad"), colour = 1)

        # the mesh window, headless on CairoMakie: sliders rebuild the mesh and Save writes what batch mode writes
        ctrls = (inv = joinpath(pkgdir(MTGeophysics), "examples", "ctrl", "2D", "InvCtrl.GN"),
                 vfsa = joinpath(pkgdir(MTGeophysics), "examples", "ctrl", "2D", "InvCtrl.VFSA"))
        w = MTGeophysics._makemesh2d_window(load_data2d(dpath), ReadTopo2D(tpath), lake, MTGeophysics._MAKEMESH2D_DEFAULTS,
                                            joinpath(dir, "gui"), ctrls)
        MTGeophysics.CairoMakie.set_close_to!(w.grid.sliders[1], 0.25)
        MTGeophysics.CairoMakie.set_close_to!(w.grid.sliders[3], 8)
        @test w.params().cell_width_frac == 0.25 && w.params().n_pad == 8
        @test occursin("250 m", w.info.text[])
        w.save_button.clicks[] += 1
        b = MakeMesh2D(dpath; out_dir = joinpath(dir, "batch"), topo_path = tpath, water = lake, inv_ctrl = ctrls.inv,
                       vfsa_ctrl = ctrls.vfsa, cell_width_frac = 0.25, n_pad = 8)
        for f in ("model.start", "model.prior", "cov.ctrl", "mask.ctrl", "fwd.ctrl", "inv.ctrl", "vfsa.ctrl", "data.dat")
            @test read(joinpath(dir, "gui", f)) == read(joinpath(dir, "batch", f))
        end
    end
end

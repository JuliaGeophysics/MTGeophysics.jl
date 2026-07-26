# Topography / bathymetry extraction, air and water masks, and depth-weighted
# RBF placement tests (headless, no solver). Builds tiny synthetic WS3D models
# through the file round-trip so grid conventions match load_ws3d_model exactly.

using Random

@testset "Bathymetry and water masks (Mask3D)" begin

    nx, ny, nz = 6, 5, 4
    dx = fill(1000.0, nx); dy = fill(1000.0, ny)
    dz = [100.0, 200.0, 400.0, 800.0]        # bottoms at 100, 300, 700, 1500 m
    A = fill(2.0, nx, ny, nz)
    # water: columns i=1:2 (all j) wet through layer 2; column (3,1) wet layer 1
    A[1:2, :, 1:2] .= -0.5
    A[3, 1, 1] = -0.5

    path = joinpath(mktempdir(), "tiny.ws")
    write_ws3d_model(path, dx, dy, dz, A)
    m = load_ws3d_model(path)

    expected = falses(nx, ny, nz)
    expected[1:2, :, 1:2] .= true
    expected[3, 1, 1] = true

    #---------- threshold mask ----------
    mask_thr = water_mask_from_model(m; water_log10=0.5)
    @test mask_thr == expected
    @test count(mask_thr) == 21                     # 2x5x2 block + 1 cell

    #---------- bathymetry extraction ----------
    bathy = extract_bathymetry(m; water_log10=0.5)
    @test length(bathy.depth) == 11                 # 10 columns at 300 m + 1 at 100 m
    @test count(==(300.0), bathy.depth) == 10
    @test count(==(100.0), bathy.depth) == 1

    #---------- file round trip and mask equivalence ----------
    bpath = joinpath(mktempdir(), "bathy.dat")
    write_bathymetry(bpath, bathy; water_log10=0.5)
    bathy2 = read_bathymetry(bpath)
    @test bathy2.x ≈ bathy.x
    @test bathy2.y ≈ bathy.y
    @test bathy2.depth ≈ bathy.depth
    @test water_mask_from_bathymetry(m, bathy2) == expected

    #---------- rbf placement: exclusion and depth weighting ----------
    rng = MersenneTwister(7)
    rbfmap = build_rbf_map(m, 1:nx, 1:ny, 40, rng;
                           depth_power=1.0, exclude=expected)
    M = length(rbfmap.ci)
    @test M == 40
    @test all(t -> !expected[rbfmap.ci[t], rbfmap.cj[t], rbfmap.ck[t]], 1:M)

    # requesting more controls than non-excluded cells warns and truncates
    rbfmap_full = @test_logs (:warn, r"only") match_mode=:any build_rbf_map(
        m, 1:nx, 1:ny, nx*ny*nz, rng; exclude=expected)
    @test length(rbfmap_full.ci) == nx*ny*nz - 21

    # depth_power = 0, no exclusion: uniform placement still works
    rbfmap_u = build_rbf_map(m, 1:nx, 1:ny, 25, rng)
    @test length(rbfmap_u.ci) == 25
end

@testset "Topography and air masks (Mask3D)" begin

    nx, ny, nz = 6, 5, 4
    dx = fill(1000.0, nx); dy = fill(1000.0, ny)
    dz = [100.0, 200.0, 400.0, 800.0]        # bottoms at 100, 300, 700, 1500 m
    A = fill(2.0, nx, ny, nz)
    # air: columns i=1:2 (all j) down through layer 2; column (3,1) layer 1 only
    A[1:2, :, 1:2] .= 17.0
    A[3, 1, 1] = 17.0

    origin = [-4000.0, -3000.0, -30.0]       # nonzero z datum: cz != depth
    path = joinpath(mktempdir(), "topo.ws")
    write_ws3d_model(path, dx, dy, dz, A, origin)
    m = load_ws3d_model(path)
    @test m.origin[3] == -30.0

    expected = falses(nx, ny, nz)
    expected[1:2, :, 1:2] .= true
    expected[3, 1, 1] = true

    #---------- air is tagged NaN by the loader ----------
    mask_air = air_mask_from_model(m)
    @test mask_air == expected
    @test count(mask_air) == 21
    @test all(isnan, m.A[expected])

    #---------- topography extraction ----------
    topo = extract_topography(m)
    @test length(topo.depth) == 11                  # 10 columns to 300 m + 1 to 100 m
    @test count(==(300.0), topo.depth) == 10
    @test count(==(100.0), topo.depth) == 1
    # base_z carries origin[3] so the surface plots against m.z / m.cz
    @test all(topo.z .≈ topo.depth .+ origin[3])

    #---------- file round trip and mask equivalence ----------
    tpath = joinpath(mktempdir(), "topography.dat")
    write_topography(tpath, topo; origin=m.origin)
    topo2 = read_topography(tpath)
    @test topo2.x ≈ topo.x
    @test topo2.y ≈ topo.y
    @test topo2.depth ≈ topo.depth
    @test topo2.z ≈ topo.z
    # regression: reconstruction must not mix the grid frame with origin[3]
    @test air_mask_from_topography(m, topo2) == expected

    #---------- padding blends leave tagged cells alone ----------
    ix, iy = 3:4, 2:4
    kz = 1:nz
    before = copy(m.A)
    MTGeophysics.smooth_padding_decay_z!(m, ix, iy, kz, 2.0, expected)
    MTGeophysics.smooth_padding_decay_xy!(m, ix, iy, 2.0, 10.0, expected)
    @test all(isnan, m.A[expected])
    @test m.A[expected] |> length == 21
    # unprotected padding cells were still blended
    @test any(k -> m.A[1, 1, k] != before[1, 1, k], 3:nz)

    #---------- no tags: masks empty, blend unchanged ----------
    Ap = fill(2.0, nx, ny, nz)
    ppath = joinpath(mktempdir(), "plain.ws")
    write_ws3d_model(ppath, dx, dy, dz, Ap, origin)
    mp = load_ws3d_model(ppath)
    @test count(air_mask_from_model(mp)) == 0
    @test isempty(extract_topography(mp).depth)
    m1 = load_ws3d_model(ppath); m2 = load_ws3d_model(ppath)
    MTGeophysics.smooth_padding_decay_xy!(m1, ix, iy, 2.0, 10.0)
    MTGeophysics.smooth_padding_decay_xy!(m2, ix, iy, 2.0, 10.0, falses(nx, ny, nz))
    @test m1.A == m2.A
end

@testset "Padding blends: frozen sources and below-core thirding" begin

    nx, ny, nz = 6, 5, 6
    dx = fill(1000.0, nx); dy = fill(1000.0, ny)
    dz = [100.0, 200.0, 400.0, 3000.0, 9000.0, 27000.0]   # steep grading below the core
    A = fill(0.0, nx, ny, nz)
    A[:, :, 1] .= -0.5                       # a conductive sheet over the whole top layer
    path = joinpath(mktempdir(), "blend.ws")
    write_ws3d_model(path, dx, dy, dz, A)
    m = load_ws3d_model(path)

    ix, iy, kz = 3:4, 2:4, 1:3
    bg = 2.0
    protected = falses(nx, ny, nz)
    protected[:, :, 1] .= true               # the sheet is frozen everywhere

    #---------- a frozen source must not bleed into the padding ----------
    # seawater is finite, so isfinite alone lets the ocean edge paint a
    # conductive halo across dry padding cells the bathymetry never covered
    m2 = load_ws3d_model(path)
    m2.A[1, 1, 1] = bg                       # a dry padding cell next to the sheet
    protected[1, 1, 1] = false
    MTGeophysics.smooth_padding_decay_xy!(m2, ix, iy, bg, 10.0, protected)
    @test m2.A[1, 1, 1] ≈ bg
    protected[1, 1, 1] = true

    # an unprotected source still blends as before
    m3 = load_ws3d_model(path)
    MTGeophysics.smooth_padding_decay_xy!(m3, ix, iy, bg, 10.0)
    @test m3.A[1, 1, 1] < bg

    #---------- below-core carry-down keeps a third per layer, to the model base ----------
    # the schedule is in layers, not metres: dz grades 400 -> 27000 m here, so a
    # physical e-fold would spend itself entirely on the first layer
    m4 = load_ws3d_model(path)
    m4.A[:, :, 3] .= 0.0
    MTGeophysics.smooth_padding_decay_z!(m4, ix, iy, kz, bg)
    for (k, w) in zip(4:6, (1/3, 1/9, 1/27))
        @test m4.A[3, 2, k] ≈ bg * (1 - w)
    end

    # a frozen column is left to the background rather than carried down
    m5 = load_ws3d_model(path)
    m5.A[:, :, 3] .= 0.0
    p5 = falses(nx, ny, nz)
    p5[3, 2, 3] = true
    MTGeophysics.smooth_padding_decay_z!(m5, ix, iy, kz, bg, p5)
    @test all(k -> m5.A[3, 2, k] ≈ bg, 4:6)
    @test m5.A[4, 2, 4] ≈ bg * (1 - 1/3)

    #---------- one wild core-edge cell must not paint a whole padding row ----------
    m6 = load_ws3d_model(path)
    m6.A[:, :, 3] .= 0.0
    m6.A[4, 3, 3] = 9.0                      # a lone spike on the core edge
    MTGeophysics.smooth_padding_decay_xy!(m6, ix, iy, bg, 10.0)
    @test m6.A[4, 3, 3] ≈ 9.0                # the core itself is untouched
    @test all(i -> m6.A[i, 3, 3] < bg, 5:nx) # median source ignores the spike

    #---------- distance is physical, not index steps times a median width ----------
    dxg = [1000.0, 1000.0, 1000.0, 1000.0, 40000.0, 40000.0]   # padding grades hard
    Ag = fill(0.0, nx, ny, nz)
    gpath = joinpath(mktempdir(), "graded.ws")
    write_ws3d_model(gpath, dxg, dy, dz, Ag)
    mg = load_ws3d_model(gpath)
    MTGeophysics.smooth_padding_decay_xy!(mg, 2:4, iy, bg, 1.0)
    # L = 1 core cell = 1000 m, but i=5 is 20.5 km from the core edge, so it must
    # already be at background. index steps x median width would call it 1 cell
    # out, weight exp(-1), and leave it at 1.26 -- a third of the way to the source
    @test mg.A[5, 3, 3] ≈ bg atol=1e-6
    @test mg.A[6, 3, 3] ≈ bg atol=1e-6
end

@testset "Depth-scaled RBF widths" begin
    nx, ny, nz = 8, 7, 6
    dx = fill(1000.0, nx); dy = fill(1000.0, ny)
    dz = [50.0, 100.0, 300.0, 900.0, 2700.0, 8100.0]
    A = fill(2.0, nx, ny, nz)
    path = joinpath(mktempdir(), "tiny2.ws")
    write_ws3d_model(path, dx, dy, dz, A)
    m = load_ws3d_model(path)

    keep = trues(nx, ny, nz)
    keep[4, 4, 1] = false
    keep[4, 4, 6] = false
    rng = MersenneTwister(11)
    r = build_rbf_map(m, 1:nx, 1:ny, 2, rng;
                      sigma_scale=1.0, sigma_scale_deep=3.0, exclude=keep)
    @test length(r.ci) == 2
    q_top = findfirst(==(1), r.ck)
    q_bot = findfirst(==(6), r.ck)
    @test q_top !== nothing && q_bot !== nothing

    N = nx * ny * nz
    for row in 1:N
        @test sum(r.wts[r.ptr[row]:r.ptr[row+1]-1]) ≈ 1.0
    end

    rowid(i, j, k) = i + (j - 1) * nx + (k - 1) * nx * ny
    wof(row, q) = begin
        p = r.ptr[row]:r.ptr[row+1]-1
        t = findfirst(==(q), r.nbrs[p])
        t === nothing ? 0.0 : r.wts[p][t]
    end
    @test wof(rowid(4, 4, 4), q_bot) > 0.9
    @test wof(rowid(4, 4, 2), q_top) > wof(rowid(4, 4, 2), q_bot)

    dv = zeros(nx, ny, nz)
    apply_rbf_map!(dv, r, fill(3.7, 2))
    @test all(x -> isapprox(x, 3.7; atol=1e-12), dv)
end

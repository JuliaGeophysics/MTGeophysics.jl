# Bathymetry extraction / water-mask / depth-weighted RBF placement tests
# (headless, no solver). Builds a tiny synthetic WS3D model through the file
# round-trip so grid conventions match load_ws3d_model exactly.

using Random

@testset "Bathymetry and water masks (Bathymetry3D)" begin

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

# 3D air and water masks
# Author: @pankajkmishra
# Ensures bathymetry and topography are extracted from tagged WS3D models, round-trip through their files and rebuild
# the same masks, that RBF placement avoids masked cells, and that padding blends leave tagged cells alone
# Tiny synthetic WS3D models go through the file round trip, so grid conventions match load_ws3d_model exactly

using Test, Random

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

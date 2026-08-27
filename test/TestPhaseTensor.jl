@testset "Phase tensors and induction vectors" begin

    @testset "invariants" begin
        # 1-D earth: Z = a(1+i) off-diagonal only, so Phi is the identity
        a  = 3.0
        pt = phase_tensor(0.0 + 0im, a*(1 + 1im), -a*(1 + 1im), 0.0 + 0im)
        @test pt.phimin ≈ 1.0
        @test pt.phimax ≈ 1.0
        @test pt.phi1   ≈ 1.0
        @test pt.phi2   ≈ 1.0
        @test pt.beta   ≈ 0.0 atol = 1e-12
        @test pt.ellipticity ≈ 0.0 atol = 1e-12
        @test pt.anisotropy  ≈ 1.0

        # 2-D earth: Phi = diag(2, 1), principal axes along the coordinate axes
        pt2 = phase_tensor(0.0 + 0im, a*(1 + 1im), -a*(1 + 2im), 0.0 + 0im)
        @test pt2.phimax ≈ 2.0
        @test pt2.phimin ≈ 1.0
        @test pt2.beta ≈ 0.0 atol = 1e-12
        @test pt2.azimuth ≈ 0.0 atol = 1e-12
        @test pt2.anisotropy ≈ 2.0
        @test pt2.ellipticity ≈ 1/3
        @test pt2.phidiff ≈ 1.0

        # 3-D: a non-zero skew survives
        pt3 = phase_tensor(0.2 + 0.4im, a*(1 + 1im), -a*(1 + 2im), 0.1 - 0.3im)
        @test abs(pt3.beta) > 1e-6

        @test isnothing(phase_tensor(0.0 + 0im, 0.0 + 1im, 0.0 + 1im, 0.0 + 0im))  # singular Re(Z)
        @test isnothing(phase_tensor(NaN + 0im, 1 + 1im, -1 - 1im, 0.0 + 0im))
    end

    @testset "induction vectors" begin
        Tzx, Tzy = 0.3 + 0.1im, -0.4 + 0.2im
        p = induction_vector(Tzx, Tzy)                          # :parkinson
        w = induction_vector(Tzx, Tzy; convention = :wiese)

        @test p.re == (-real(Tzy), -real(Tzx))                  # (east, north)
        @test p.im == (-imag(Tzy), -imag(Tzx))
        @test w.re == (real(Tzy), real(Tzx))
        @test p.re_mag ≈ hypot(real(Tzx), real(Tzy)) ≈ w.re_mag
        @test p.im_mag ≈ hypot(imag(Tzx), imag(Tzy))
        @test abs(p.re_azim - w.re_azim) ≈ 180.0
        @test isnothing(induction_vector(NaN + 0im, 0.1 + 0im))
    end

    @testset "whole data sets" begin
        d = make_nan_data()
        d.T = [1.0, 10.0]
        d.f = 1.0 ./ d.T
        d.nf = 2
        d.site = ["S01", "S02", "S03"]
        d.ns = 3
        d.loc = [60.0 24.0 0.0; 60.1 24.2 0.0; 60.2 24.1 0.0]
        d.Z = Array{ComplexF64,3}(undef, d.nf, 4, d.ns)
        d.tip = Array{ComplexF64,3}(undef, d.nf, 2, d.ns)
        for ip in 1:d.nf, is in 1:d.ns
            d.Z[ip, 1, is] = 0.0 + 0im
            d.Z[ip, 2, is] = 2.0 + 2.0im
            d.Z[ip, 3, is] = -2.0 - 4.0im
            d.Z[ip, 4, is] = 0.0 + 0im
            d.tip[ip, 1, is] = 0.1 * ip + 0.0im
            d.tip[ip, 2, is] = 0.05 * is + 0.0im
        end

        PT = phase_tensors_from_data(d)
        IV = induction_vectors_from_data(d)
        @test size(PT) == (d.nf, d.ns)
        @test size(IV) == (d.nf, d.ns)
        @test all(!isnothing, PT)
        @test PT[1, 1].phimax ≈ 2.0
        @test has_tipper_data(d)
        @test all(!isnothing, IV)

        d.tip = Array{ComplexF64}(undef, 0, 0, 0)
        @test !has_tipper_data(d)
        @test all(isnothing, induction_vectors_from_data(d))
    end

    @testset "symbol geometry" begin
        xs, ys = MTGeophysics._ellipse_ring(1.0, 2.0, 0.5, 0.25, 0.0; n = 16)
        @test length(xs) == 17
        @test (xs[1], ys[1]) == (xs[end], ys[end])              # closed ring
        @test maximum(ys) ≈ 2.5                                 # major axis points north
        @test maximum(xs) ≈ 1.25

        xs2, _ = MTGeophysics._ellipse_ring(1.0, 2.0, 0.5, 0.25, 0.0; n = 16, kx = 2.0)
        @test maximum(xs2) ≈ 1.5                                # east offsets stretched

        pt = phase_tensor(0.0 + 0im, 2 + 2im, -2 - 4im, 0.0 + 0im)
        A, B = MTGeophysics._ellipse_semiaxes(pt, 1.0, 4.0)
        @test A ≈ 2.0                                           # scale * Lref / 2
        @test B ≈ A * pt.phimin / pt.phimax
        @test MTGeophysics._ellipse_semiaxes((phimax = 0.0,), 1.0, 4.0) == (0.0, 0.0)

        (sx, sy), (hx, hy) = MTGeophysics._arrow_parts(0.0, 0.0, 0.0, 1.0;
                                                       head_frac = 0.25, head_width = 0.4)
        @test sx == [0.0, 0.0]
        @test sy ≈ [0.0, 0.75]                                  # shaft stops at the head base
        @test hy[1] ≈ 1.0                                       # tip at the vector end
        @test length(hx) == 4 && hx[1] == hx[end]               # closed triangle

        ox, oy = MTGeophysics._arrow_outline(0.0, 0.0, 0.0, 1.0)
        @test length(ox) == length(oy) == 6
        @test MTGeophysics._arrow_parts(1.0, 1.0, 0.0, 0.0) == ((Float64[], Float64[]),
                                                                (Float64[], Float64[]))
        @test MTGeophysics._arrow_outline(1.0, 2.0, 0.0, 0.0) == ([1.0], [2.0])

        @test MTGeophysics._median_site_spacing([0.0], [0.0], 1.0) == 1.0
        @test MTGeophysics._median_site_spacing([0.0, 2.0, 4.0], [0.0, 0.0, 0.0], 1.0) ≈ 2.0
        @test MTGeophysics._median_site_spacing([0.0, 2.0], [0.0, 0.0], 0.5) ≈ 1.0
    end

    @testset "GIS export" begin
        d = make_nan_data()
        d.T = [1.0, 10.0]
        d.f = 1.0 ./ d.T
        d.nf = 2
        d.site = ["S01", "S02"]
        d.ns = 2
        d.loc = [60.0 24.0 0.0; 60.1 24.2 0.0]
        d.x = [0.0, 10000.0]
        d.y = [0.0, 10000.0]
        d.z = [0.0, 0.0]
        d.zrot = zeros(d.nf, d.ns)
        d.trot = zeros(d.nf, d.ns)
        d.origin = [60.05, 24.1, 0.0]
        d.responses = ["ZXX", "ZXY", "ZYX", "ZYY", "TX", "TY"]
        d.nr = 6
        d.Z = Array{ComplexF64,3}(undef, d.nf, 4, d.ns)
        d.Zerr = fill(ComplexF64(0.1, 0.0), d.nf, 4, d.ns)
        d.tip = Array{ComplexF64,3}(undef, d.nf, 2, d.ns)
        d.tiperr = fill(ComplexF64(0.02, 0.0), d.nf, 2, d.ns)
        for ip in 1:d.nf, is in 1:d.ns
            d.Z[ip, 1, is] = 0.05 + 0.02im
            d.Z[ip, 2, is] = 2.0 + 2.0im
            d.Z[ip, 3, is] = -2.0 - 4.0im
            d.Z[ip, 4, is] = -0.03 + 0.01im
            d.tip[ip, 1, is] = 0.2 + 0.05im
            d.tip[ip, 2, is] = -0.1 + 0.02im
        end
        d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)

        mktempdir() do dir
            datafile = joinpath(dir, "ptiv.dat")
            write_data_modem(datafile, d)
            outdir = joinpath(dir, "gis")

            @test write_ptiv_gis(datafile; output_dir = outdir) == outdir
            files = readdir(outdir)
            @test count(f -> endswith(f, ".shp"), files) == 2 * d.nf   # PT and IV per period
            @test count(f -> endswith(f, ".prj"), files) == 2 * d.nf
            @test "README.txt" in files

            pt_shp = joinpath(outdir, first(sort(filter(f -> startswith(f, "ptiv_PT_") &&
                                                             endswith(f, ".shp"), files))))
            tbl = MTGeophysics.Shapefile.Table(pt_shp)
            @test length(tbl.site) == d.ns
            @test Set(tbl.site) == Set(d.site)
            @test all(isfinite, tbl.beta)
            @test all(v -> 0 <= v <= 90, tbl.phimax)               # angles, not tangents
            @test occursin("WGS 84", read(joinpath(outdir, "ptiv_PT_T0001.0000.prj"), String))

            iv_tbl = MTGeophysics.Shapefile.Table(joinpath(outdir, "ptiv_IV_T0001.0000.shp"))
            @test Set(iv_tbl.part) == Set(["real", "imag"])
            @test all(==("parkinson"), iv_tbl.conv)

            # raw tangents instead of angles, and a skew cut that drops every site
            raw = write_ptiv_gis(datafile; output_dir = joinpath(dir, "gis_raw"),
                                 as_angle = false)
            raw_tbl = MTGeophysics.Shapefile.Table(joinpath(raw, "ptiv_PT_T0001.0000.shp"))
            @test maximum(raw_tbl.phimax) > 1.0

            empty_dir = write_ptiv_gis(datafile; output_dir = joinpath(dir, "gis_skip"),
                                       skip_beta_above = 0.0)
            @test count(f -> startswith(f, "ptiv_PT_"), readdir(empty_dir)) == 0
        end
    end
end

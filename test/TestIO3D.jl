@testset "3D ModEM I/O" begin

    @testset "Data structures" begin
        data = make_nan_data()
        @test isa(data, Data)
        @test data.ns == 0
        @test data.nf == 0
    end

    @testset "ModEM data writer round-trip" begin
        d = make_nan_data()
        d.T = [0.01, 1.0, 100.0]
        d.f = 1.0 ./ d.T
        d.nf = 3
        d.site = ["S01", "S02"]
        d.ns = 2
        d.loc = [60.0 24.0 100.0; 61.0 25.0 200.0]
        d.x = [1000.0, 2000.0]
        d.y = [-500.0, 500.0]
        d.z = [100.0, 200.0]
        d.responses = ["ZXX", "ZXY", "ZYX", "ZYY", "TX", "TY"]
        d.nr = 6
        d.zrot = zeros(d.nf, d.ns)
        d.trot = d.zrot
        d.origin = [60.5, 24.5, 0.0]

        d.Z = Array{ComplexF64,3}(undef, d.nf, 4, d.ns)
        d.Zerr = Array{ComplexF64,3}(undef, d.nf, 4, d.ns)
        d.tip = Array{ComplexF64,3}(undef, d.nf, 2, d.ns)
        d.tiperr = Array{ComplexF64,3}(undef, d.nf, 2, d.ns)
        for is in 1:d.ns, ip in 1:d.nf
            for ic in 1:4
                d.Z[ip, ic, is] = ComplexF64(ip + 0.1 * ic + 0.01 * is, -(ip + 0.2 * ic))
                d.Zerr[ip, ic, is] = ComplexF64(0.05 * ip, 0.0)
            end
            for ic in 1:2
                d.tip[ip, ic, is] = ComplexF64(0.3 * ip + 0.01 * ic, 0.1 * is)
                d.tiperr[ip, ic, is] = ComplexF64(0.02 * ip, 0.0)
            end
        end
        d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)

        mktempdir() do dir
            file = joinpath(dir, "roundtrip.dat")
            @test write_data_modem(file, d) == file

            d2 = load_data_modem(file)
            @test d2.site == d.site
            @test d2.ns == d.ns
            @test d2.nf == d.nf
            @test d2.nr == d.nr
            @test d2.T ≈ d.T
            @test d2.loc ≈ d.loc
            @test d2.x ≈ d.x
            @test d2.y ≈ d.y
            @test d2.origin ≈ d.origin
            @test d2.Z ≈ d.Z
            @test d2.tip ≈ d.tip
            @test real.(d2.Zerr) ≈ abs.(d.Zerr)
            @test real.(d2.tiperr) ≈ abs.(d.tiperr)
            @test all(isfinite, d2.ρ)
            @test all(isfinite, d2.φ)

            nan_idx = (2, 2, 1)
            d.Z[nan_idx...] = ComplexF64(NaN, NaN)
            file2 = joinpath(dir, "with_nan.dat")
            write_data_modem(file2, d)
            d3 = load_data_modem(file2)
            @test isnan(real(d3.Z[nan_idx...]))
            @test count(z -> isfinite(real(z)), d3.Z) == length(d3.Z) - 1

            file3 = joinpath(dir, "imp_only.dat")
            write_data_modem(file3, d; include_tipper = false)
            d4 = load_data_modem(file3)
            @test all(z -> isnan(real(z)), d4.tip)
            @test d4.responses == ["ZXX", "ZXY", "ZYX", "ZYY"]
        end
    end

end

@testset "3D galvanic distortion" begin

    # small synthetic dataset, 3 periods x 2 sites, full impedance + tipper
    function _make_pred_data()
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
        return d
    end

    # apply a per-site real 2x2 distortion to the impedance rows
    function _distort!(d, Cs)
        for is in 1:d.ns
            c11, c12, c21, c22 = Cs[is]
            for ip in 1:d.nf
                zxx, zxy, zyx, zyy = d.Z[ip, 1, is], d.Z[ip, 2, is], d.Z[ip, 3, is], d.Z[ip, 4, is]
                d.Z[ip, 1, is] = c11 * zxx + c12 * zyx
                d.Z[ip, 2, is] = c11 * zxy + c12 * zyy
                d.Z[ip, 3, is] = c21 * zxx + c22 * zyx
                d.Z[ip, 4, is] = c21 * zxy + c22 * zyy
            end
        end
        return d
    end

    C_true = [(1.2, 0.3, -0.1, 0.8), (0.9, -0.2, 0.15, 1.1)]

    mktempdir() do dir
        pred_file = joinpath(dir, "dpred.dat")
        obs_file = joinpath(dir, "dobs.dat")
        dp = _make_pred_data()
        write_data_modem(pred_file, dp)
        dobs = _distort!(_make_pred_data(), C_true)
        write_data_modem(obs_file, dobs)

        @testset "known C recovery" begin
            fit = chi2_and_rms_distorted(obs_file, pred_file; damping=1e-8)
            @test isa(fit, DistortionFit)
            for is in 1:2
                c11, c12, c21, c22 = C_true[is]
                @test fit.C[1, 1, is] ≈ c11 atol = 1e-3
                @test fit.C[1, 2, is] ≈ c12 atol = 1e-3
                @test fit.C[2, 1, is] ≈ c21 atol = 1e-3
                @test fit.C[2, 2, is] ≈ c22 atol = 1e-3
            end
            # impedance residuals vanish after correction; only tipper misfit remains (zero here)
            @test fit.rms < 1e-4
        end

        @testset "damping=Inf reproduces chi2_and_rms" begin
            fit = chi2_and_rms_distorted(obs_file, pred_file; damping=Inf)
            @test all(fit.C[1, 1, :] .== 1.0)
            @test all(fit.C[2, 2, :] .== 1.0)
            @test all(fit.C[1, 2, :] .== 0.0)
            @test all(fit.C[2, 1, :] .== 0.0)
            @test fit.penalty == 0.0
            s = chi2_and_rms(obs_file, pred_file)
            @test fit.rms == s.rms
            @test fit.χ² == s.χ²
        end

        @testset "finite damping lowers misfit" begin
            fit_inf = chi2_and_rms_distorted(obs_file, pred_file; damping=Inf)
            fit_low = chi2_and_rms_distorted(obs_file, pred_file; damping=0.01)
            @test fit_low.rms <= fit_inf.rms
            @test fit_low.penalty >= 0.0
        end

        @testset "diagnostics file" begin
            fit = chi2_and_rms_distorted(obs_file, pred_file; damping=1e-8)
            out = joinpath(dir, "distortion_best.txt")
            @test write_distortion_file(out, fit; iter=7, rms=fit.rms, damping=1e-8) == out
            lines = readlines(out)
            @test any(l -> occursin("iter=7", l), lines)
            body = filter(l -> !startswith(l, "#"), lines)
            @test length(body) == 2
            # sorted by ||C-I||_F descending: site 1 (C further from I) first
            @test startswith(body[1], "S01")
        end
    end

    @testset "nan handling" begin
        mktempdir() do dir
            dp = _make_pred_data()
            dobs = _distort!(_make_pred_data(), C_true)
            # knock out one obs entry and one pred entry; site 2 impedance all missing
            dobs.Z[1, 2, 1] = ComplexF64(NaN, NaN)
            dp.Z[2, 3, 1] = ComplexF64(NaN, NaN)
            dobs.Z[:, :, 2] .= ComplexF64(NaN, NaN)
            pred_file = joinpath(dir, "dpred.dat")
            obs_file = joinpath(dir, "dobs.dat")
            write_data_modem(pred_file, dp)
            write_data_modem(obs_file, dobs)
            fit = chi2_and_rms_distorted(obs_file, pred_file; damping=1e-8)
            @test all(isfinite, fit.C[:, :, 1])
            # all-nan site falls back to identity
            @test fit.C[:, :, 2] == [1.0 0.0; 0.0 1.0]
            @test isfinite(fit.rms)
        end
    end
end

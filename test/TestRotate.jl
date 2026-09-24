# Data rotation
# Author: @pankajkmishra
# Ensures rotate_data turns the impedance tensor and tipper as R Z Rᵀ and T Rᵀ, keeps the rotation invariants, turns
# station x, y and the header angle for a mesh turn but not for a declination correction (one angle or one per site),
# drops incomplete tensors, and writes and reads back the rotation history of the file

using Test, LinearAlgebra, Random

@testset "Data rotation" begin
    rng = MersenneTwister(11)
    nf, ns = 3, 3
    d = make_nan_data()
    d.T = [0.01, 1.0, 100.0]
    d.f = 1 ./ d.T
    d.nf, d.ns = nf, ns
    d.site = ["JK01", "JK02", "JK03"]
    d.loc = [62.20 25.70 0.0; 62.25 25.75 0.0; 62.30 25.80 0.0]
    d.x, d.y, d.z = [-5000.0, 0.0, 5000.0], [-2000.0, 1000.0, 4000.0], zeros(ns)
    d.origin = [62.25, 25.75, 0.0]
    d.responses, d.nr = ["ZXX", "ZXY", "ZYX", "ZYY", "TX", "TY"], 6
    d.Z = randn(rng, ComplexF64, nf, 4, ns)
    d.Zerr = complex.(0.05 .* abs.(d.Z))
    d.tip = 0.2 .* randn(rng, ComplexF64, nf, 2, ns)
    d.tiperr = fill(complex(0.02), nf, 2, ns)
    d.zrot = zeros(nf, ns)
    d.trot = zeros(nf, ns)
    d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)

    tensor(d, ip, is) = permutedims(reshape(d.Z[ip, :, is], 2, 2))
    R(θ) = [cosd(θ) sind(θ); -sind(θ) cosd(θ)]

    @testset "mesh turn" begin
        r = rotate_data(d, 30.0)
        for is in 1:ns, ip in 1:nf
            Z = tensor(d, ip, is)
            @test tensor(r, ip, is) ≈ R(30) * Z * R(30)'
            @test det(tensor(r, ip, is)) ≈ det(Z)                                  # rotation invariants
            @test r.Z[ip, 2, is] - r.Z[ip, 3, is] ≈ d.Z[ip, 2, is] - d.Z[ip, 3, is]
            @test r.tip[ip, :, is] ≈ vec(transpose(d.tip[ip, :, is]) * R(30)')
        end
        @test hypot.(r.x, r.y) ≈ hypot.(d.x, d.y) && r.x ≈ cosd(30) .* d.x .+ sind(30) .* d.y
        @test r.loc == d.loc && all(==(30.0), r.zrot)
        @test r.rotations == [RotationStep(:mesh, 30.0)]
        back = rotate_data(r, -30.0)
        @test back.Z ≈ d.Z && back.tip ≈ d.tip && back.x ≈ d.x && all(abs.(back.zrot) .< 1e-12)
        @test length(back.rotations) == 2
        @test d.rotations == RotationStep[] && all(iszero, d.zrot)                 # the input is untouched
    end

    @testset "declination" begin
        D = [9.0, 9.5, 10.0]
        r = rotate_data(d, D; kind = :declination)
        @test all(tensor(r, ip, is) ≈ R(-D[is]) * tensor(d, ip, is) * R(-D[is])' for is in 1:ns, ip in 1:nf)
        @test r.x == d.x && r.y == d.y && r.zrot == d.zrot                          # positions and header stay
        @test r.rotations == [RotationStep(:declination, D, d.site)]
        @test_logs (:warn, r"already corrected") rotate_data(r, 1.0; kind = :declination)
        @test rotate_data(d, 9.5; kind = :declination).rotations == [RotationStep(:declination, 9.5)]
        @test_throws ArgumentError rotate_data(d, D)                                # per site is declination only
        @test_throws ArgumentError rotate_data(d, [1.0, 2.0]; kind = :declination)
        @test_throws ArgumentError rotate_data(d, 1.0; kind = :north)
    end

    @testset "incomplete tensors" begin
        gap = deepcopy(d)
        gap.Z[1, 1, 1] = complex(NaN, NaN)
        r = @test_logs (:warn, r"dropped") rotate_data(gap, 20.0)
        @test all(isnan, r.Z[1, :, 1]) && all(isfinite, r.Z[2, :, 1]) && all(isfinite, r.tip[1, :, 1])
        offdiag = deepcopy(d)
        offdiag.Z[:, [1, 4], :] .= complex(NaN, NaN)
        offdiag.tip .= complex(NaN, NaN)
        @test_throws ErrorException rotate_data(offdiag, 20.0)
    end

    @testset "files and history" begin
        mktempdir() do dir
            path = redirect_stdout(() -> write_data_modem(joinpath(dir, "data.obs"), d; sign = -1, units = "[mV/km]/[nT]"), devnull)
            out = redirect_stdout(() -> rotate_data(path, 30.0), devnull)
            @test out == joinpath(dir, "data-r.obs")
            lines = readlines(out)
            @test any(==("> 30.00"), lines) && any(==("> exp(-i\\omega t)"), lines)
            @test any(l -> startswith(l, "#") && endswith(l, "| rotated: mesh +30.00"), lines)
            loaded = redirect_stdout(() -> load_data_modem(out; warn_rotation = false), devnull)
            @test loaded.rotations == [RotationStep(:mesh, 30.0)] && all(==(30.0), loaded.zrot)
            @test loaded.Z ≈ rotate_data(d, 30.0).Z rtol = 1e-5
            # a second rotation of data-r.obs rewrites it and extends the history; the header keeps the mesh angle
            again = redirect_stdout(() -> rotate_data(out, [9.0, 9.5, 10.0]; kind = :declination), devnull)
            @test again == out
            twice = redirect_stdout(() -> load_data_modem(again; warn_rotation = false), devnull)
            @test twice.rotations == [RotationStep(:mesh, 30.0), RotationStep(:declination, [9.0, 9.5, 10.0], d.site)]
            @test all(==(30.0), twice.zrot) && any(==("> 30.00"), readlines(again))
            # the 2D reader keeps the history, and a strike turn adds to it
            @test load_data2d(again).rotations == twice.rotations
        end
    end
end

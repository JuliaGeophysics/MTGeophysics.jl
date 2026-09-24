# 2D strike and data rotation
# Author: @pankajkmishra
# Ensures the phase tensor strike is recovered and resolved across the station line, that rotation makes a 2D tensor
# diagonal-free with errors, positions and header angle carried, that rotating again is a no-op, and that data
# without ZXX and ZYY keep their frame or refuse a manual strike

using Test

@testset "2D strike and rotation" begin
    # a 2D tensor with strike 30°, stations along 120° (across strike), seen in the north frame
    f = [0.1, 1.0, 10.0]
    nf, ns = length(f), 5
    t = collect(-2000.0:1000.0:2000.0)
    R(θ) = [cosd(θ) sind(θ); -sind(θ) cosd(θ)]
    Zs = [[0 2(1 + 1im) * k; -(1 + 2im) * k 0] for k in f, _ in 1:ns]
    Zn = [R(30)' * Z * R(30) for Z in Zs]
    comp(i, j) = [Z[i, j] for Z in Zn]
    err = fill(0.05, nf, ns)
    rp = fill(NaN, nf, ns)
    data = DataFile2D(title = "strike", periods = 1 ./ f, frequencies = f, site_names = ["TK" * lpad(i, 2, '0') for i in 1:ns],
                      receivers = t .* sind(120), x_positions = t .* cosd(120), z_positions = zeros(ns),
                      z_xx = comp(1, 1), z_xx_error = copy(err), z_xy = comp(1, 2), z_xy_error = copy(err),
                      z_yx = comp(2, 1), z_yx_error = copy(err), z_yy = comp(2, 2), z_yy_error = copy(err),
                      rho_xy = rp, phase_xy = rp, rho_yx = rp, phase_yx = rp,
                      latitudes = fill(62.25, ns), longitudes = fill(25.75, ns), origin = [62.25, 25.75])

    e = EstimateStrike2D(data)
    @test e.strike ≈ 30 atol = 1e-8
    @test e.consistency ≈ 1 && e.skew < 1e-8 && e.count == nf * ns
    r = StrikeData2D(data)
    @test r.rotation == 30 && r.data.rotation == 30
    @test r.data.z_xy ≈ [Z[1, 2] for Z in Zs] && r.data.z_yx ≈ [Z[2, 1] for Z in Zs]
    @test maximum(abs, r.data.z_xx) < 1e-12 && maximum(abs, r.data.z_yy) < 1e-12
    @test r.data.receivers ≈ t && maximum(abs, r.data.x_positions) < 1e-9
    @test r.data.z_xy_error ≈ err                                          # c⁴ + s⁴ + 2c²s² = 1
    @test StrikeData2D(data, 30.0).data.z_xy ≈ r.data.z_xy                # manual strike
    @test StrikeData2D(data, 0.0).data === data                           # the file frame already

    # stations along strike instead: the other axis is taken, so the profile stays across strike
    along = DataFile2D(; (k => getfield(data, k) for k in fieldnames(DataFile2D) if k ∉ (:receivers, :x_positions))...,
                       receivers = t .* sind(30), x_positions = t .* cosd(30))
    @test abs(MTGeophysics._wrap_angle(EstimateStrike2D(along).strike - 120, 180)) < 1e-8

    # no ZXX, ZYY: auto keeps the frame, a manual strike cannot rotate
    offdiag = DataFile2D(; (k => getfield(data, k) for k in fieldnames(DataFile2D) if k ∉ (:z_xx, :z_yy))...,
                         z_xx = fill(complex(NaN), nf, ns), z_yy = fill(complex(NaN), nf, ns))
    @test EstimateStrike2D(offdiag) === nothing && StrikeData2D(offdiag).rotation == 0
    @test_throws ErrorException StrikeData2D(offdiag, 30.0)

    mktempdir() do dir
        path = write_data2d(joinpath(dir, "data.obs"), data; full_tensor = true)
        @test occursin("ZXX", read(path, String))
        w = RotateToStrike2D(path)
        @test w.path == joinpath(dir, "data-r.obs") && w.rotation == 30
        @test any(==("> 30.00"), readlines(w.path))                        # the header rotation line
        back = load_data2d(w.path)
        @test back.rotation == 30
        @test back.rotations == [RotationStep(:strike, 30.0)]                  # the history records the turn
        @test back.z_xy ≈ r.data.z_xy rtol = 1e-5
        @test back.receivers ≈ t atol = 1e-3
        @test abs(StrikeData2D(back).rotation) < 1e-6                     # rotating again is a no-op
    end
end

# 3D ASA helpers
# Author: @pankajkmishra
# Checks the pieces of ASA3DMT that do not need the forward solver: the annealing schedule and its index, the
# transpose of the RBF map, the Rademacher probe file, the reader of the solver's J' output and the reheat factors
# A synthetic J' file is written with the same record layout, so the reader is tested without running the solver

using Test, Random, Statistics

@testset "ASA 3D helpers" begin

    D, k_f = 16, 120.0
    m_scale = log(1 / MTGeophysics._ASA_TEMP_RATIO_SCALE)
    c = m_scale * exp(-log(k_f) / D)
    @test MTGeophysics._asa_temperature(1.0, c, 0, D) == 1.0
    @test MTGeophysics._asa_temperature(1.0, c, k_f, D) ≈ MTGeophysics._ASA_TEMP_RATIO_SCALE rtol=1e-10
    T = MTGeophysics._asa_temperature(2.0, c, 37.5, D)
    @test MTGeophysics._asa_index(2.0, T, c, D) ≈ 37.5 rtol=1e-8
    @test MTGeophysics._asa_index(2.0, 3.0, c, D) == 0.0

    nx, ny, nz = 7, 6, 5
    path = joinpath(mktempdir(), "asa.ws")
    write_ws3d_model(path, fill(1000.0, nx), fill(1000.0, ny), fill(500.0, nz), fill(2.0, nx, ny, nz))
    m = load_ws3d_model(path)
    rng = MersenneTwister(3)
    rbf = build_rbf_map(m, 1:nx, 1:ny, 12, rng; sigma_scale=1.5)
    d = randn(rng, length(rbf.ci))
    g = randn(rng, nx, ny, nz)
    field = zeros(nx, ny, nz)
    apply_rbf_map!(field, rbf, d)
    @test sum(field .* g) ≈ sum(d .* MTGeophysics._asa_rbf_transpose(rbf, g)) rtol=1e-10

    dir = mktempdir()
    dobs = joinpath(dir, "obs.dat")
    write(dobs, """
    # header
    > Full_Impedance
    1.0E+01 S01 1.0 2.0 0.0 0.0 0.0 ZXY    1.0E+00    2.0E+00    2.0E-01
    1.0E+01 S01 1.0 2.0 0.0 0.0 0.0 ZYX   -1.0E+00   -2.0E+00    5.0E-01
    """)
    probe = MTGeophysics._asa_write_probe(dobs, joinpath(dir, "probe.dat"), MersenneTwister(1))
    lines = readlines(probe)
    @test lines[1] == "# header" && lines[2] == "> Full_Impedance"
    v3 = parse.(Float64, split(lines[3])[9:11])
    v4 = parse.(Float64, split(lines[4])[9:11])
    @test abs.(v3[1:2]) ≈ [5.0, 5.0] && v3[3] ≈ 0.2
    @test abs.(v4[1:2]) ≈ [2.0, 2.0] && v4[3] ≈ 0.5

    jt = joinpath(dir, "jt.bin")
    fields = [randn(nx, ny, nz) for _ in 1:3]
    open(jt, "w") do io
        rec(b) = (write(io, Int32(length(b))); write(io, b); write(io, Int32(length(b))))
        rec(Vector{UInt8}(rpad("JT multi_Tx vectors", 80)))
        rec(Vector{UInt8}(rpad("JT multi_Tx vectors", 80)))
        rec(reinterpret(UInt8, Int32[3]))
        for f in fields
            rec(Vector{UInt8}(rpad("LOGE", 80)))
            rec(reinterpret(UInt8, Int32[nx, ny, nz]))
            rec(reinterpret(UInt8, fill(1000.0, nx)))
            rec(reinterpret(UInt8, fill(1000.0, ny)))
            rec(reinterpret(UInt8, fill(500.0, nz)))
            rec(reinterpret(UInt8, [1e-10]))
            rec(reinterpret(UInt8, vec(f)))
        end
    end
    @test MTGeophysics._asa_read_jt(jt, nx, ny, nz) ≈ sum(fields)
    @test_throws ErrorException MTGeophysics._asa_read_jt(jt, nx + 1, ny, nz)

    h = vcat(collect(range(1.0, 1000.0; length=99)), [0.0])
    w, href = MTGeophysics._asa_reheat(h)
    @test href ≈ quantile(h[1:99], MTGeophysics._ASA_REF_QUANTILE)
    @test all(1.0 .<= w .<= MTGeophysics._ASA_REHEAT_MAX)
    @test w[end] == MTGeophysics._ASA_REHEAT_MAX
    @test w[1] == min(href / 1.0, MTGeophysics._ASA_REHEAT_MAX)
end

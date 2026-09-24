# 2D control files
# Author: @pankajkmishra
# Ensures fwd.ctrl, inv.ctrl, the VFSA control, cov.ctrl and mask.ctrl round-trip, that each file holds only its own
# algorithm's keys, that unknown, duplicate, missing and misplaced keys are errors, and that the shipped controls read

using Test

@testset "2D control files" begin
    mktempdir() do dir
        fwd = FwdCtrl2D(mode = :TE, air_layers = 7, air_thickness = 30_000.0, air_growth = 1.5,
                        air_resistivity = 1e8, write_frechet = true)
        @test ReadFwdCtrl2D(WriteFwdCtrl2D(joinpath(dir, "fwd.ctrl"), fwd)) == fwd
        @test fwd.strike === nothing && occursin(r"Strike \(deg\)\s+: auto", read(joinpath(dir, "fwd.ctrl"), String))
        fixed = FwdCtrl2D(mode = :TETM, air_layers = 7, air_thickness = 30_000.0, air_growth = 1.5,
                          air_resistivity = 1e8, strike = 32.5)
        @test ReadFwdCtrl2D(WriteFwdCtrl2D(joinpath(dir, "fixed.ctrl"), fixed)).strike == 32.5
        inv = InvCtrl2D(algorithm = :nlcg, lambda = 3.0, target_rms = 1.2, max_iter = 40,
                        smooth_y = 2.0, nlcg_precondition = false)
        @test ReadInvCtrl2D(WriteInvCtrl2D(joinpath(dir, "inv.ctrl"), inv)) == inv
        # bounds are a VFSA setting: GN and NLCG are unbounded and reject the key
        @test !occursin("bounds", read(joinpath(dir, "inv.ctrl"), String))
        write(joinpath(dir, "bounded.ctrl"), read(joinpath(dir, "inv.ctrl"), String) * "Log10 resistivity bounds : 0 4\n")
        @test_throws ErrorException ReadInvCtrl2D(joinpath(dir, "bounded.ctrl"))
        vfsa = VFSACtrl2D(target_rms = 1.1, max_iter = 9, seed = 7, chains = 3, rbf_top = 1.5, core_expansion = 2, log_bounds = (0.5, 3.5))
        @test ReadVFSACtrl2D(WriteVFSACtrl2D(joinpath(dir, "vfsa.ctrl"), vfsa)) == vfsa
        mask = [1 1 0 1; 1 9 9 1; 1 1 1 0]
        @test ReadMask2D(WriteMask2D(joinpath(dir, "mask.ctrl"), mask)) == mask
        @test readlines(joinpath(dir, "mask.ctrl"))[1] == "4 3"
        # each file holds only its own algorithm's keys
        @test !occursin("VFSA", read(joinpath(dir, "inv.ctrl"), String)) && !occursin("GN damping", read(joinpath(dir, "inv.ctrl"), String))
        @test !any(k -> occursin(k, read(joinpath(dir, "vfsa.ctrl"), String)), ("Smallness", "Algorithm", "lambda"))
        cov = Cov2D(sy = [0.1, 0.2, 0.3], sz = 0.4, n_smooth = 2, exceptions = [(2, 1, 0.0)], mask = [1 1 0 1; 1 2 2 1; 1 1 1 1])
        back = ReadCov2D(WriteCov2D(joinpath(dir, "cov.ctrl"), cov))
        @test back.sy == cov.sy && back.sz == cov.sz && back.n_smooth == 2
        @test back.exceptions == cov.exceptions && back.mask == cov.mask
        @test !any(f -> occursin('#', read(joinpath(dir, f), String)), ("fwd.ctrl", "inv.ctrl", "vfsa.ctrl", "cov.ctrl", "mask.ctrl"))

        # only the required keys, the rest take defaults
        path = joinpath(dir, "min.ctrl")
        write(path, "Algorithm : GN\nInitial damping factor lambda : 1\nExit search when rms is less than : 1.05\n" *
                    "Maximum number of iterations : 5   # trailing comment\n")
        c = ReadInvCtrl2D(path)
        @test c.algorithm == :gn && c.max_iter == 5 && c.target_rms == 1.05 && c.mode == :TETM

        write(path, "Algorithm : VFSA\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\n" *
                    "Maximum number of iterations : 5\n")
        @test_throws ErrorException ReadInvCtrl2D(path)                       # vfsa has its own control
        write(path, "Exit search when rms is less than : 1\nMaximum number of iterations : 5\nVFSA chains : 2\n")
        @test_throws ErrorException ReadVFSACtrl2D(path)                      # old prefixed keys
        write(path, "1 2\n1\n")
        @test_throws ErrorException ReadMask2D(path)                          # one value short
        write(path, "Algorithm : GN\nMaximum number of iterations : 5\n")
        @test_throws ErrorException ReadInvCtrl2D(path)                       # missing required keys
        write(path, "Mode : TE\nAir layers : 5\nAir thickness (m) : 1e4\nAir growth factor : 2\n" *
                    "Air resistivity (ohm m) : 1e9\nAir colour : blue\n")
        @test_throws ErrorException ReadFwdCtrl2D(path)                       # unknown key
        write(path, "Mode : TE\nAir layers : 5\nAir thickness (m) : 1e4\nAir growth factor : 2\n")
        @test_throws ErrorException ReadFwdCtrl2D(path)                       # no default air resistivity
        write(path, "Mode : TE\nMode : TM\n")
        @test_throws ErrorException ReadFwdCtrl2D(path)                       # duplicate key
    end

    # the shipped controls read cleanly, one per algorithm
    ctrl_dir = joinpath(dirname(@__DIR__), "examples", "ctrl", "2D")
    @test ReadFwdCtrl2D(joinpath(ctrl_dir, "FwdCtrl")).mode == :TETM
    @test ReadFwdCtrl2D(joinpath(ctrl_dir, "FwdCtrl")).strike === nothing
    @test [ReadInvCtrl2D(joinpath(ctrl_dir, "InvCtrl.$a")).algorithm for a in ("GN", "NLCG")] == [:gn, :nlcg]
    @test ReadVFSACtrl2D(joinpath(ctrl_dir, "InvCtrl.VFSA")).chains >= 1
end

# Data masks and the data dashboard
# Author: @pankajkmishra
# Ensures a mask written from a keep array reads back to the same keep array (single data, whole periods, whole sites),
# that applying it to a ModEM file drops exactly the masked lines and rewrites the block counts, that EDIs come back
# masked in <dir>_masked_<stamp> with their ModEM file, and that the headless dashboard builds with and without the map
# and the full tensor and masks what a box covers

using Test, Random

@testset "Data masks" begin
    rng = MersenneTwister(5)
    nf, ns = 4, 3
    d = make_nan_data()
    d.T = [0.01, 0.1, 1.0, 10.0]
    d.f = 1 ./ d.T
    d.nf, d.ns = nf, ns
    d.site = ["MK01", "MK02", "MK03"]
    d.loc = [62.20 25.70 10.0; 62.25 25.75 20.0; 62.30 25.80 30.0]
    d.x, d.y, d.z = [-5000.0, 0.0, 5000.0], [-2000.0, 1000.0, 4000.0], [10.0, 20.0, 30.0]
    d.origin = [62.25, 25.75, 0.0]
    d.responses, d.nr = ["ZXX", "ZXY", "ZYX", "ZYY", "TX", "TY"], 6
    d.Z = 1e-3 .* randn(rng, ComplexF64, nf, 4, ns)
    d.Zerr = complex.(0.05 .* abs.(d.Z))
    d.tip = 0.2 .* randn(rng, ComplexF64, nf, 2, ns)
    d.tiperr = fill(complex(0.02), nf, 2, ns)
    d.zrot, d.trot = zeros(nf, ns), zeros(nf, ns)
    d.ρ, d.φ, d.ρerr, d.φerr = calc_rho_pha(d.Z, d.Zerr, d.T)

    keep = trues(nf, 6, ns)
    keep[2, 3, 1] = false          # one datum
    keep[4, :, 2] .= false         # a whole period
    keep[:, :, 3] .= false         # a whole site
    masked = count(!, keep)

    mktempdir() do dir
        mask = joinpath(dir, "mask.txt")
        write_data_mask(mask, d, keep)
        m = read_data_mask(mask)
        @test length(m.entries) == 3
        @test m.zrot == 0.0
        @test mask_keep(d, mask) == keep
        bad = joinpath(dir, "bad.txt")
        write(bad, "MK01 0.1 Zqq\n")
        @test_throws ErrorException read_data_mask(bad)

        e = deepcopy(d)
        @test apply_data_mask!(e, mask) == masked
        @test !isfinite(e.Z[2, 3, 1]) && isfinite(e.Z[2, 2, 1]) && !any(isfinite, e.tip[:, :, 3])

        modem = joinpath(dir, "survey.dat")
        write_data_modem(modem, d)
        out = apply_data_mask(modem, mask; output = joinpath(dir, "survey_masked.dat"))
        a, b = load_data_modem(modem), load_data_modem(out)
        @test count(isfinite, a.Z) + count(isfinite, a.tip) - masked == count(isfinite, b.Z) + count(isfinite, b.tip)
        @test b.ns == 2 && b.site == ["MK01", "MK02"]
        ok = isfinite.(b.Z[:, :, 1])
        @test count(!, ok) == 1 && !ok[2, 3]
        @test b.Z[:, :, 1][ok] ≈ a.Z[:, :, 1][ok]

        edis = ModEMToEDI(modem, joinpath(dir, "EDI"))
        @test length(edis) == ns
        r = apply_data_mask(joinpath(dir, "EDI"), mask)
        @test occursin(r"EDI_masked_\d{8}_\d{6}$", r.edi_dir)
        @test sort(readdir(r.edi_dir)) == ["MK01.edi", "MK02.edi"]
        c = load_data_modem(r.modem)
        @test count(isfinite, c.Z) + count(isfinite, c.tip) == count(isfinite, b.Z) + count(isfinite, b.tip)

        w = DataDashboard(modem; mask, interactive = false, figsize = (1200, 700))
        @test count(w.state.keep) == length(keep) - masked
        @test w.map_axis() !== nothing && length(w.tipper_axes()) == 2
        for (map, tip) in ((false, true), (false, false), (true, false), (true, true))
            w.set_map!(map); w.set_tipper!(tip)
            @test (w.map_axis() !== nothing) == map
            @test length(w.tipper_axes()) == (tip ? 2 : 0)
        end
        w.goto!(1)
        # the diagonals are hidden, yet a band masks all four impedances; Zyx at 0.1 s was masked already
        @test w.mask_band!(w.rho_axis, 0.005, 20.0, false) == 4 * nf - 1
        @test !any(w.state.keep[:, 1:4, 1]) && all(w.state.keep[:, 5:6, 1])
        w.full_toggle.active[] = true
        @test w.mask_band!(w.phase_axis, 0.05, 0.5, true) == 4                            # all four back at 0.1 s
        @test all(w.state.keep[2, 1:4, 1])
        @test w.mask_band!(w.tipper_axes()[2], 0.005, 0.05, false) == 2                  # Tzx and Tzy at 0.01 s

        # a click on a data panel, through the window's own mouse handler, masks nothing
        M = MTGeophysics.Makie
        M.update_state_before_display!(w.fig)
        h = first(x for x in w.hits if x.ax === w.rho_axis && x.ic == 2)
        p = MTGeophysics._db_px(w.rho_axis, h.x[2], h.y[2])
        before = copy(w.state.keep)
        ev = M.events(w.fig)
        ev.mouseposition[] = (p[1] + 1.0, p[2])
        ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.press)
        ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.release)
        @test w.state.keep == before

        # a drag through the same handler masks all four impedances in the band, with the diagonals hidden
        w.full_toggle.active[] = false
        w.state.keep[:, :, 1] .= true
        M.update_state_before_display!(w.fig)
        vp = w.rho_axis.scene.viewport[]
        ev.mouseposition[] = (vp.origin[1] + 2.0, vp.origin[2] + vp.widths[2] / 2)
        ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.press)
        ev.mouseposition[] = (vp.origin[1] + vp.widths[1] - 2.0, vp.origin[2] + 10.0)
        ev.mousebutton[] = M.MouseButtonEvent(M.Mouse.left, M.Mouse.release)
        @test !any(w.state.keep[:, 1:4, 1]) && all(w.state.keep[:, 5:6, 1])

        # every write gets a new time-stamped name beside the data
        dat = MTGeophysics._db_write_modem(w.state, dir, "survey")
        edi = MTGeophysics._db_write_edi(w.state, dir, "survey")
        @test occursin(r"survey_\d{8}_\d{6}\.dat$", dat) && occursin(r"survey_EDI_\d{8}_\d{6}$", edi)
        wd = load_data_modem(dat)
        @test count(isfinite, wd.Z) == count(w.state.keep[:, 1:4, :])
        ref = a.Zerr[[findfirst(t -> isapprox(t, T; rtol = 1e-5), a.T) for T in wd.T], :, indexin(wd.site, a.site)]
        @test abs.(wd.Zerr[isfinite.(wd.Z)]) ≈ abs.(ref[isfinite.(wd.Z)]) rtol = 1e-5                # recorded errors, no floor
        @test sort(readdir(edi)) == ["MK01.edi", "MK02.edi"]

        # a click on the map goes to the site under it, through the window's own mouse handler
        w.set_map!(true)
        MTGeophysics.Makie.update_state_before_display!(w.fig)
        q = MTGeophysics._db_px(w.map_axis(), d.loc[3, 2], d.loc[3, 1])
        ev = MTGeophysics.Makie.events(w.fig)
        ev.mouseposition[] = (q[1] + 2.0, q[2] - 2.0)
        ev.mousebutton[] = MTGeophysics.Makie.MouseButtonEvent(MTGeophysics.Makie.Mouse.left, MTGeophysics.Makie.Mouse.press)
        ev.mousebutton[] = MTGeophysics.Makie.MouseButtonEvent(MTGeophysics.Makie.Mouse.left, MTGeophysics.Makie.Mouse.release)
        @test w.site() == 3
    end
end

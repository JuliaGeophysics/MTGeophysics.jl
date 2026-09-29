# MT data dashboard window
# Author: @pankajkmishra
# The window of DataDashboard (panels and hit tests in DataDashboard.jl): a header to page through the sites, a row of
# tools, and the body in five columns: the site map, the button that collapses it, ρa over φ, the button that
# collapses the tipper, and Tzx over Tzy. The ρa and φ axes and both buttons stay; the map and tipper axes are made
# when their side opens and deleted when it closes. No axis zooms or pans, so every click and drag edits the mask

const _DB_BAND = (:grey45, 0.18)
const _DB_EDGE_PX = 22

_db_freeze!(ax) = foreach(k -> deactivate_interaction!(ax, k), collect(keys(interactions(ax))))

function _dashboard_window(st, stem::AbstractString; export_dir::AbstractString, mask_path::AbstractString,
                           figsize, show_map::Bool, show_tipper::Bool, full_tensor::Bool)
    d = st.obs
    tipper = any(isfinite, d.tip)
    fig = Figure(size = figsize)

    header = fig[1, 1] = GridLayout(halign = :left, tellwidth = false)
    b_first = Button(header[1, 1], label = "|<")
    b_prev = Button(header[1, 2], label = "< Prev")
    b_next = Button(header[1, 3], label = "Next >")
    b_last = Button(header[1, 4], label = ">|")
    menu_sites = Menu(header[1, 5], options = [(s, i) for (i, s) in enumerate(d.site)], default = d.site[1], width = 170)
    info = Label(header[1, 6], "", fontsize = 13, halign = :left)
    items = Any[MarkerElement(; _MT_MARKER..., color = _DB_ZCOLOURS[ic]) for ic in (2, 3, 1, 4)]
    labels = ["Zxy", "Zyx", "Zxx", "Zyy"]
    if tipper
        append!(items, [MarkerElement(; _MT_MARKER..., color = _DB_RE.colour), MarkerElement(; _MT_MARKER..., color = _DB_IM.colour)])
        append!(labels, ["Re T", "Im T"])
    end
    st.pred === nothing || (push!(items, LineElement(color = :black, linewidth = _MT_LINEWIDTH)); push!(labels, "predicted"))
    push!(items, MarkerElement(; _DB_MASKED..., marker = :circle)); push!(labels, "masked")
    Legend(header[1, 7], items, labels; orientation = :horizontal, framevisible = false, labelsize = 12, patchsize = (14, 12))

    tools = fig[2, 1] = GridLayout(halign = :left, tellwidth = false)
    full_toggle = Toggle(tools[1, 1], active = full_tensor)
    Label(tools[1, 2], "Full tensor", fontsize = 12)
    b_site = Button(tools[1, 3], label = "Mask site")
    b_reset = Button(tools[1, 4], label = "Restore site")
    b_mask = Button(tools[1, 5], label = "Save mask")
    b_modem = Button(tools[1, 6], label = "Write ModEM")
    b_edi = Button(tools[1, 7], label = "Write EDI")
    b_png = Button(tools[1, 8], label = "Save PNG")
    hint = "Drag across a panel to mask a band of periods; Shift-drag restores"
    status = Label(tools[1, 9], hint, fontsize = 11, halign = :left, width = 560)

    body = fig[3, 1] = GridLayout(2, 5)
    axρ = _db_period_axis(body[1, 3]; yscale = log10, ylabel = "Apparent resistivity (Ω·m)",
                          xticklabelsvisible = false, xlabelvisible = false)
    axφ = _db_period_axis(body[2, 3]; ylabel = "Phase (°)", yticks = -180:90:180)
    edge(col) = Button(body[1:2, col]; label = "", width = _DB_EDGE_PX, height = Relative(1), tellheight = false,
                       cornerradius = 2, buttoncolor = RGBf(0.93, 0.93, 0.93), fontsize = 26, padding = (0, 0, 0, 0))
    b_map, b_tip = edge(2), edge(4)
    foreach(_db_freeze!, (axρ, axφ))

    ok = isfinite.(d.loc[:, 1]) .& isfinite.(d.loc[:, 2])
    lat0 = any(ok) ? mean(d.loc[ok, 1]) : 0.0
    mx, my = d.loc[:, 2], d.loc[:, 1]
    # the sites plus 8 % on each side; the frame takes the shape of that extent at true scale
    function pad(v)
        lo, hi = isempty(v) ? (-1.0, 1.0) : extrema(v)
        m = max(0.08 * (hi - lo), 0.01)
        (lo - m, hi + m)
    end
    lon_lim, lat_lim = pad(mx[ok]), pad(my[ok])
    map_aspect = (lon_lim[2] - lon_lim[1]) * cosd(lat0) / (lat_lim[2] - lat_lim[1])
    site_fill = Observable(fill(RGBAf(0.2, 0.2, 0.2, 1), d.ns))
    current = Observable(Point2f(mx[1], my[1]))

    cur = Ref(1)
    hits = Any[]
    map_ax = Ref{Any}(nothing)
    tips = Any[]
    map_open, tip_open = Ref(show_map && any(ok)), Ref(show_tipper)

    function refresh_map!()
        map_ax[] === nothing && return
        site_fill[] = [any(st.keep[:, :, is]) ? RGBAf(0.2, 0.2, 0.2, 1) : RGBAf(1, 1, 1, 0) for is in 1:d.ns]
        current[] = Point2f(mx[cur[]], my[cur[]])
    end

    function refresh_info!()
        is = cur[]
        lat, lon, elev = d.loc[is, :]
        nT = count(ip -> any(isfinite, d.Z[ip, :, is]) || any(isfinite, d.tip[ip, :, is]), 1:d.nf)
        nkept, nall = count(st.keep[:, :, is]), count(isfinite, d.Z[:, :, is]) + count(isfinite, d.tip[:, :, is])
        s = @sprintf("%s  (%d/%d)   %.4f°, %.4f°, %.0f m   %d periods   kept %d/%d", d.site[is], is, d.ns, lat, lon, elev, nT, nkept, nall)
        st.pred === nothing || (s *= @sprintf("   RMS %.2f (all %.2f)", _db_site_rms(st, is).rms, _db_total_rms(st)))
        tipper || (s *= "   no tipper in this survey")
        info.text[] = s
    end

    function redraw!()
        empty!(hits)
        _db_impedance!(axρ, axφ, st, cur[], full_toggle.active[] ? (2, 3, 1, 4) : (2, 3), hits)
        for (j, ax) in enumerate(tips)
            _db_tipper!(ax, st, cur[], j, hits)
        end
        refresh_map!()
        refresh_info!()
    end

    # open or close the side panels; the columns stay, a closed one has no width
    function sides!()
        map_ax[] === nothing || (delete!(map_ax[]); map_ax[] = nothing)
        foreach(delete!, tips); empty!(tips)
        if map_open[]
            m = _mt_axis(body[1:2, 1]; tellheight = false, valign = :top, aspect = AxisAspect(map_aspect),
                         limits = (lon_lim, lat_lim),
                         xlabel = "Longitude (°)", ylabel = "Latitude (°)", xticks = WilkinsonTicks(3), yticks = WilkinsonTicks(4),
                         xtickformat = _plain_tickformat, ytickformat = _plain_tickformat)
            _db_freeze!(m)
            scatter!(m, mx, my; color = site_fill, strokecolor = :black, strokewidth = 0.8, markersize = 7)
            scatter!(m, current; color = :transparent, strokecolor = :red, strokewidth = 2.5, markersize = 18)
            map_ax[] = m
        end
        if tip_open[]
            tzx = _db_period_axis(body[1, 5]; ylabel = "Tzx", xticklabelsvisible = false, xlabelvisible = false)
            tzy = _db_period_axis(body[2, 5]; ylabel = "Tzy")
            foreach(_db_freeze!, (tzx, tzy))
            append!(tips, [tzx, tzy])
        end
        # map, ρa/φ and tipper share the width 1 : 2 : 2, or equally when one side is closed
        colsize!(body, 1, map_open[] ? Auto(tip_open[] ? 0.5 : 1.0) : Fixed(0))
        colsize!(body, 3, Auto(1))
        colsize!(body, 5, tip_open[] ? Auto(1) : Fixed(0))
        b_map.label[] = map_open[] ? "‹" : "›"
        b_tip.label[] = tip_open[] ? "›" : "‹"
        redraw!()
    end

    syncing = Ref(false)
    function goto!(is::Integer)
        cur[] = clamp(is, 1, d.ns)
        syncing[] = true
        menu_sites.i_selected[] = cur[]
        syncing[] = false
        status.text[] = hint
        redraw!()
        cur[]
    end

    on(b_first.clicks) do _; goto!(1); end
    on(b_prev.clicks) do _; goto!(cur[] == 1 ? d.ns : cur[] - 1); end
    on(b_next.clicks) do _; goto!(cur[] % d.ns + 1); end
    on(b_last.clicks) do _; goto!(d.ns); end
    on(menu_sites.selection) do i
        syncing[] || i === nothing || i == cur[] || goto!(i)
    end
    on(_ -> redraw!(), full_toggle.active)
    on(b_map.clicks) do _
        any(ok) ? (map_open[] = !map_open[]; sides!()) : (status.text[] = "no site has coordinates")
    end
    on(b_tip.clicks) do _; tip_open[] = !tip_open[]; sides!(); end

    function set_site!(value::Bool)
        _db_set!(st, cur[], [(1, ip) for ip in 1:d.nf], value, 1:6)
        redraw!()
    end
    on(_ -> set_site!(false), b_site.clicks)
    on(_ -> set_site!(true), b_reset.clicks)
    on(b_mask.clicks) do _
        path = isempty(mask_path) ? joinpath(pwd(), "$(stem)_mask.txt") : mask_path
        try
            write_data_mask(path, d, st.keep; source = basename(d.name))
            status.text[] = "Wrote $path"
        catch e
            status.text[] = "Save failed: $(sprint(showerror, e))"
        end
    end
    for (button, write, what) in ((b_modem, _db_write_modem, "ModEM"), (b_edi, _db_write_edi, "EDI"))
        on(button.clicks) do _
            try
                status.text[] = "Wrote " * write(st, export_dir, stem)
            catch e
                status.text[] = "$what export failed: $(sprint(showerror, e))"
            end
        end
    end
    on(b_png.clicks) do _
        path = joinpath(pwd(), "$(stem)-$(d.site[cur[]]).png")
        try
            save(path, fig; px_per_unit = 2, backend = CairoMakie)
            status.text[] = "Saved $path"
        catch e
            status.text[] = "Save failed: $(sprint(showerror, e))"
        end
    end

    # the axes a band spans: ρa and φ together, or both tipper panels
    # on ρa or φ a mask always takes all four impedances, whether the diagonals are shown or not
    group(ax) = ax === axρ || ax === axφ ? [axρ, axφ] : copy(tips)
    comps(ax) = ax === axρ || ax === axφ ? (1:4) : nothing
    function mask_band!(ax, lo, hi, value::Bool)
        picks = _db_in_band(hits, group(ax), lo, hi)
        comps(ax) === nothing || (picks = unique((1, ip) for (_, ip) in picks))
        n = _db_set!(st, cur[], picks, value, comps(ax))
        redraw!()
        n
    end

    # a press in a data panel starts a band; a release close to the press is a click and masks nothing
    drag = Ref{Any}(nothing)
    on(events(fig).mousebutton, priority = 10) do ev
        ev.button == Mouse.left || return Consume(false)
        p = Point2f(events(fig).mouseposition[])
        if ev.action == Mouse.press
            m = map_ax[]
            if m !== nothing && is_mouseinside(m.scene)
                dist = [ok[is] ? norm(_db_px(m, mx[is], my[is]) - p) : Inf for is in 1:d.ns]
                i = argmin(dist)
                dist[i] <= _DB_CLICK_PX && goto!(i)
                return Consume(true)
            end
            k = findfirst(ax -> is_mouseinside(ax.scene), [axρ, axφ, tips...])
            k === nothing && return Consume(false)
            ax = [axρ, axφ, tips...][k]
            x0 = _db_data(ax, p)[1]
            band = Observable((x0, x0))
            plots = [vspan!(a, lift(first, band), lift(last, band); color = _DB_BAND) for a in group(ax)]
            drag[] = (ax = ax, start = p, band = band, plots = plots)
            return Consume(true)
        end
        (ev.action == Mouse.release && drag[] !== nothing) || return Consume(false)
        g = drag[]
        drag[] = nothing
        for (a, pl) in zip(group(g.ax), g.plots)
            pl in a.scene.plots && delete!(a, pl)
        end
        abs(p[1] - g.start[1]) < _DB_DRAG_PX ||
            mask_band!(g.ax, g.band[]..., ispressed(fig, Keyboard.left_shift | Keyboard.right_shift))
        Consume(true)
    end
    on(events(fig).mouseposition, priority = 10) do p
        p = Point2f(p)
        g = drag[]
        g === nothing && return Consume(false)
        g.band[] = (g.band[][1], _db_data(g.ax, p)[1])
        Consume(false)
    end
    on(events(fig).keyboardbutton) do ev
        ev.action in (Keyboard.press, Keyboard.repeat) || return Consume(false)
        ev.key == Keyboard.right && (goto!(cur[] % d.ns + 1); return Consume(true))
        ev.key == Keyboard.left && (goto!(cur[] == 1 ? d.ns : cur[] - 1); return Consume(true))
        Consume(false)
    end

    sides!()
    (fig = fig, goto! = goto!, site = () -> cur[], state = st, hits = hits, full_toggle = full_toggle,
     rho_axis = axρ, phase_axis = axφ, map_axis = () -> map_ax[], tipper_axes = () -> copy(tips),
     set_map! = v -> (map_open[] = v && any(ok); sides!()), set_tipper! = v -> (tip_open[] = v; sides!()),
     mask_band! = mask_band!)
end

"""
    DataDashboard(observed; predicted=nothing, mask=nothing, mask_path="",
                  show_map=true, show_tipper=nothing, full_tensor=false, export_dir=dirname(observed),
                  figsize=(1750, 950), maximize=true,
                  interactive=true, block=!isinteractive(), snapshot_dir="", snapshot_sites=nothing)

Page through the sites of a ModEM data file (convert EDIs with `EDIToModEM` first) and
mask data. `predicted` is a ModEM response, matched to the observed sites by name and to
their periods within 2 %, and turned into the observed frame when the rotations differ.

Each site shows apparent resistivity above phase for Zxy (red) and Zyx (blue), with Zxx
and Zyy (paler) in the same panels when "Full tensor" is on (or `full_tensor = true`);
phases are as recorded, on a fixed -200..200°. The real (red) and imaginary (blue) Tzx
above Tzy sit on the right and a map of the sites in longitude and latitude on the left;
each collapses with the button at its edge (`show_tipper`, default: open when the survey
has a tipper, and `show_map`). The columns share the width 1 : 2 : 2, or equally when one
side is closed. A click on the map picks a site.

No panel zooms, and "Full tensor" only changes the view. Masking is by drag only: a drag
across ρa or φ selects a band of periods over both and masks all four impedances in it (a
drag across the tipper masks Tzx and Tzy), or restores them with Shift held; a click on a
data panel does nothing; "Restore
site" restores the whole site. "Save mask" writes the mask to `mask_path` (default
`mask`, else `<name>_mask.txt` in the working directory), which `apply_data_mask`
applies to the EDIs or to any ModEM file of the survey. `mask` starts from a saved mask.
"Write ModEM" writes the kept data to `<name>_<date_time>.dat` and "Write EDI" as one EDI
per site to `<name>_EDI_<date_time>/`, both with the errors as read and in `export_dir`
(default: beside the data file), so no write replaces an earlier one. The site RMS
against a response uses the same errors; no floor is applied anywhere. Keys: ←/→ change
site.

The window opens maximised to the screen (`maximize = false` keeps `figsize`) and the
layout follows it when resized. `interactive = false` needs no display: it builds the window with CairoMakie and, when
`snapshot_dir` is given, writes one PNG per site (or per name in `snapshot_sites`).
Returns the window's NamedTuple (`fig`, `goto!`, `state`, `set_map!`, `set_tipper!`, ...).
"""
function DataDashboard(observed::AbstractString;
    predicted::Union{Nothing, AbstractString} = nothing,
    mask::Union{Nothing, AbstractString} = nothing,
    mask_path::AbstractString = something(mask, ""),
    show_map::Bool = true,
    show_tipper::Union{Nothing, Bool} = nothing,
    full_tensor::Bool = false,
    export_dir::AbstractString = dirname(abspath(observed)),
    figsize = (1750, 950),
    maximize::Bool = true,
    interactive::Bool = true,
    block::Bool = !isinteractive(),
    snapshot_dir::AbstractString = "",
    snapshot_sites = nothing)

    (isdir(observed) || occursin(r"\.edi$"i, observed)) &&
        error("DataDashboard reads ModEM data; convert the EDIs first: EDIToModEM(\"$observed\")")
    isfile(observed) || error("not found: $observed")
    obs = redirect_stdout(() -> load_data_modem(observed; warn_rotation = false), devnull)
    for A in (obs.Zerr, obs.tiperr), i in eachindex(A)
        abs(A[i]) > 1e10 && (A[i] = complex(NaN, NaN))
    end
    stem = replace(splitext(basename(observed))[1], r"[^A-Za-z0-9]" => "")
    println("Loading ModEM data: $observed")
    @printf("  %d sites, %d periods (%.4g .. %.4g s), components %s\n", obs.ns, obs.nf, minimum(obs.T), maximum(obs.T), join(obs.responses, " "))
    pred = if predicted === nothing
        nothing
    else
        isfile(predicted) || error("predicted file not found: $predicted")
        println("Loading predicted response: $predicted")
        _db_align_predicted(obs, redirect_stdout(() -> load_data_modem(predicted; warn_rotation = false), devnull))
    end

    keep = BitArray(undef, obs.nf, 6, obs.ns)
    keep[:, 1:4, :] .= isfinite.(obs.Z)
    keep[:, 5:6, :] .= isfinite.(obs.tip)
    if mask !== nothing
        keep .&= mask_keep(obs, mask)
        println("  mask $mask: $(count(!, keep) - count(!isfinite, obs.Z) - count(!isfinite, obs.tip)) data masked")
    end
    Trange = (minimum(obs.T) / 1.5, maximum(obs.T) * 1.5)
    st = (obs = obs, pred = pred, keep = keep, Trange = Trange)
    pred === nothing || @printf("  RMS %.3f over all sites (recorded errors)\n", _db_total_rms(st))

    println(interactive ? "Opening the editor" : "Building the editor window headless")
    if interactive
        isdefined(@__MODULE__, :GLMakie) || error("the dashboard window needs GLMakie and a display; interactive = false renders with CairoMakie")
        GLMakie.activate!(title = "MTGeophysics data dashboard: $(basename(observed))")
    else
        CairoMakie.activate!()
    end
    w = _dashboard_window(st, stem; export_dir, mask_path, figsize, show_map, full_tensor,
                          show_tipper = something(show_tipper, any(isfinite, obs.tip)))

    if !isempty(snapshot_dir)
        mkpath(snapshot_dir)
        which = snapshot_sites === nothing ? (1:obs.ns) : [findfirst(==(s), obs.site) for s in snapshot_sites]
        for is in which
            is === nothing && continue
            w.goto!(is)
            save(joinpath(snapshot_dir, "$(stem)-$(obs.site[is]).png"), w.fig; px_per_unit = 1.5, backend = CairoMakie)
        end
        println("  $(length(which)) snapshot(s) in $snapshot_dir")
        w.goto!(1)
    end
    if interactive
        screen = display(w.fig)
        if maximize
            try
                GLMakie.GLFW.MaximizeWindow(screen.glscreen)
            catch e
                @warn "could not maximise the window; it keeps figsize" exception = e
            end
        end
        println("Dashboard ready: ←/→ site, click the map to pick a site, click or drag across a panel to mask")
        block && wait(screen)
    end
    w
end

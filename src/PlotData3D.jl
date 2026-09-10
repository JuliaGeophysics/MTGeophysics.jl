# Interactive viewers for MT data sets 
# Author: @pankajkmishra
# This file holds the period-stepped phase tensor / induction vector map.
# The invariants, the symbol geometry and the shapefile export are headless and
# live in PhaseTensor.jl; only the window and its widgets are here

# ---------- clipped shapefile overlays ----------

"""
    _draw_clipped_coords!(ax, coords, ct, out, lim; color, linewidth, alpha, visible)

Draw one GeoInterface coordinate tree, clipped to `lim`.
In:  axis, nested coordinates, CRS transform, handle sink, clip box, style.
Out: nothing; plot handles are appended to `out`.
"""
function _draw_clipped_coords!(ax, coords, ct, out, lim;
                               color, linewidth, alpha, point_size, visible)
    if _is_xy_tuple(coords) || _is_xy_vector(coords)
        x, y = ct(Float64(coords[1]), Float64(coords[2]))
        _inside_box(x, y, lim) || return nothing
        push!(out, scatter!(ax, [x], [y]; color = (color, alpha),
                            markersize = point_size, visible = visible))
        return nothing
    end
    (coords isa AbstractVector && !isempty(coords)) || return nothing
    if all(c -> _is_xy_tuple(c) || _is_xy_vector(c), coords)
        xy = [ct(Float64(c[1]), Float64(c[2])) for c in coords]
        length(xy) < 2 && return nothing
        for (px, py) in _clip_polyline_to_box(first.(xy), last.(xy), lim)
            push!(out, lines!(ax, px, py; color = (color, alpha),
                              linewidth = linewidth, visible = visible))
        end
        return nothing
    end
    for c in coords
        _draw_clipped_coords!(ax, c, ct, out, lim;
                              color = color, linewidth = linewidth, alpha = alpha,
                              point_size = point_size, visible = visible)
    end
    return nothing
end

"""
    _draw_clipped_shapefiles!(ax, loaded_shapefiles, lim; color, linewidth, alpha, point_size, visible)

Draw every loaded overlay, clipped to `lim` so an oversized shapefile cannot
zoom the map out. Unlike `_draw_all_shapefiles!` it clips and returns handles.
In:  axis, prepared shapefiles, clip box and the style overrides.
Out: vector of plot handles, for restyling or deletion.
"""
function _draw_clipped_shapefiles!(ax, loaded_shapefiles, lim;
                                   color, linewidth, alpha, point_size, visible = true)
    out = Any[]
    for shp in loaded_shapefiles, g in shp.geoms
        coords = try GeoInterface.coordinates(g) catch; nothing end
        isnothing(coords) && continue
        _draw_clipped_coords!(ax, coords, shp.coord_transform, out, lim;
                              color = color, linewidth = linewidth, alpha = alpha,
                              point_size = point_size, visible = visible)
    end
    return out
end

# ---------- phase tensor / induction vector map ----------

"""
    PlotPTIVMap(data_file; crs, shapefiles, ..., interactive)

Interactive map of phase tensor ellipses and induction arrows, one period at a
time, with shapefile overlays, PNG export and GIS export of every period.

In:  path to a ModEM data file (impedance, optionally tipper) and the display
     options; see `examples/plot_PTIV_map.jl` for the full set with defaults.
Out: the Figure, after the window closes.

Needs GLMakie and a display; `write_ptiv_gis` writes the same shapefiles headless.

Ellipses are normalised so every site is the same size and colored by the chosen
invariant. Arrows follow `iv_convention`; Parkinson points towards conductors.
"""
function PlotPTIVMap(data_file::AbstractString;
    crs::AbstractString = "EPSG:4326",
    shapefiles = [],
    shapefile_color = :grey30,
    shapefile_alpha::Real = 0.9,
    shapefile_line_width::Real = 1.2,
    shapefile_point_size::Real = 6,
    pt_fill::Symbol = :beta,
    pt_as_angle::Bool = true,
    pt_scale::Real = 0.66,
    pt_scale_step::Real = 1.25,
    pt_colormap = :Spectral,
    pt_range = nothing,
    pt_stroke = :black,
    pt_strokewidth::Real = 1.1,
    skip_beta_above::Union{Nothing, Real} = nothing,
    iv_convention::Symbol = :parkinson,
    iv_scale::Real = 1.2,
    iv_scale_step::Real = 1.25,
    iv_real_color = :black,
    iv_imag_color = :grey45,
    iv_linewidth::Real = 1.4,
    iv_head_frac::Real = 0.30,
    iv_head_width::Real = 0.42,
    iv_max_magnitude::Real = 1.0,
    show_phase_tensor::Bool = true,
    show_iv_real::Bool = true,
    show_iv_imag::Bool = false,
    show_sites::Bool = true,
    show_site_labels::Bool = false,
    site_color = :grey20,
    site_markersize::Real = 4,
    map_pad::Real = 0.06,
    viewer_figsize = (1250, 950),
    export_dpi::Int = 3,
    export_figsize = (1150, 950),
    gis_output_dir::AbstractString = "",
    interactive::Bool = true)

    isempty(data_file) && error("data_file is required")
    isfile(data_file)  || error("Data file not found: $data_file")
    haskey(PT_FILL_OPTIONS, pt_fill) ||
        error("pt_fill must be one of $(PT_FILL_ORDER), got :$pt_fill")

    println("Loading ModEM data: $data_file")
    d = load_data_modem(data_file)
    println("  sites   : $(d.ns)")
    println("  periods : $(d.nf)  ($(minimum(d.T)) .. $(maximum(d.T)) s)")

    fr = _ptiv_frame(d, crs; pt_scale = pt_scale, iv_scale = iv_scale,
                     map_pad = map_pad, convention = iv_convention)
    fr.has_tipper || println("  tipper  : absent, induction vectors disabled")
    @printf("  median site spacing: %.4g map units\n", fr.Lref)
    n_pt = count(!isnothing, fr.PT)
    @printf("  phase tensors: %d of %d site-periods (%.0f%%)\n",
            n_pt, d.nf*d.ns, 100*n_pt/(d.nf*d.ns))

    data_name = splitext(basename(data_file))[1]
    gis_dir(dir) = isempty(dir) ? joinpath(pwd(), "$(data_name)-PTIV-GIS") : String(dir)
    export_gis(; as_angle, pt_sc, iv_sc) = _export_ptiv_gis(d, fr;
        output_dir = gis_dir(gis_output_dir), data_name = data_name, crs = crs,
        source_file = data_file, as_angle = as_angle, pt_scale = pt_sc, iv_scale = iv_sc,
        iv_convention = iv_convention, iv_max_magnitude = iv_max_magnitude,
        iv_head_frac = iv_head_frac, iv_head_width = iv_head_width,
        skip_beta_above = skip_beta_above)

    loaded_shapefiles = Any[]
    if !isempty(shapefiles)
        println("\nShapefile overlay:")
        append!(loaded_shapefiles, prepare_shapefiles(shapefiles, crs))
    end
    shp_def(path) = (enabled = true, path = String(path), color = shapefile_color,
                     alpha = shapefile_alpha, point_size = shapefile_point_size,
                     line_width = shapefile_line_width)

    site_x, site_y, Lref, kx = fr.site_x, fr.site_y, fr.Lref, fr.lon_stretch
    map_limits = fr.limits
    # same axis wording as the model viewers in GeoRef3D._build_model_in_crs
    xlabel, ylabel = uppercase(strip(crs)) == "MODEL" ? ("Y (m)", "X (m)") :
                     fr.geographic ? ("Longitude (°)", "Latitude (°)") :
                     ("Easting (m)", "Northing (m)")
    keep_site(pt) = isnothing(skip_beta_above) || abs(pt.beta) <= skip_beta_above

    fig = Figure(size = viewer_figsize)
    title_str = Observable("")
    Label(fig[0, 1:2], title_str, fontsize = 18, font = :bold)
    ax = Axis(fig[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        aspect = AxisAspect(fr.aspect),
        limits = map_limits,
        xgridvisible = false,
        ygridvisible = false,
        xtickformat = _plain_tickformat,
        ytickformat = _plain_tickformat)

    cur_period = Observable(1)
    opt_pt     = Observable(show_phase_tensor)
    opt_ivr    = Observable(show_iv_real && fr.has_tipper)
    opt_ivi    = Observable(show_iv_imag && fr.has_tipper)
    opt_sites  = Observable(show_sites)
    opt_labels = Observable(show_site_labels)
    opt_fill   = Observable(pt_fill)
    opt_angle  = Observable(pt_as_angle)
    info_str   = Observable("")

    ell_polys  = Observable(Vector{Point2f}[])
    ell_colors = Observable(Float32[])
    cmap_obs   = Observable(pt_colormap)
    ivr_lines  = Observable(Point2f[])
    ivi_lines  = Observable(Point2f[])
    ivr_heads  = Observable(Vector{Point2f}[])
    ivi_heads  = Observable(Vector{Point2f}[])
    pt_scale_obs = Observable(Float64(pt_scale))
    iv_scale_obs = Observable(Float64(iv_scale))
    crange     = Observable((0.0f0, 1.0f0))
    cbar_label = Observable("")

    opt_shp   = Observable(true)
    shp_color = Observable(shapefile_color)
    shp_width = Observable(Float64(shapefile_line_width))
    shp_plots = Any[]

    function redraw_overlays!()
        for pl in shp_plots
            try; delete!(ax, pl); catch; end
        end
        empty!(shp_plots)
        append!(shp_plots, _draw_clipped_shapefiles!(ax, loaded_shapefiles, map_limits;
                                                     color = shp_color[], linewidth = shp_width[],
                                                     alpha = shapefile_alpha,
                                                     point_size = shapefile_point_size,
                                                     visible = opt_shp))
        limits!(ax, map_limits...)
    end

    poly!(ax, ell_polys; color = ell_colors, colormap = cmap_obs, colorrange = crange,
          strokecolor = pt_stroke, strokewidth = pt_strokewidth, visible = opt_pt)
    lines!(ax, ivr_lines; color = iv_real_color, linewidth = iv_linewidth, visible = opt_ivr)
    poly!(ax, ivr_heads; color = iv_real_color, strokewidth = 0, visible = opt_ivr)
    lines!(ax, ivi_lines; color = iv_imag_color, linewidth = iv_linewidth, visible = opt_ivi)
    poly!(ax, ivi_heads; color = iv_imag_color, strokewidth = 0, visible = opt_ivi)
    scatter!(ax, site_x, site_y; color = site_color, markersize = site_markersize,
             visible = opt_sites)
    text!(ax, site_x, site_y; text = d.site, fontsize = 7, color = site_color,
          align = (:left, :bottom), offset = (3, 3), visible = opt_labels)
    Colorbar(fig[1, 2], colormap = cmap_obs, limits = crange, label = cbar_label, width = 14)

    function refresh!()
        ip, key, asa = cur_period[], opt_fill[], opt_angle[]

        polys = Vector{Point2f}[]; cols = Float32[]; vals = Float64[]
        for is in 1:d.ns
            pt = fr.PT[ip, is]
            (isnothing(pt) || !keep_site(pt)) && continue
            v = _pt_fill_value(pt, key, asa)
            isfinite(v) || continue
            a, b = _ellipse_semiaxes(pt, pt_scale_obs[], Lref)
            a > 0 || continue
            xs, ys = _ellipse_ring(site_x[is], site_y[is], a, b, pt.azimuth; kx = kx)
            push!(polys, [Point2f(xs[k], ys[k]) for k in eachindex(xs)])
            push!(cols, Float32(v)); push!(vals, v)
        end

        if !isnothing(pt_range)
            crange[] = (Float32(pt_range[1]), Float32(pt_range[2]))
        elseif !isempty(vals)
            lo, hi = quantile(vals, 0.02), quantile(vals, 0.98)
            if PT_FILL_OPTIONS[key].symmetric
                m = max(abs(lo), abs(hi)); lo, hi = -m, m
            end
            lo == hi && (lo -= 0.5; hi += 0.5)
            crange[] = (Float32(lo), Float32(hi))
        end
        cmap_obs[]   = pt_colormap
        cbar_label[] = PT_FILL_OPTIONS[key].label *
                       (asa && PT_FILL_OPTIONS[key].angle ? " (deg)" : "")

        lr = Point2f[]; li = Point2f[]
        hr = Vector{Point2f}[]; hi = Vector{Point2f}[]
        if fr.has_tipper
            sc = iv_scale_obs[] * Lref
            for is in 1:d.ns
                iv = fr.IV[ip, is]
                isnothing(iv) && continue
                for (lstore, hstore, vec, mag) in ((lr, hr, iv.re, iv.re_mag),
                                                   (li, hi, iv.im, iv.im_mag))
                    (isfinite(mag) && 0 < mag <= iv_max_magnitude) || continue
                    (sx, sy), (hx, hy) = _arrow_parts(site_x[is], site_y[is],
                                                      vec[1]*sc, vec[2]*sc;
                                                      head_frac = iv_head_frac,
                                                      head_width = iv_head_width, kx = kx)
                    isempty(sx) && continue
                    for k in eachindex(sx); push!(lstore, Point2f(sx[k], sy[k])); end
                    push!(lstore, Point2f(NaN, NaN))
                    push!(hstore, [Point2f(hx[k], hy[k]) for k in eachindex(hx)])
                end
            end
        end
        ivr_lines[] = lr; ivi_lines[] = li
        ivr_heads[] = hr; ivi_heads[] = hi
        ell_polys[] = polys; ell_colors[] = cols

        T = d.T[ip]
        title_str[] = @sprintf("Period %d / %d   |   T = %.4g s   f = %.4g Hz",
                               ip, d.nf, T, T > 0 ? 1/T : NaN)
        nb = count(is -> !isnothing(fr.PT[ip,is]) && abs(fr.PT[ip,is].beta) > 3, 1:d.ns)
        info_str[] = @sprintf("%d ellipses  |  |beta| > 3 deg at %d sites  |  fill: %s  |  PT x%.2f  IV x%.2f",
                              length(vals), nb, PT_FILL_OPTIONS[key].label,
                              pt_scale_obs[] / pt_scale, iv_scale_obs[] / iv_scale)
        return nothing
    end

    refresh!()
    redraw_overlays!()

    function export_figure()
        ex = Figure(size = export_figsize, fontsize = 16)
        Label(ex[1, 1:2], title_str[], fontsize = 20, font = :bold)
        exax = Axis(ex[2, 1], xlabel = xlabel, ylabel = ylabel,
                    aspect = AxisAspect(fr.aspect), limits = map_limits,
                    xtickformat = _plain_tickformat, ytickformat = _plain_tickformat)
        opt_shp[] && _draw_clipped_shapefiles!(exax, loaded_shapefiles, map_limits;
                                               color = shp_color[], linewidth = shp_width[],
                                               alpha = shapefile_alpha,
                                               point_size = shapefile_point_size)
        poly!(exax, ell_polys[]; color = ell_colors[], colormap = cmap_obs[],
              colorrange = crange[], strokecolor = pt_stroke, strokewidth = pt_strokewidth,
              visible = opt_pt[])
        if opt_ivr[]
            lines!(exax, ivr_lines[]; color = iv_real_color, linewidth = iv_linewidth)
            poly!(exax, ivr_heads[]; color = iv_real_color, strokewidth = 0)
        end
        if opt_ivi[]
            lines!(exax, ivi_lines[]; color = iv_imag_color, linewidth = iv_linewidth)
            poly!(exax, ivi_heads[]; color = iv_imag_color, strokewidth = 0)
        end
        opt_sites[] && scatter!(exax, site_x, site_y; color = site_color,
                                markersize = site_markersize)
        Colorbar(ex[2, 2], colormap = cmap_obs[], limits = crange[],
                 label = cbar_label[], width = 14)
        limits!(exax, map_limits...)
        Label(ex[3, 1:2], info_str[], fontsize = 11)
        fn = @sprintf("%s_PTIV_T%09.4f.png", data_name, d.T[cur_period[]])
        # CairoMakie writes the file: letting GLMakie render an off-screen figure
        # creates and destroys a second GL context, taking the viewer window with it
        _save_figure_headless(fn, ex, export_dpi)
        info_str[] = "Figure exported: $(basename(fn))"
        println("Figure exported: $fn")
        return fn
    end

    slider_grid = fig[2, 1:2] = GridLayout()
    btn_first = Button(slider_grid[1, 1], label = "|<< First")
    btn_prev  = Button(slider_grid[1, 2], label = "<< Prev")
    sl        = Slider(slider_grid[1, 3], range = 1:d.nf, startvalue = 1, width = 400)
    per_lbl   = Observable("1 / $(d.nf)")
    Label(slider_grid[1, 4], per_lbl, fontsize = 14)
    btn_next  = Button(slider_grid[1, 5], label = "Next >>")
    btn_last  = Button(slider_grid[1, 6], label = "Last >>|")

    toggle_grid = fig[3, 1:2] = GridLayout()
    tg_pt   = Toggle(toggle_grid[1, 1],  active = opt_pt[])
    Label(toggle_grid[1, 2],  "Phase tensors", fontsize = 12)
    tg_ivr  = Toggle(toggle_grid[1, 3],  active = opt_ivr[])
    Label(toggle_grid[1, 4],  "IV real", fontsize = 12)
    tg_ivi  = Toggle(toggle_grid[1, 5],  active = opt_ivi[])
    Label(toggle_grid[1, 6],  "IV imag", fontsize = 12)
    tg_site = Toggle(toggle_grid[1, 7],  active = opt_sites[])
    Label(toggle_grid[1, 8],  "Sites", fontsize = 12)
    tg_lab  = Toggle(toggle_grid[1, 9],  active = opt_labels[])
    Label(toggle_grid[1, 10], "Labels", fontsize = 12)
    tg_ang  = Toggle(toggle_grid[1, 11], active = opt_angle[])
    Label(toggle_grid[1, 12], "Phi as angle", fontsize = 12)
    Label(toggle_grid[1, 13], "Fill:", fontsize = 12)
    mn_fill = Menu(toggle_grid[1, 14],
                   options = [(PT_FILL_OPTIONS[k].label, k) for k in PT_FILL_ORDER],
                   default = PT_FILL_OPTIONS[pt_fill].label, width = 170)
    Label(toggle_grid[1, 15], "Size:", fontsize = 12)
    btn_pt_dec = Button(toggle_grid[1, 16], label = "PT -")
    btn_pt_inc = Button(toggle_grid[1, 17], label = "PT +")
    btn_iv_dec = Button(toggle_grid[1, 18], label = "IV -")
    btn_iv_inc = Button(toggle_grid[1, 19], label = "IV +")

    button_grid = fig[4, 1:2] = GridLayout()
    tg_shp = Toggle(button_grid[1, 1], active = true)
    Label(button_grid[1, 2], "Overlays", fontsize = 12)
    btn_browse   = Button(button_grid[1, 3], label = "Add Shapefile")
    btn_clear    = Button(button_grid[1, 4], label = "Clear")
    mn_shp_color = Menu(button_grid[1, 5],
                        options = ["grey30", "black", "white", "firebrick", "steelblue"],
                        default = "grey30", width = 115)
    mn_shp_width = Menu(button_grid[1, 6], options = ["thin", "medium", "thick"],
                        default = "medium", width = 95)
    shp_count = Observable(isempty(loaded_shapefiles) ? "none loaded" :
                           "$(length(loaded_shapefiles)) loaded")
    Label(button_grid[1, 7], shp_count, fontsize = 12)
    btn_reset  = Button(button_grid[1, 8],  label = "Reset Zoom")
    btn_export = Button(button_grid[1, 9],  label = "Export Figure")
    btn_gis    = Button(button_grid[1, 10], label = "Export GIS")

    info_grid = fig[5, 1:2] = GridLayout()
    Label(info_grid[1, 1], info_str, fontsize = 12)

    on(sl.value) do v
        cur_period[] = v
        per_lbl[] = "$v / $(d.nf)"
        refresh!()
    end
    on(btn_prev.clicks)  do _; set_close_to!(sl, max(1, sl.value[] - 1)); end
    on(btn_next.clicks)  do _; set_close_to!(sl, min(d.nf, sl.value[] + 1)); end
    on(btn_first.clicks) do _; set_close_to!(sl, 1); end
    on(btn_last.clicks)  do _; set_close_to!(sl, d.nf); end
    on(btn_reset.clicks) do _; limits!(ax, map_limits...); end

    bump!(obs, factor) = (obs[] = clamp(obs[] * factor, 0.02, 20.0); refresh!())
    on(btn_pt_dec.clicks) do _; bump!(pt_scale_obs, 1 / pt_scale_step); end
    on(btn_pt_inc.clicks) do _; bump!(pt_scale_obs, pt_scale_step); end
    on(btn_iv_dec.clicks) do _; bump!(iv_scale_obs, 1 / iv_scale_step); end
    on(btn_iv_inc.clicks) do _; bump!(iv_scale_obs, iv_scale_step); end

    on(tg_pt.active)   do v; opt_pt[]     = v; end
    on(tg_ivr.active)  do v; opt_ivr[]    = v && fr.has_tipper; end
    on(tg_ivi.active)  do v; opt_ivi[]    = v && fr.has_tipper; end
    on(tg_site.active) do v; opt_sites[]  = v; end
    on(tg_lab.active)  do v; opt_labels[] = v; end
    on(tg_ang.active)  do v; opt_angle[]  = v; refresh!(); end
    on(mn_fill.selection) do v; isnothing(v) || (opt_fill[] = v; refresh!()); end

    on(tg_shp.active)  do v; opt_shp[] = v; end
    on(mn_shp_color.selection) do v
        isnothing(v) && return
        shp_color[] = Symbol(v); redraw_overlays!()
    end
    on(mn_shp_width.selection) do v
        isnothing(v) && return
        shp_width[] = v == "thin" ? 0.6 : v == "thick" ? 2.4 : 1.2
        redraw_overlays!()
    end
    on(btn_browse.clicks) do _
        path = String(strip(pick_shapefile()))
        if isempty(path)
            info_str[] = "No file chosen (or no file chooser on this system)."
            return
        elseif !isfile(path)
            info_str[] = "Not a file: $path"
            return
        end
        try
            added = prepare_shapefiles([shp_def(path)], crs)
            if isempty(added)
                info_str[] = "Nothing loaded from $(basename(path))"
            else
                append!(loaded_shapefiles, added)
                redraw_overlays!()
                shp_count[] = "$(length(loaded_shapefiles)) loaded"
                info_str[] = "Added $(basename(path)) ($(length(added[1].geoms)) features)"
            end
        catch e
            showerror(stderr, e, catch_backtrace()); println(stderr)
            info_str[] = "Import failed: $(sprint(showerror, e))"
        end
    end
    on(btn_clear.clicks) do _
        empty!(loaded_shapefiles)
        redraw_overlays!()
        shp_count[] = "none loaded"
        info_str[] = "Overlays cleared."
    end

    on(btn_export.clicks) do _
        try
            export_figure()
        catch e
            showerror(stderr, e, catch_backtrace()); println(stderr)
            info_str[] = "Export failed: $(sprint(showerror, e))"
        end
    end
    on(btn_gis.clicks) do _
        try
            dir = export_gis(; as_angle = opt_angle[], pt_sc = pt_scale_obs[],
                               iv_sc = iv_scale_obs[])
            info_str[] = "GIS exported: $(basename(dir))/"
        catch e
            showerror(stderr, e, catch_backtrace()); println(stderr)
            info_str[] = "GIS export failed: $(sprint(showerror, e))"
        end
    end

    println("\nViewer ready!")
    println("  - Coordinate system: $crs")
    println("  - Use the slider or the Prev/Next buttons to step through periods")
    println("  - PT +/- and IV +/- resize the symbols; the GIS export follows them")
    println("  - Click 'Export Figure' for a high-resolution PNG")
    isempty(loaded_shapefiles) || println("  - $(length(loaded_shapefiles)) shapefile(s) overlaid")

    if interactive
        Makie.update_state_before_display!(fig)
        screen = GLMakie.Screen(fig.scene)
        println("\nClose the figure window to exit...")
        wait(screen)
    end
    return fig
end

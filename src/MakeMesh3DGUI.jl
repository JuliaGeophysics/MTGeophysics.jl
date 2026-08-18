#----- rectangle outline for the core box -----
corner_points(yl, yh, xl, xh) = [Point2f(yl, xl), Point2f(yh, xl), Point2f(yh, xh),
                                 Point2f(yl, xh), Point2f(yl, xl)]

#----- check the mesh and advise what to fix -----
function suggestions(m)
    msgs = String[]
    m.maxper > 1 && push!(msgs,
        "• More than one site falls in some cells — decrease 'Core cell N–S (m)' / 'Core cell E–W (m)'.")
    m.pad_ok || push!(msgs,
        "• Lateral padding is too small — increase 'Pad cells N–S' / 'Pad cells E–W' or 'Pad factor'; " *
        "the boundary should sit at least one skin depth at the longest period (δ(Tmax)) from the core.")
    m.depth_km < m.δmax_km && push!(msgs,
        "• Model is too shallow — increase 'Depth × δ(Tmax)' so the mesh spans the depth of investigation.")
    !isempty(m.site_surface.impossible) && push!(msgs,
        "• Some sites still sit above the shallowest available earth cell — reduce the top air thickness or revise the topography/datum so the earth model starts above those site elevations.")
    m.site_air.ok || push!(msgs, m.site_air.text * " " * m.site_air.fix)
    ok = isempty(msgs)
    text = ok ? "✓ Mesh looks well-sized for ModEM inversion." : join(msgs, "    ")
    return (text = text, ok = ok)
end

const PARAM_HELP = Dict(
    "Core cell N–S (m)" =>
        "Core cell width in the N–S (X) direction, in metres. Use about half the " *
        "station spacing so no cell holds more than one site and structure between " *
        "stations is resolved. Smaller = finer resolution but more cells and slower.",
    "Core cell E–W (m)" =>
        "Core cell width in the E–W (Y) direction, in metres. Use about half the " *
        "station spacing so no cell holds more than one site and structure between " *
        "stations is resolved. Smaller = finer resolution but more cells and slower.",
    "Core cells N–S" =>
        "Number of core cells in the N–S (X) direction. Independent of 'Core cell N–S (m)': " *
        "the core spans count × cell width, centred on the stations.",
    "Core cells E–W" =>
        "Number of core cells in the E–W (Y) direction. Independent of 'Core cell E–W (m)': " *
        "the core spans count × cell width, centred on the stations.",
    "ρ background (Ω·m)" =>
        "Half-space resistivity that fills the starting model. It also sets the skin " *
        "depths that size the vertical mesh and padding. The default is suggested from " *
        "the median off-diagonal apparent resistivity of the data.",
    "Pad factor" =>
        "Geometric growth ratio of the padding cells (~1.4–1.6). A larger ratio expands " *
        "the grid quickly with few cells so the boundaries sit far from the sites.",
    "Pad cells N–S" =>
        "Number of padding cells added on each side in N–S (X). More cells push the " *
        "boundary farther out; it should reach ≥ one skin depth at the longest period " *
        "δ(Tmax) so boundary conditions do not affect the core.",
    "Pad cells E–W" =>
        "Number of padding cells added on each side in E–W (Y). More cells push the " *
        "boundary farther out; it should reach ≥ one skin depth at the longest period " *
        "δ(Tmax) so boundary conditions do not affect the core.",
    "First layer (m)" =>
        "Thickness of the top vertical cell. Set from the shortest-period skin depth " *
        "(≈ δ(Tmin) / 4) so the shallowest structure the data sees is resolved.",
    "Vertical factor" =>
        "Geometric growth of layer thickness with depth (~1.1–1.3). Larger = fewer " *
        "layers but coarser resolution at depth.",
    "Depth × δ(Tmax)" =>
        "Model bottom depth as a multiple of the longest-period skin depth δ(Tmax). " *
        "Use ≥ 1 so the model spans the depth of investigation; ~1.5 is typical.",
    "Covariance" =>
        "Recursive autoregression smoothing value written into the ModEM covariance file. " *
        "Use 0 to turn smoothing off, ~0.1 for weak smoothing, and ~0.3 for moderate smoothing.",
)

function _make_mesh3D_gui(;
    data_file, out_model, out_covariance, topo_file,
    cov_value, cov_apply, n_pad, pad_factor, vertical_factor, depth_mult,
    air_layers, air_factor, air_first_div,
    colormap, resistivity_range, site_color, site_size_full, site_size_core,
    grid_color, grid_linewidth, fig_size, interactive,
    d, sx, sy, Tobs, spacing, ρ_bg_data,
    ρ_bg0, dx_core0, dy_core0, nx_core0, ny_core0, z_first0, topo_ctx)

    use_topo = !isempty(strip(topo_file))

#----- draw the map view of the mesh -----
function draw2d!(ax, m)
    heatmap!(ax, m.y_edges_km, m.x_edges_km, fill(log10(m.ρ_bg), m.ny, m.nx);
             colormap = colormap, colorrange = resistivity_range)
    vlines!(ax, m.y_edges_km; color = grid_color, linewidth = grid_linewidth)
    hlines!(ax, m.x_edges_km; color = grid_color, linewidth = grid_linewidth)
    lines!(ax, corner_points(m.core_y0_km, m.core_y1_km, m.core_x0_km, m.core_x1_km);
           color = :black, linewidth = 2.5)
    scatter!(ax, sy ./ 1000, sx ./ 1000; color = site_color, marker = :circle,
             markersize = site_size)
    return ax
end

fig = Figure(size = fig_size, figure_padding = (40, 14, 14, 14))

status = Observable("Outputs → $(basename(out_model)), $(basename(out_covariance))")
help_obs = Observable("Click  ?  next to a parameter for an explanation.")
sug_obs = Observable("")
sug_color = Observable(RGBf(0.15, 0.5, 0.25))
show_full = Observable(true)
site_size = Observable(site_size_full)

Colorbar(fig[1, 3]; colormap = colormap, colorrange = resistivity_range,
    label = "log₁₀ ρ (Ω·m)", width = 16)

left = GridLayout(fig[1, 1]; valign = :top)

infotext = join([
    @sprintf("datafile  >  %s", basename(data_file)),
    @sprintf("sites     >  %d", length(sx)),
    @sprintf("periods   >  %d   (%.3g – %.3g s)", length(Tobs), minimum(Tobs), maximum(Tobs)),
    @sprintf("spacing   >  %.1f km", spacing.median / 1000),
    @sprintf("rho(data) >  %.0f ohm-m", ρ_bg_data),
    @sprintf("topo      >  %s", use_topo ? basename(topo_file) : "off"),
    @sprintf("cov(mode) >  topography + site elevations"),
    @sprintf("cov(out)  >  %s", basename(out_covariance)),
], "\n")
Label(left[1, 1:3], infotext; font = "Consolas", fontsize = 13, halign = :left,
    justification = :left, color = :gray25, tellwidth = false)

next_row = Ref(1)

function add_param!(name, default; isint = false)
    r = (next_row[] += 1)
    Label(left[r, 1], name; halign = :right, fontsize = 13)
    s = isint ? string(Int(round(default))) : string(default)
    tb = Textbox(left[r, 2]; stored_string = s, validator = Float64, width = 84)
    ib = Button(left[r, 3]; label = "?", fontsize = 13, width = 26)
    on(ib.clicks) do _
        help_obs[] = get(PARAM_HELP, name, "")
    end
    return tb
end

tb_dx   = add_param!("Core cell N–S (m)", dx_core0; isint = true)
tb_dy   = add_param!("Core cell E–W (m)", dy_core0; isint = true)
tb_nxc  = add_param!("Core cells N–S", nx_core0; isint = true)
tb_nyc  = add_param!("Core cells E–W", ny_core0; isint = true)
tb_rho  = add_param!("ρ background (Ω·m)", ρ_bg0; isint = true)
tb_npx  = add_param!("Pad cells N–S", n_pad; isint = true)
tb_npy  = add_param!("Pad cells E–W", n_pad; isint = true)
tb_pf   = add_param!("Pad factor", pad_factor)
tb_zf   = add_param!("First layer (m)", z_first0; isint = true)
tb_zfac = add_param!("Vertical factor", vertical_factor)
tb_dm   = add_param!("Depth × δ(Tmax)", depth_mult)
tb_cov  = add_param!("Covariance", cov_value)

updatebtn = Button(left[(next_row[] += 1), 1:3]; label = "Update mesh", fontsize = 14)

viewgrid = GridLayout(left[(next_row[] += 1), 1:3])
corefullbtn = Button(viewgrid[1, 1]; label = "Show core", fontsize = 14, tellwidth = false)
resetbtn    = Button(viewgrid[1, 2]; label = "Reset zoom", fontsize = 14, tellwidth = false)
colsize!(viewgrid, 1, Relative(0.5)); colsize!(viewgrid, 2, Relative(0.5))
colgap!(viewgrid, 8)

actiongrid = GridLayout(left[(next_row[] += 1), 1:3])
savebtn   = Button(actiongrid[1, 1]; label = "Save model + cov + data", fontsize = 14, tellwidth = false)
exportbtn = Button(actiongrid[1, 2]; label = "Export figure", fontsize = 14, tellwidth = false)
colsize!(actiongrid, 1, Relative(0.5)); colsize!(actiongrid, 2, Relative(0.5))
colgap!(actiongrid, 8)

rowgap!(left, 7)
colgap!(left, 8)
colsize!(left, 2, Fixed(86))
colsize!(left, 3, Fixed(26))

Label(fig[2, 1:3], sug_obs; fontsize = 13, halign = :left, justification = :left,
    color = sug_color, word_wrap = true, tellwidth = false)
Label(fig[3, 1:3], help_obs; fontsize = 12.5, halign = :left, justification = :left,
    color = :gray30, word_wrap = true, tellwidth = false)
Label(fig[4, 1:3], status; fontsize = 11, halign = :left, color = :gray45,
    word_wrap = true, tellwidth = false)

colsize!(fig.layout, 1, Fixed(290))
colsize!(fig.layout, 3, Fixed(70))
rowsize!(fig.layout, 1, Relative(0.82))
for r in 2:4
    rowsize!(fig.layout, r, Auto())
end
rowgap!(fig.layout, 6)

tbget(tb, dflt) = begin
    s = tb.displayed_string[]
    (s === nothing || isempty(strip(s))) && return dflt
    v = tryparse(Float64, s)
    v === nothing ? dflt : v
end

function build_mesh()
    build_mesh_bundle(d, sx, sy, Tobs;
        ρ_bg = tbget(tb_rho, Float64(ρ_bg0)),
        dx_core = tbget(tb_dx, Float64(dx_core0)),
        dy_core = tbget(tb_dy, Float64(dy_core0)),
        nx_core = Int(round(tbget(tb_nxc, Float64(nx_core0)))),
        ny_core = Int(round(tbget(tb_nyc, Float64(ny_core0)))),
        nx_pad = Int(round(tbget(tb_npx, Float64(n_pad)))),
        ny_pad = Int(round(tbget(tb_npy, Float64(n_pad)))),
        pad_factor = tbget(tb_pf, pad_factor),
        z_first = tbget(tb_zf, Float64(z_first0)),
        z_factor = tbget(tb_zfac, vertical_factor),
        depth_mult = tbget(tb_dm, depth_mult),
        cov_value = tbget(tb_cov, cov_value),
        topo_ctx = topo_ctx,
        topo_path = topo_file, n_air = air_layers, air_factor = air_factor,
        air_first_div = air_first_div)
end

mesh = Observable(build_mesh())
current_ax = Ref{Any}(nothing)

function apply_limits!(ax, m)
    if show_full[]
        xlims!(ax, m.full_y0_km, m.full_y1_km)
        ylims!(ax, m.full_x0_km, m.full_x1_km)
    else
        padx = 0.05 * (m.core_y1_km - m.core_y0_km)
        pady = 0.05 * (m.core_x1_km - m.core_x0_km)
        xlims!(ax, m.core_y0_km - padx, m.core_y1_km + padx)
        ylims!(ax, m.core_x0_km - pady, m.core_x1_km + pady)
    end
    return nothing
end

function refresh!(m)
    current_ax[] === nothing || delete!(current_ax[])
    ax = Axis(fig[1, 2]; xlabel = "Y East (km)", ylabel = "X North (km)",
        aspect = DataAspect(), halign = :right,
        xgridvisible = false, ygridvisible = false,
        xlabelfont = :bold, ylabelfont = :bold,
        xticklabelcolor = :gray55, yticklabelcolor = :gray55)
    draw2d!(ax, m)
    current_ax[] = ax
    apply_limits!(ax, m)
    s = suggestions(m)
    topo_info = m.topo_active ? @sprintf("topo %s · datum %.1f m · air %d", m.topo_name, m.datum_elev, m.nz_air) : "topo flat"
    site_info = m.site_air.ok ?
        (m.site_surface.nadjusted == 0 ? "all sites below surface" : @sprintf("all sites below surface (%d auto-adjusted column(s))", m.site_surface.nadjusted)) :
        @sprintf("%d site(s) in air", m.site_air.nbad)
    gridinfo = @sprintf("%d × %d × %d cells · core %d × %d (%.0f × %.0f m) · max sites/cell %d · depth %.0f km · cov %.2f (topography + sites) · %s · %s",
        m.nx, m.ny, m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core, m.maxper, m.depth_km, m.cov_value, topo_info, site_info)
    sug_obs[] = gridinfo * "\n" * s.text
    sug_color[] = s.ok ? RGBf(0.15, 0.5, 0.25) : RGBf(0.75, 0.2, 0.15)
    return nothing
end

on(mesh) do m
    refresh!(m)
end
refresh!(mesh[])

on(updatebtn.clicks) do _
    mesh[] = build_mesh()
end

for tb in (tb_dx, tb_dy, tb_nxc, tb_nyc, tb_rho, tb_npx, tb_npy, tb_pf, tb_zf, tb_zfac, tb_dm, tb_cov)
    on(tb.stored_string) do _
        mesh[] = build_mesh()
    end
end

on(corefullbtn.clicks) do _
    show_full[] = !show_full[]
    corefullbtn.label[] = show_full[] ? "Show core" : "Show full"
    site_size[] = show_full[] ? site_size_full : site_size_core
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(resetbtn.clicks) do _
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(savebtn.clicks) do _
    m = mesh[]
    _save_mesh_outputs(m, out_model, out_covariance; cov_apply = cov_apply)
    air_status = m.site_air.ok ?
        (m.site_surface.nadjusted == 0 ? "all sites below surface" : @sprintf("all sites below surface after %d auto-adjusted column(s)", m.site_surface.nadjusted)) :
        @sprintf("WARNING: %d site(s) in air; see guidance above", m.site_air.nbad)
    status[] = @sprintf("Saved %d×%d×%d model → %s | covariance %.2f (topography + sites) → %s | %s", m.nx, m.ny, m.nz, out_model, m.cov_value, out_covariance, air_status)
    @info status[]
end

on(exportbtn.clicks) do _
    m = mesh[]
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    efig = Figure(size = (1200, 1000))
    eax = Axis(efig[1, 1]; xlabel = "Y East (km)", ylabel = "X North (km)",
        aspect = DataAspect(), halign = :right, xgridvisible = false, ygridvisible = false)
    draw2d!(eax, m)
    apply_limits!(eax, m)
    Colorbar(efig[1, 2]; colormap = colormap, colorrange = resistivity_range,
        label = "log₁₀ ρ (Ω·m)", width = 16)
    fn = joinpath(dirname(out_model), "mesh_map_$ts.png")
    save(fn, efig; px_per_unit = 3)
    status[] = "Saved figure → $fn"
    @info status[]
end

    if interactive
        screen = GLMakie.Screen(fig.scene)
        wait(screen)
    end

    return fig, mesh
end

# Interactive ModEM mesh builder.
# Inputs (optional CLI args): [data_file (ModEM .dat)] [out_model (.rho)];
# defaults to geoenergialoikka/data_new.dat and mesh_start_model.rho next to it.
# Set MODE = "nogui" in the user controls to build+write the default mesh headless.

using MTGeophysics
using GLMakie
using Statistics
using Printf
using Dates

GLMakie.activate!()

#----- user controls (edit here) -----
const MODE            = "gui"                # "gui" = interactive window; "nogui" = headless build + write
const DATA_PATH       = length(ARGS) >= 1 ? ARGS[1] : normpath(@__DIR__, "geoenergialoikka", "data_new.dat")  # input ModEM data file
const OUT_PATH        = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(DATA_PATH), "mesh_start_model.rho")     # output ModEM model (.rho)
const CELL_WIDTH_FRAC = 0.5                  # core cell width = this × median station spacing
const N_PAD           = 12                   # padding cells per side (N–S and E–W)
const PAD_FACTOR      = 1.5                  # geometric padding growth ratio
const FIRST_LAYER_DIV = 4.0                  # first vertical layer = skin_depth(Tmin) / this
const VERTICAL_FACTOR = 1.2                  # geometric vertical growth ratio
const DEPTH_MULT      = 1.5                  # model bottom depth = this × skin_depth(Tmax)
const CMAP            = :Spectral            # resistivity colormap
const CRANGE          = (0.0, 4.0)           # log10 ρ colour range
const SITE_COLOR      = :black               # station marker colour
const SITE_SIZE_FULL  = 6                    # station dot size (px) in full-model view
const SITE_SIZE_CORE  = 5                    # station dot size (px) in core view
const GRID_COLOR      = (:grey, 0.7)         # mesh grid line colour (colour, alpha)
const GRID_WIDTH      = 1.0                  # mesh grid line thickness (≥1 px = crisp)
const FIG_SIZE        = (1500, 1000)         # viewer window size (points)
 
#----- EM skin depth from resistivity and period -----
skin_depth(ρ, T) = 503.0 * sqrt(ρ * T)

#----- find the typical station spacing -----
function nearest_neighbour_spacing(sx, sy)
    n = length(sx)
    n < 2 && return (median = 1000.0, min = 1000.0)
    dmin = fill(Inf, n)
    @inbounds for i in 1:n, j in 1:n
        i == j && continue
        d = hypot(sx[i] - sx[j], sy[i] - sy[j])
        d < dmin[i] && (dmin[i] = d)
    end
    dmin = filter(isfinite, dmin)
    return (median = median(dmin), min = minimum(dmin))
end

#----- count how many sites fall in each core cell -----
function site_occupancy(sx, sy, x_edges, y_edges, ix_core, iy_core)
    counts = Dict{Tuple{Int,Int},Int}()
    nx = length(x_edges) - 1
    ny = length(y_edges) - 1
    for k in eachindex(sx)
        ix = clamp(searchsortedlast(x_edges, sx[k]), 1, nx)
        iy = clamp(searchsortedlast(y_edges, sy[k]), 1, ny)
        counts[(ix, iy)] = get(counts, (ix, iy), 0) + 1
    end
    maxper = isempty(counts) ? 0 : maximum(values(counts))
    incore = count(p -> first(p)[1] in ix_core && first(p)[2] in iy_core, collect(counts))
    ncore = length(ix_core) * length(iy_core)
    occ = ncore == 0 ? 0.0 : incore / ncore
    return (maxper = maxper, occupied = occ)
end

#----- build the mesh from sites, padding and depth settings -----
function design_mesh(sx, sy, T; ρ_bg, dx_core, dy_core, nx_pad, ny_pad,
                     pad_factor, z_first, z_factor, depth_mult,
                     nx_core = nothing, ny_core = nothing)
    dxc = max(1.0, dx_core)
    dyc = max(1.0, dy_core)
    nx_pad = max(0, nx_pad)
    ny_pad = max(0, ny_pad)

    cx_mid = (minimum(sx) + maximum(sx)) / 2
    cy_mid = (minimum(sy) + maximum(sy)) / 2
    nx_core = nx_core === nothing ? max(1, ceil(Int, max(maximum(sx) - minimum(sx), dxc) / dxc)) : max(1, Int(nx_core))
    ny_core = ny_core === nothing ? max(1, ceil(Int, max(maximum(sy) - minimum(sy), dyc) / dyc)) : max(1, Int(ny_core))
    x0 = cx_mid - nx_core * dxc / 2
    y0 = cy_mid - ny_core * dyc / 2

    padx = [dxc * pad_factor^i for i in 1:nx_pad]
    pady = [dyc * pad_factor^i for i in 1:ny_pad]
    dx = vcat(reverse(padx), fill(dxc, nx_core), padx)
    dy = vcat(reverse(pady), fill(dyc, ny_core), pady)

    origin_x = x0 - sum(padx)
    origin_y = y0 - sum(pady)
    origin = [origin_x, origin_y, 0.0]

    Tmin, Tmax = minimum(T), maximum(T)
    δmin = skin_depth(ρ_bg, Tmin)
    δmax = skin_depth(ρ_bg, Tmax)
    z_target = depth_mult * δmax

    dz = Float64[]
    z = 0.0
    t = z_first
    while z < z_target
        push!(dz, t)
        z += t
        t *= z_factor
    end

    x_edges = origin_x .+ vcat(0.0, cumsum(dx))
    y_edges = origin_y .+ vcat(0.0, cumsum(dy))
    z_edges = vcat(0.0, cumsum(dz))

    ix_core = (nx_pad + 1):(nx_pad + nx_core)
    iy_core = (ny_pad + 1):(ny_pad + ny_core)
    diag = site_occupancy(sx, sy, x_edges, y_edges, ix_core, iy_core)

    core_x0 = x_edges[nx_pad + 1]
    core_x1 = x_edges[nx_pad + nx_core + 1]
    core_y0 = y_edges[ny_pad + 1]
    core_y1 = y_edges[ny_pad + ny_core + 1]
    pad_extent = min(sum(padx), sum(pady))

    return (
        dx = dx, dy = dy, dz = dz, origin = origin, ρ_bg = ρ_bg,
        nx = length(dx), ny = length(dy), nz = length(dz),
        nx_core = nx_core, ny_core = ny_core, dx_core = dxc, dy_core = dyc,
        x_edges_km = x_edges ./ 1000, y_edges_km = y_edges ./ 1000,
        z_edges_km = z_edges ./ 1000,
        core_x0_km = core_x0 / 1000, core_x1_km = core_x1 / 1000,
        core_y0_km = core_y0 / 1000, core_y1_km = core_y1 / 1000,
        full_x0_km = x_edges[1] / 1000, full_x1_km = x_edges[end] / 1000,
        full_y0_km = y_edges[1] / 1000, full_y1_km = y_edges[end] / 1000,
        δmin_km = δmin / 1000, δmax_km = δmax / 1000,
        depth_km = z / 1000, pad_extent_km = pad_extent / 1000,
        pad_ok = pad_extent >= δmax,
        maxper = diag.maxper, occupied = diag.occupied,
    )
end

#----- snap a value to the nearest in a range -----
snap(v, r) = collect(r)[argmin(abs.(collect(r) .- v))]

d = load_data_modem(DATA_PATH)                                     # ModEM data file (sites, periods, app. res.)
sx = collect(Float64.(d.x))                                       # site northings X (m)
sy = collect(Float64.(d.y))                                       # site eastings Y (m)
Tobs = collect(Float64.(d.T))                                     # observed periods (s)

spacing = nearest_neighbour_spacing(sx, sy)                       # nearest-neighbour station spacing (m)
ρ_off = filter(x -> isfinite(x) && x > 0, vec(d.ρ[:, [2, 3], :])) # off-diagonal apparent resistivities (Ω·m)
ρ_bg_data = isempty(ρ_off) ? 100.0 : median(ρ_off)               # data-suggested background ρ (Ω·m)
ρ_bg0 = snap(ρ_bg_data, 10:10:20000)                             # background ρ seed for the GUI
dx_core0 = snap(spacing.median * CELL_WIDTH_FRAC, 50:25:10000)    # N–S core cell width seed (m)
dy_core0 = dx_core0                                               # E–W core cell width seed (m)
span_x0 = max(maximum(sx) - minimum(sx), 1.0)                     # N–S station footprint (m)
span_y0 = max(maximum(sy) - minimum(sy), 1.0)                     # E–W station footprint (m)
nx_core0 = max(1, ceil(Int, span_x0 / dx_core0))                  # N–S core cell count seed
ny_core0 = max(1, ceil(Int, span_y0 / dy_core0))                  # E–W core cell count seed
z_first0 = snap(skin_depth(Float64(ρ_bg0), minimum(Tobs)) / FIRST_LAYER_DIV, 5:5:2000)  # first layer seed (m)

@info @sprintf("Loaded %d sites, %d periods (%.3g–%.3g s); median spacing %.0f m; suggested background ρ = %.0f Ω·m (median off-diagonal app. res.)",
               length(sx), length(Tobs), minimum(Tobs), maximum(Tobs), spacing.median, ρ_bg_data)

if MODE == "nogui"
    m = design_mesh(sx, sy, Tobs; ρ_bg = Float64(ρ_bg0),
        dx_core = Float64(dx_core0), dy_core = Float64(dy_core0), nx_pad = N_PAD, ny_pad = N_PAD,
        pad_factor = PAD_FACTOR, z_first = Float64(z_first0), z_factor = VERTICAL_FACTOR, depth_mult = DEPTH_MULT)
    @printf("grid %d×%d×%d (%d cells); core %d×%d @ %.0f×%.0f m; max sites/cell %d; occupied %.0f%%; pad %.0f km/side; depth %.0f km; pad≥δmax %s\n",
        m.nx, m.ny, m.nz, m.nx * m.ny * m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core,
        m.maxper, 100 * m.occupied, m.pad_extent_km, m.depth_km, m.pad_ok ? "yes" : "no")
    A = fill(log10(m.ρ_bg), m.nx, m.ny, m.nz)
    write_ws3d_model(OUT_PATH, m.dx, m.dy, m.dz, A, m.origin; rotation = 0.0, type_str = "LOGE")
    @printf("wrote %s\n", OUT_PATH)
    exit(0)
end

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
)

#----- draw the map view of the mesh -----
function draw2d!(ax, m)
    heatmap!(ax, m.y_edges_km, m.x_edges_km, fill(log10(m.ρ_bg), m.ny, m.nx);
             colormap = CMAP, colorrange = CRANGE)
    vlines!(ax, m.y_edges_km; color = GRID_COLOR, linewidth = GRID_WIDTH)
    hlines!(ax, m.x_edges_km; color = GRID_COLOR, linewidth = GRID_WIDTH)
    lines!(ax, corner_points(m.core_y0_km, m.core_y1_km, m.core_x0_km, m.core_x1_km);
           color = :black, linewidth = 2.5)
    scatter!(ax, sy ./ 1000, sx ./ 1000; color = SITE_COLOR, marker = :circle,
             markersize = site_size)
    return ax
end

fig = Figure(size = FIG_SIZE, figure_padding = (40, 14, 14, 14))  # (left, right, bottom, top)

status = Observable("Output → $(basename(OUT_PATH))")              # save/export status line
help_obs = Observable("Click  ?  next to a parameter for an explanation.")  # parameter help text
sug_obs = Observable("")                                           # grid summary + suggestions text
sug_color = Observable(RGBf(0.15, 0.5, 0.25))                      # suggestion colour (green/red)
show_full = Observable(true)                                       # map shows full grid vs core only
site_size = Observable(SITE_SIZE_FULL)                             # station dot size (px), set by view mode

Colorbar(fig[1, 3]; colormap = CMAP, colorrange = CRANGE,
    label = "log₁₀ ρ (Ω·m)", width = 16)

left = GridLayout(fig[1, 1]; valign = :top)

infotext = join([
    @sprintf("datafile  >  %s", basename(DATA_PATH)),
    @sprintf("sites     >  %d", length(sx)),
    @sprintf("periods   >  %d   (%.3g – %.3g s)", length(Tobs), minimum(Tobs), maximum(Tobs)),
    @sprintf("spacing   >  %.1f km", spacing.median / 1000),
    @sprintf("rho(data) >  %.0f ohm-m", ρ_bg_data),
], "\n")
Label(left[1, 1:3], infotext; font = "Consolas", fontsize = 13, halign = :left,
    justification = :left, color = :gray25, tellwidth = false)

next_row = Ref(1)
#----- add a labelled input box with a help button -----
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

tb_dx   = add_param!("Core cell N–S (m)", dx_core0; isint = true)   # N–S core cell width (m)
tb_dy   = add_param!("Core cell E–W (m)", dy_core0; isint = true)   # E–W core cell width (m)
tb_nxc  = add_param!("Core cells N–S", nx_core0; isint = true)      # N–S core cell count (independent)
tb_nyc  = add_param!("Core cells E–W", ny_core0; isint = true)      # E–W core cell count (independent)
tb_rho  = add_param!("ρ background (Ω·m)", ρ_bg0; isint = true)     # half-space resistivity (Ω·m)
tb_npx  = add_param!("Pad cells N–S", N_PAD; isint = true)          # padding cells per side, N–S
tb_npy  = add_param!("Pad cells E–W", N_PAD; isint = true)          # padding cells per side, E–W
tb_pf   = add_param!("Pad factor", PAD_FACTOR)                      # geometric padding growth
tb_zf   = add_param!("First layer (m)", z_first0; isint = true)     # top vertical cell thickness (m)
tb_zfac = add_param!("Vertical factor", VERTICAL_FACTOR)            # geometric vertical growth
tb_dm   = add_param!("Depth × δ(Tmax)", DEPTH_MULT)                 # bottom depth in skin depths

updatebtn = Button(left[(next_row[] += 1), 1:3]; label = "Update mesh", fontsize = 14)

viewgrid = GridLayout(left[(next_row[] += 1), 1:3])
corefullbtn = Button(viewgrid[1, 1]; label = "Show core", fontsize = 14, tellwidth = false)
resetbtn    = Button(viewgrid[1, 2]; label = "Reset zoom", fontsize = 14, tellwidth = false)
colsize!(viewgrid, 1, Relative(0.5)); colsize!(viewgrid, 2, Relative(0.5))
colgap!(viewgrid, 8)

actiongrid = GridLayout(left[(next_row[] += 1), 1:3])
savebtn   = Button(actiongrid[1, 1]; label = "Save model (.rho)", fontsize = 14, tellwidth = false)
exportbtn = Button(actiongrid[1, 2]; label = "Export figure", fontsize = 14, tellwidth = false)
colsize!(actiongrid, 1, Relative(0.5)); colsize!(actiongrid, 2, Relative(0.5))
colgap!(actiongrid, 8)

rowgap!(left, 7)
colgap!(left, 8)
colsize!(left, 2, Fixed(86))
colsize!(left, 3, Fixed(26))

#----- bottom panel (spans controls + map): mesh quality and parameter info -----
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

#----- read a number from a textbox -----
tbget(tb, d) = begin
    s = tb.displayed_string[]
    (s === nothing || isempty(strip(s))) && return d
    v = tryparse(Float64, s)
    v === nothing ? d : v
end

#----- rebuild the mesh from the current inputs (each field independent) -----
function build_mesh()
    design_mesh(sx, sy, Tobs;
        ρ_bg = tbget(tb_rho, Float64(ρ_bg0)),
        dx_core = tbget(tb_dx, Float64(dx_core0)),
        dy_core = tbget(tb_dy, Float64(dy_core0)),
        nx_core = Int(round(tbget(tb_nxc, Float64(nx_core0)))),
        ny_core = Int(round(tbget(tb_nyc, Float64(ny_core0)))),
        nx_pad = Int(round(tbget(tb_npx, Float64(N_PAD)))),
        ny_pad = Int(round(tbget(tb_npy, Float64(N_PAD)))),
        pad_factor = tbget(tb_pf, PAD_FACTOR),
        z_first = tbget(tb_zf, Float64(z_first0)),
        z_factor = tbget(tb_zfac, VERTICAL_FACTOR),
        depth_mult = tbget(tb_dm, DEPTH_MULT))
end

mesh = Observable(build_mesh())
current_ax = Ref{Any}(nothing)

#----- zoom the axis to core or full extent -----
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

#----- redraw the map and update the messages -----
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
    gridinfo = @sprintf("%d × %d × %d cells · core %d × %d (%.0f × %.0f m) · max sites/cell %d · depth %.0f km",
        m.nx, m.ny, m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core, m.maxper, m.depth_km)
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

for tb in (tb_dx, tb_dy, tb_nxc, tb_nyc, tb_rho, tb_npx, tb_npy, tb_pf, tb_zf, tb_zfac, tb_dm)
    on(tb.stored_string) do _          # pressing Enter in any field applies and re-syncs
        mesh[] = build_mesh()
    end
end

on(corefullbtn.clicks) do _
    show_full[] = !show_full[]
    corefullbtn.label[] = show_full[] ? "Show core" : "Show full"
    site_size[] = show_full[] ? SITE_SIZE_FULL : SITE_SIZE_CORE
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(resetbtn.clicks) do _
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(savebtn.clicks) do _
    m = mesh[]
    A = fill(log10(m.ρ_bg), m.nx, m.ny, m.nz)
    write_ws3d_model(OUT_PATH, m.dx, m.dy, m.dz, A, m.origin; rotation = 0.0, type_str = "LOGE")
    status[] = @sprintf("Saved %d×%d×%d model → %s", m.nx, m.ny, m.nz, OUT_PATH)
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
    Colorbar(efig[1, 2]; colormap = CMAP, colorrange = CRANGE,
        label = "log₁₀ ρ (Ω·m)", width = 16)
    fn = joinpath(dirname(OUT_PATH), "mesh_map_$ts.png")
    save(fn, efig; px_per_unit = 3)
    status[] = "Saved figure → $fn"
    @info status[]
end

screen = display(fig)
if !isinteractive()
    wait(screen)
end

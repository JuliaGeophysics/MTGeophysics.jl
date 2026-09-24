# 2D mesh tool
# Author: @pankajkmishra
# Inversion inputs from a data file, as MakeMesh3D does in 3D: a padded profile mesh sized by the station spacing
# and skin depths, topography and water cut in from topo.dat, and model.start, model.prior, cov.ctrl, mask.ctrl,
# fwd.ctrl, inv.ctrl and vfsa.ctrl written next to each other; mode = :gui opens the GLMakie window
# (MakeMesh2DGUI.jl), :batch does not

using Printf
using Statistics

const _MAKEMESH2D_DEFAULTS = (
    cell_width_frac = 0.5, core_margin_cells = 4, n_pad = 12, pad_factor = 1.3,
    first_layer_div = 5.0, vertical_factor = 1.1, depth_mult = 4.0, background_resistivity = 0.0,
    air_layers = 10, air_thickness = 50_000.0, air_growth = 2.0, air_resistivity = 1e9, dipole_length = 100.0,
    cov_smoothing = 0.3, n_smooth = 1, fixed_below_m = Inf, water_resistivity = 100.0,
)

# median apparent resistivity of the off-diagonal data, the default background
function _makemesh2d_background(data)
    v = filter(x -> isfinite(x) && x > 0, vcat(vec(data.rho_xy), vec(data.rho_yx)))
    isempty(v) ? 100.0 : 10^median(log10.(v))
end

# the mesh, model, mask and data of one parameter set; pure, the gui calls it on every change
function _makemesh2d_build(data, topo, water, p)
    y = sort(data.receivers)
    spacing = length(y) > 1 ? median(diff(y)) : 1000.0
    dy = p.cell_width_frac * spacing
    half = ceil(Int, (maximum(abs.(y)) + p.core_margin_cells * dy) / dy)
    core = fill(dy, 2half)
    pad = dy .* p.pad_factor .^ (1:p.n_pad)
    dys = vcat(reverse(pad), core, pad)
    ρbg = p.background_resistivity > 0 ? p.background_resistivity : _makemesh2d_background(data)
    dz = mt2d_geometric_layers(data.frequencies; background_resistivity = ρbg, first_layer_div = p.first_layer_div,
                               vertical_factor = p.vertical_factor, depth_mult = p.depth_mult)
    model = ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = dys, z_cell_sizes = dz,
                        resistivity = fill(ρbg, length(dz), length(dys)), n_air_cells = 0,
                        origin = [0.0, -sum(dys) / 2, 0.0], rotation = 0.0, format = "LOGE")
    wet = nothing
    if topo !== nothing
        t = Topography2D(model, data, topo; water, water_resistivity = p.water_resistivity)
        model, data, wet = t.model, t.data, t.mask .== MT2D_MASK_WATER
    end
    mask = Mask2D(model; water = wet, fixed_below_m = p.fixed_below_m)
    fwd = FwdCtrl2D(mode = :TETM, air_layers = p.air_layers, air_thickness = p.air_thickness, air_growth = p.air_growth,
                    air_resistivity = p.air_resistivity, dipole_length = p.dipole_length)
    mesh, _ = Mesh2DFromInputs(model, data, fwd; warn = false)

    # checks, as MakeMesh3D advises
    δmax = mt2d_skin_depth(ρbg, minimum(data.frequencies))
    δmin = mt2d_skin_depth(ρbg, maximum(data.frequencies))
    percell = maximum(values(Dict(c => count(==(c), mt2d_receiver_columns(mesh)) for c in mt2d_receiver_columns(mesh))))
    reach = sum(pad)
    notes = String[]
    percell > 1 && push!(notes, "more than one station in a cell: lower the cell width")
    reach < δmax && push!(notes, @sprintf("padding reaches %.1f km, less than δ(f_min) = %.1f km: add padding cells or growth", reach / 1000, δmax / 1000))
    dz[1] > δmin / 3 && push!(notes, @sprintf("first layer %.0f m is thicker than δ(f_max)/3 = %.0f m", dz[1], δmin / 3))
    sum(dz) < δmax && push!(notes, "the model is shallower than δ(f_min): raise the depth multiplier")
    off = filter(o -> abs(o.offset) > o.tolerance, mt2d_station_offsets(mesh, data))
    isempty(off) || push!(notes, "$(length(off)) station(s) snap by more than half a cell")
    (; model, mask, data, fwd, mesh, ρbg, dy, δmax, δmin, notes,
       summary = @sprintf("%d × %d cells, core %d × %.0f m, padding %.1f km, depth %.1f km, first layer %.1f m, ρ %.0f Ω·m%s",
                          length(dz), length(dys), 2half, dy, reach / 1000, sum(dz) / 1000, dz[1], ρbg,
                          topo === nothing ? "" : @sprintf(", %d air and %d water cells", count(==(0), mask), count(==(9), mask))))
end

function _makemesh2d_write(b, out_dir, p, ctrls, topo_given)
    mkpath(out_dir)
    path(f) = joinpath(out_dir, f)
    m = b.model
    paths = (
        start_model_path = WriteModel2D(path("model.start"), m.y_cell_sizes, m.z_cell_sizes, m.resistivity),
        prior_model_path = WriteModel2D(path("model.prior"), m.y_cell_sizes, m.z_cell_sizes, m.resistivity),
        cov_path = WriteCov2D(path("cov.ctrl"), Cov2D(sy = fill(p.cov_smoothing, length(m.z_cell_sizes)), sz = p.cov_smoothing,
                                                      n_smooth = p.n_smooth, mask = b.mask)),
        mask_path = WriteMask2D(path("mask.ctrl"), b.mask),
        fwd_path = WriteFwdCtrl2D(path("fwd.ctrl"), b.fwd),
        inv_path = (cp(ctrls.inv, path("inv.ctrl"); force = true); path("inv.ctrl")),
        vfsa_path = (cp(ctrls.vfsa, path("vfsa.ctrl"); force = true); path("vfsa.ctrl")),
        data_path = topo_given ? write_data2d(path("data.dat"), b.data) : "",
    )
    plot_mt2d_mesh(b.mesh; output_path = path("Mesh.png"), background_resistivity = b.ρbg)
    paths
end

"""
    MakeMesh2D(data_path; out_dir=dirname(data_path), topo_path="", water=[], mode=:batch,
               inv_ctrl=examples/ctrl/2D/InvCtrl.GN, vfsa_ctrl=examples/ctrl/2D/InvCtrl.VFSA,
               kwargs...) -> (; paths, summary, notes)

Inversion inputs for a 2D data file, as `MakeMesh3D` does in 3D:
- lateral: a uniform core of `cell_width_frac` × the median station spacing over the
  stations plus `core_margin_cells`, centred on y = 0, then `n_pad` cells growing by
  `pad_factor` on each side;
- vertical: `mt2d_geometric_layers` (first layer δ(f_max)/`first_layer_div`, growth
  `vertical_factor`, down to `depth_mult` δ(f_min)) for `background_resistivity`
  (0 = the median apparent resistivity of the data);
- with `topo_path` (`topo.dat`, WGS84 lat lon elevation) and `water` lakes
  `(y_range, level)`: the topography cut in by `Topography2D`, and `data.dat` rewritten
  with Z = depth below the model top;
- `cov.ctrl` (GN, NLCG) from `Mask2D` (air 0, water 9, `fixed_below_m`), smoothing
  `cov_smoothing`, and the same mask as `mask.ctrl` (VFSA);
- `fwd.ctrl` with the air (`air_layers`, `air_thickness`, `air_growth`, `air_resistivity`)
  and `dipole_length`; `inv.ctrl` copied from `inv_ctrl`, `vfsa.ctrl` from `vfsa_ctrl`.

Writes `model.start`, `model.prior` (the background), `cov.ctrl`, `mask.ctrl`, `fwd.ctrl`,
`inv.ctrl`, `vfsa.ctrl`, `data.dat` (with topography) and `Mesh.png` into `out_dir`, and prints the mesh summary
and advice. `mode = :gui` opens an interactive window (needs GLMakie and a display).
"""
function MakeMesh2D(data_path::AbstractString; out_dir::AbstractString = dirname(abspath(data_path)),
                    topo_path::AbstractString = "", water = NamedTuple[], mode::Symbol = :batch,
                    inv_ctrl::AbstractString = joinpath(dirname(@__DIR__), "examples", "ctrl", "2D", "InvCtrl.GN"),
                    vfsa_ctrl::AbstractString = joinpath(dirname(@__DIR__), "examples", "ctrl", "2D", "InvCtrl.VFSA"),
                    kwargs...)
    p = merge(_MAKEMESH2D_DEFAULTS, values(kwargs))
    unknown = setdiff(keys(p), keys(_MAKEMESH2D_DEFAULTS))
    isempty(unknown) || throw(ArgumentError("unknown MakeMesh2D settings $unknown"))
    data = load_data2d(data_path)
    topo = isempty(topo_path) ? nothing : ReadTopo2D(topo_path)
    ReadInvCtrl2D(inv_ctrl); ReadVFSACtrl2D(vfsa_ctrl)
    ctrls = (inv = inv_ctrl, vfsa = vfsa_ctrl)
    if mode == :gui
        isdefined(@__MODULE__, :GLMakie) || error("the MakeMesh2D window needs GLMakie and a display; use mode = :batch")
        return _make_mesh2D_gui(data, topo, water, p, out_dir, ctrls)
    end
    mode == :batch || throw(ArgumentError("mode must be :batch or :gui"))
    b = _makemesh2d_build(data, topo, water, p)
    paths = _makemesh2d_write(b, out_dir, p, ctrls, topo !== nothing)
    println(b.summary)
    foreach(n -> println("  • ", n), b.notes)
    (; paths, summary = b.summary, notes = b.notes, mesh = b.mesh)
end

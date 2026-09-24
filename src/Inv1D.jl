# 1D MT inversion
# Author: @pankajkmishra
# Every site of a data file inverted on its own, on the layering MakeMesh1D gives it, from its background
# halfspace: Gauss-Newton through the shared Inv2D driver or VFSA2D. One small control file, the rest fixed

using CairoMakie
using Printf

#---------- control ----------

"""
    InvCtrl1D

1D inversion control, read from `InvCtrl.GN` or `InvCtrl.VFSA` (`examples/ctrl/1D`).
- `algorithm`: `:gn` or `:vfsa`
- `mode`: impedance fitted, `:XY`, `:YX`, `:XYYX` or `:DET`
- `target_rms`, `max_iter`, `log_bounds`: stopping controls and the log10 ρ box
- `lambda`: GN regularization weight (vertical smoothing, small smallness), fixed
- `chains`, `control_points`: VFSA chains and RBF control points per chain

Fixed inside: GN damping 0.01, smallness 0.01, max log10 step 0.5; VFSA temperature 1,
cooling ratio 0.001, step scale 0.11, RBF width 2.5 layers.
"""
Base.@kwdef struct InvCtrl1D
    algorithm::Symbol
    target_rms::Float64
    max_iter::Int
    mode::Symbol = :XYYX
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    lambda::Float64 = 1.0
    chains::Int = 10
    control_points::Int = 25
end

_ctrl_algorithm1d(s) = (a = Symbol(lowercase(strip(s))); a in (:gn, :vfsa) ? a :
                        throw(ArgumentError("1D algorithm must be GN or VFSA")))

const _INV1D_CTRL_SPEC = (
    "Algorithm"                          => (:algorithm, _ctrl_algorithm1d, true),
    "Mode"                               => (:mode, _ctrl_mode1d, false),
    "Exit search when rms is less than"  => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"       => (:max_iter, _ctrl_int, true),
    "Log10 resistivity bounds"           => (:log_bounds, _ctrl_parse_pair, false),
    "Initial damping factor lambda"      => (:lambda, _ctrl_float, false),
    "Number of chains"                   => (:chains, _ctrl_int, false),
    "Control points"                     => (:control_points, _ctrl_int, false),
)

function _validate_ctrl(c::InvCtrl1D)
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c.log_bounds[1] < c.log_bounds[2] || throw(ArgumentError("log10 bounds must be increasing"))
    c.lambda >= 0 || throw(ArgumentError("lambda must be nonnegative"))
    c.chains >= 1 && c.control_points >= 1 || throw(ArgumentError("chains and control points must be at least 1"))
    c
end

"""
    ReadInvCtrl1D(path) -> InvCtrl1D

Read a 1D inversion control file. GN files may not carry the VFSA keys, nor VFSA files lambda.
"""
function ReadInvCtrl1D(path::AbstractString)
    v = _read_ctrl(path, _INV1D_CTRL_SPEC, "1D inversion control")
    other = v[:algorithm] == :gn ? (:chains, :control_points) : (:lambda,)
    foreign = filter(k -> haskey(v, k), other)
    isempty(foreign) || error("$path: $(join(foreign, ", ")) do not apply to $(uppercase(string(v[:algorithm])))")
    _validate_ctrl(InvCtrl1D(; v...))
end

"""
    WriteInvCtrl1D(path, ctrl::InvCtrl1D) -> path

Write a 1D inversion control file with its algorithm's keys.
"""
function WriteInvCtrl1D(path::AbstractString, c::InvCtrl1D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g(x) = @sprintf("%.6g", x)
    rows = [("Algorithm", uppercase(string(c.algorithm))), ("Mode", string(c.mode)),
            ("Exit search when rms is less than", g(c.target_rms)), ("Maximum number of iterations", string(c.max_iter)),
            ("Log10 resistivity bounds", "$(g(c.log_bounds[1])) $(g(c.log_bounds[2]))")]
    c.algorithm == :gn ? push!(rows, ("Initial damping factor lambda", g(c.lambda))) :
        append!(rows, [("Number of chains", string(c.chains)), ("Control points", string(c.control_points))])
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

#---------- inversion ----------

"""
    Invert1D(data_path, inv_path, meshes=MakeMesh1D(data_path); sites=nothing, run_dir=nothing)

1D inversion of every site of the data file (or those named in `sites`) on its own, on its
`MakeMesh1D` layering, from its background halfspace (also the GN reference model). The
control (`InvCtrl1D`) picks GN or VFSA and the impedance. Writes into `run_dir` (default
`run_YYYYmmdd_HHMMSS/` next to the data) one folder per site with `model.start`,
`model.rho`, `data.pred` and `History.csv` (GN) or `vfsa/` (VFSA, `model.rho` = ensemble
mean), plus `data.pred` of all sites, `Summary.txt` and the inputs. Returns the run
directory, the per-site results and the overall rms.
"""
function Invert1D(data_path::AbstractString, inv_path::AbstractString, meshes = nothing;
                  sites = nothing, run_dir::Union{Nothing, AbstractString} = nothing)
    observed = load_data2d(data_path)
    ctrl = ReadInvCtrl1D(inv_path)
    meshes = something(meshes, MakeMesh1D(observed; mode = ctrl.mode))
    mode = _mt1d_driver_mode(ctrl.mode)
    picked = sites === nothing ? collect(eachindex(observed.site_names)) :
             [something(findfirst(==(String(s)), observed.site_names), 0) for s in sites]
    all(>(0), picked) || error("unknown sites $(sites[picked .== 0])")
    lookup = Dict(m.site => m for m in meshes)
    missing_mesh = setdiff(observed.site_names[picked], keys(lookup))
    isempty(missing_mesh) || error("no mesh for sites $missing_mesh")
    lo, hi = ctrl.log_bounds
    dir = _inv2d_open_run(run_dir, data_path, (data_path, inv_path))

    results = map(picked) do i
        site = mt1d_site_data(observed, i; mode = ctrl.mode)
        name = site.site_names[1]
        m = lookup[name]
        mesh = Mesh1D(m.thicknesses, observed.frequencies)
        ρ0 = fill(10^clamp(log10(m.background), lo, hi), length(m.thicknesses), 1)
        sdir = mkpath(joinpath(dir, name))
        WriteModel2D(joinpath(sdir, "model.start"), [1.0], m.thicknesses, ρ0)
        history, vfsa = nothing, nothing
        if ctrl.algorithm == :gn
            options = Inv2DOptions(mode = mode, max_iter = ctrl.max_iter, beta = ctrl.lambda, smallness = 0.01,
                                   smooth_z = 1.0, log_bounds = ctrl.log_bounds, max_step = 0.5, target_rms = ctrl.target_rms)
            r = Invert2D(mesh, ρ0, site; algorithm = GaussNewton2DConfig(damping = 1e-2), options)
            final, response, history, reason, converged = r.resistivity, r.response, r.history, r.reason, r.converged
            _inv2d_write_history(joinpath(sdir, "History.csv"), history)
        else
            config = VFSA2DConfig(n_chains = ctrl.chains, n_ctrl = min(ctrl.control_points, length(m.thicknesses)),
                                  max_iter = ctrl.max_iter, log_bounds = ctrl.log_bounds, target_rms = ctrl.target_rms,
                                  rbf_sigma_scale_z = 2.5, mode = mode)
            vfsa = VFSA2D(mesh, ρ0, site; config, run_dir = sdir)
            final, response = vfsa.resistivity, vfsa.response
            converged = vfsa.rms <= ctrl.target_rms
            reason = converged ? :target_rms : :max_iter
        end
        predicted = _inv2d_predicted(response, site, mode)
        WriteModel2D(joinpath(sdir, "model.rho"), mesh, final)
        write_data2d(joinpath(sdir, "data.pred"), predicted)
        fit = chi2_rms2d(site, predicted; components = _inv2d_components(mode))
        (; site = name, index = i, run_dir = sdir, mesh, observed = site, predicted, start = ρ0, prior = ρ0, final,
           active = trues(size(ρ0)), history, vfsa, fit, rms = fit.rms, reason, converged, ctrl, water = falses(size(ρ0)))
    end

    # all sites in one file, at their own positions
    cat2(k) = reduce(hcat, [getproperty(r.predicted, k) for r in results])
    predicted = data_from_response2d(
        MT2DResponse(frequencies = observed.frequencies, periods = observed.periods, receivers = observed.receivers[picked],
                     rho_xy = cat2(:rho_xy), phase_xy = cat2(:phase_xy), z_xy = cat2(:z_xy),
                     rho_yx = cat2(:rho_yx), phase_yx = cat2(:phase_yx), z_yx = cat2(:z_yx));
        z_xy_error = cat2(:z_xy_error), z_yx_error = cat2(:z_yx_error), site_names = observed.site_names[picked],
        x_positions = observed.x_positions[picked], z_positions = observed.z_positions[picked],
        latitudes = isempty(observed.latitudes) ? Float64[] : observed.latitudes[picked],
        longitudes = isempty(observed.longitudes) ? Float64[] : observed.longitudes[picked], origin = observed.origin)
    mode == :TE && (predicted.z_yx .= NaN)
    mode == :TM && (predicted.z_xy .= NaN)
    write_data2d(joinpath(dir, "data.pred"), predicted)
    chi2, ndata = sum(r.fit.chi2 for r in results), sum(r.fit.count for r in results)
    rms = sqrt(chi2 / ndata)
    open(joinpath(dir, "Summary.txt"), "w") do io
        println(io, "Algorithm: ", uppercase(string(ctrl.algorithm)))
        println(io, "Impedance: ", ctrl.mode)
        @printf(io, "RMS: %.6f\n", rms)
        println(io, "Real data count: ", ndata)
        for r in results
            @printf(io, "  %-12s RMS %.6f  %s  %d layers, start %.1f ohm m%s\n", r.site, r.rms, r.reason,
                    length(r.mesh.z_cell_sizes), r.start[1],
                    r.history === nothing ? "" : "  $(length(r.history) - 1) iterations")
        end
    end
    (; run_dir = dir, algorithm = ctrl.algorithm, ctrl, sites = results, predicted, rms)
end

#---------- plots ----------

"""
    plot_mt1d_model(z_cell_sizes, models; output_path, band=nothing, maximum_depth=nothing) -> path

Resistivity-depth steps of one or more 1D models on log axes. `models` holds
`(resistivity, label, color, linestyle)` tuples; `band = (low, high)` shades a range,
e.g. the VFSA 5-95% range. The last layer is drawn down to `maximum_depth` (default
twice its top).
"""
function plot_mt1d_model(z_cell_sizes::AbstractVector{<:Real}, models; output_path::AbstractString,
                         band = nothing, maximum_depth = nothing)
    CairoMakie.activate!()
    tops = vcat(0.0, cumsum(z_cell_sizes[1:end-1]))
    bottom = something(maximum_depth, 2 * max(tops[end], z_cell_sizes[1]))
    ztop = z_cell_sizes[1] / 2
    edges(z) = vcat(max(ztop, 1e-3), repeat(max.(z[2:end], ztop), inner = 2), bottom)
    steps(ρ) = repeat(vec(ρ), inner = 2)
    figure = Figure(size = (600, 750))
    axis = _mt_axis(figure[1, 1]; xlabel = "Resistivity (Ω·m)", ylabel = "Depth (m)", xscale = log10, yscale = log10,
                    yreversed = true)
    if band !== nothing
        z = edges(tops)
        band!(axis, Point2f.(steps(band[1]), z), Point2f.(steps(band[2]), z); color = (:steelblue, 0.25), label = "5-95%")
    end
    for (ρ, label, color, style) in models
        lines!(axis, steps(ρ), edges(tops); color, linestyle = style, linewidth = 2, label)
    end
    ylims!(axis, bottom, ztop)
    axislegend(axis, position = :lb, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

"""
    PlotModel1D(model_path; output_path, true_model_path=nothing, maximum_depth=nothing) -> path

Plot a one-column model file, with the true model when given.
"""
function PlotModel1D(model_path::AbstractString; output_path::AbstractString,
                     true_model_path::Union{Nothing, AbstractString} = nothing, maximum_depth = nothing)
    m = ReadModel2D(model_path)
    models = Any[(m.resistivity[:, 1], "model", :black, :solid)]
    if true_model_path !== nothing
        t = ReadModel2D(true_model_path)
        push!(models, (_mt1d_resample(t, m.z_cell_sizes), "true", :firebrick, :solid))
    end
    plot_mt1d_model(m.z_cell_sizes, models; output_path, maximum_depth)
end

# a 1D model's resistivity at the layer centres of another layering
function _mt1d_resample(model::ModelFile2D, dz::AbstractVector{<:Real})
    tops = cumsum(model.z_cell_sizes)
    centres = cumsum(dz) .- dz ./ 2
    [model.resistivity[min(searchsortedfirst(tops, z), length(tops)), 1] for z in centres]
end

"""
    PlotInversion1D(run; true_model_path=nothing, maximum_depth=nothing) -> paths

Standard plots of a run returned by `Invert1D`, into each site's `plots/`: the start,
final and (when given) true models, with the VFSA median and 5-95% range; the data fit;
the convergence.
"""
function PlotInversion1D(run; true_model_path::Union{Nothing, AbstractString} = nothing, maximum_depth = nothing)
    paths = String[]
    truth = true_model_path === nothing ? nothing : ReadModel2D(true_model_path)
    for r in run.sites
        path(n) = joinpath(r.run_dir, "plots", n)
        dz = r.mesh.z_cell_sizes
        models = Any[(r.start[:, 1], "start", :gray, :dash), (r.final[:, 1], r.vfsa === nothing ? "final" : "mean", :black, :solid)]
        band = nothing
        if r.vfsa !== nothing
            e = r.vfsa.ensemble
            push!(models, (10 .^ e.median[:, 1], "median", :steelblue, :solid))
            band = (10 .^ e.p05[:, 1], 10 .^ e.p95[:, 1])
        end
        truth === nothing || push!(models, (_mt1d_resample(truth, dz), "true", :firebrick, :solid))
        push!(paths, plot_mt1d_model(dz, models; output_path = path("Model.png"), band, maximum_depth))
        push!(paths, plot_mt2d_data_fit(r.observed, r.predicted; output_path = path("DataFit.png"),
                                        names = _mt1d_plot_names(r.ctrl.mode)))
        r.history === nothing || push!(paths, plot_inv2d_convergence(r.history; output_path = path("Convergence.png"),
                                                                      target_rms = r.ctrl.target_rms))
        r.vfsa === nothing || push!(paths, plot_vfsa2d_convergence([c.history for c in r.vfsa.chains];
                                                                    output_path = path("ConvergenceVFSA.png"),
                                                                    target_rms = r.ctrl.target_rms))
    end
    paths
end

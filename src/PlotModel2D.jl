# 2D MT model plots
# Author: @pankajkmishra
# Resistivity sections, mesh layout with skin-depth core, and inversion convergence

using CairoMakie

"""
    plot_mt2d_model(mesh, resistivity; output_path, show_air=false, show_grid=false, show_padding=true,
                    maximum_depth_km=Inf, resistivity_log10_range=(0.0, 4.0), annotation=nothing) -> path

Resistivity section with the stations on top. `annotation` is a short label drawn
inside the section, top left, e.g. an iteration number; there is no title.
"""
function plot_mt2d_model(
    mesh::MT2DMesh,
    resistivity::AbstractMatrix{<:Real};
    output_path::AbstractString,
    show_air::Bool = false,
    show_grid::Bool = false,
    show_padding::Bool = true,
    maximum_depth_km::Real = Inf,
    resistivity_log10_range::Tuple{Float64, Float64} = (0.0, 4.0),
    annotation::Union{Nothing, AbstractString} = nothing,
)
    CairoMakie.activate!()

    core_y = _core_range(mesh.y_cell_sizes)
    column_range = show_padding ? (1:size(resistivity, 2)) : core_y
    y_edges = mesh.y_nodes[first(column_range):(last(column_range) + 1)] ./ 1000
    row_range = show_air ? (1:size(resistivity, 1)) : ((mesh.n_air_cells + 1):size(resistivity, 1))
    z_edges = show_air ? (mesh.z_nodes ./ 1000) : (mesh.z_nodes[(mesh.n_air_cells + 1):end] ./ 1000)
    rho_plot = log10.(resistivity[row_range, column_range])
    show_air && (rho_plot[1:mesh.n_air_cells, :] .= NaN)      # air drawn as nan_color

    figure = Figure(size = (1100, 650))
    axis = _mt_axis(figure[1, 1]; xlabel = "Offset (km)", ylabel = "Depth (km)", yreversed = true)
    heatmap = heatmap!(axis, y_edges, z_edges, rho_plot', colormap = :Spectral, colorrange = resistivity_log10_range,
                       nan_color = :aliceblue)
    Colorbar(figure[1, 2], heatmap, label = "log10 ρ (Ω·m)", labelfont = :regular)
    xlims!(axis, minimum(y_edges), maximum(y_edges))

    scatter!(
        axis,
        mesh.receiver_positions ./ 1000,
        fill(0.0, length(mesh.receiver_positions));
        marker = :dtriangle,
        markersize = 12,
        color = :black,
    )

    if show_grid
        for edge in y_edges
            vlines!(axis, [edge], color = (:black, 0.15), linewidth = 1)
        end
        for edge in z_edges
            hlines!(axis, [edge], color = (:black, 0.15), linewidth = 1)
        end
    end

    if !show_air
        depth_limit_km = isfinite(Float64(maximum_depth_km)) ? min(Float64(maximum_depth_km), maximum(z_edges)) : maximum(z_edges)
        ylims!(axis, depth_limit_km, 0.0)
    end
    annotation === nothing || text!(axis, 0.01, 0.98; text = annotation, space = :relative,
                                    align = (:left, :top), fontsize = 14)
    _mt_save(output_path, figure)
end

"""
    PlotModel2D(model_path; output_path, show_grid=false, show_padding=true, maximum_depth_km=Inf,
                resistivity_log10_range=(0.0, 4.0)) -> path

Plot a model file, ModEM layout or the older layout with air, earth cells only.
"""
function PlotModel2D(
    model_path::AbstractString;
    output_path::AbstractString,
    show_grid::Bool = false,
    show_padding::Bool = true,
    maximum_depth_km::Real = Inf,
    resistivity_log10_range::Tuple{Float64, Float64} = (0.0, 4.0),
)
    model = ReadModel2D(model_path)
    mesh = _mt2d_earth_mesh(model)
    plot_mt2d_model(
        mesh,
        model.resistivity;
        output_path = output_path,
        show_grid = show_grid,
        show_padding = show_padding,
        maximum_depth_km = maximum_depth_km,
        resistivity_log10_range = resistivity_log10_range,
    )
end

#---------- mesh ----------

"""
    plot_mt2d_mesh(mesh; output_path, region=:full, background_resistivity=100.0) -> path

Cell edges with the air shaded, the core outlined and the skin depths of the lowest and
highest frequency in `background_resistivity` marked. The core is found as in 3D:
laterally `core_indices` on the cell widths, and down to the layer boundary nearest the
skin depth of the lowest frequency, as `z_indices_for_max_depth`; `region = :core` shows
it only.
"""
function plot_mt2d_mesh(
    mesh::MT2DMesh;
    output_path::AbstractString,
    region::Symbol = :full,
    background_resistivity::Real = 100.0,
)
    region in (:full, :core) || error("region must be :full or :core")
    CairoMakie.activate!()

    #---------- core extent ----------
    na = mesh.n_air_cells
    core_y = _core_range(mesh.y_cell_sizes)
    δ_max = mt2d_skin_depth(background_resistivity, minimum(mesh.frequencies)) / 1000
    δ_min = mt2d_skin_depth(background_resistivity, maximum(mesh.frequencies)) / 1000
    z_bottom = na + last(_depth_range(mesh.z_cell_sizes[na+1:end], 1000δ_max)) + 1
    y_core = (mesh.y_nodes[first(core_y)], mesh.y_nodes[last(core_y)+1]) ./ 1000
    z_core = mesh.z_nodes[z_bottom] / 1000

    #---------- visible edges ----------
    ys, zs = mesh.y_nodes ./ 1000, mesh.z_nodes ./ 1000
    if region == :core
        ys = ys[first(core_y):last(core_y)+1]
        zs = zs[na+1:z_bottom]
    end

    figure = Figure(size = (1100, 650))
    axis = _mt_axis(figure[1, 1]; xlabel = "Offset (km)", ylabel = "Depth (km)", yreversed = true)

    region == :full && poly!(axis, Rect(ys[1], zs[1], ys[end] - ys[1], -zs[1]), color = :aliceblue)
    vlines!(axis, ys, color = (:black, 0.35), linewidth = 0.6)
    hlines!(axis, zs, color = (:black, 0.35), linewidth = 0.6)
    region == :full && lines!(axis, [y_core[1], y_core[2], y_core[2], y_core[1], y_core[1]],
                              [0, 0, z_core, z_core, 0], color = :firebrick, linewidth = 2, label = "core")
    hlines!(axis, [δ_max], color = :darkorange, linestyle = :dash, linewidth = 2,
            label = @sprintf("δ(f_min) = %.1f km", δ_max))
    hlines!(axis, [δ_min], color = :teal, linestyle = :dot, linewidth = 2,
            label = @sprintf("δ(f_max) = %.2f km", δ_min))
    hlines!(axis, [0.0], color = :black, linewidth = 1.5)
    scatter!(axis, mesh.receiver_positions ./ 1000, zeros(length(mesh.receiver_positions));
             marker = :dtriangle, markersize = 12, color = :black)

    xlims!(axis, ys[1], ys[end])
    ylims!(axis, max(zs[end], 1.02δ_max), min(zs[1], 0.0))       # the core may stop just above δ(f_min)
    Legend(figure[1, 2], axis, @sprintf("ρ_bg = %g Ω·m", background_resistivity), framevisible = false,
           titlefont = :regular, labelfont = :regular)
    _mt_save(output_path, figure)
end

#---------- inversion convergence ----------

"""
    plot_inv2d_convergence(history; output_path, target_rms=0.0) -> path

Rms (with the target line when `target_rms > 0`) and the objective terms per iteration.
"""
function plot_inv2d_convergence(history::AbstractVector; output_path::AbstractString, target_rms::Real = 0.0)
    CairoMakie.activate!()

    it = [h.iteration for h in history]
    figure = Figure(size = (1100, 450))
    ax1 = _mt_axis(figure[1, 1]; xlabel = "Iteration", ylabel = "RMS")
    scatterlines!(ax1, it, [h.rms for h in history], color = :navy)
    target_rms > 0 && hlines!(ax1, [target_rms], color = :gray, linestyle = :dash)

    ax2 = _mt_axis(figure[1, 2]; xlabel = "Iteration", ylabel = "Objective terms", yscale = log10)
    scatterlines!(ax2, it, [h.objective for h in history], label = "objective")
    scatterlines!(ax2, it, [h.chi2 / 2 for h in history], label = "χ²/2")
    # regularization is exactly zero at the start model, skip it on a log axis
    reg = [(h.iteration, h.regularization) for h in history if h.regularization > 0]
    isempty(reg) || scatterlines!(ax2, first.(reg), last.(reg), label = "‖R(m - m_ref)‖²/2")
    axislegend(ax2, position = :rt, framevisible = false, labelfont = :regular)
    _mt_save(output_path, figure)
end

#---------- inversion run ----------

# earth-only mesh for plotting a model file on its own grid
_mt2d_earth_mesh(model::ModelFile2D; receivers = Float64[]) = MT2DMesh(
    y_nodes = model.origin[2] .+ vcat(0.0, cumsum(model.y_cell_sizes)), z_nodes = vcat(0.0, cumsum(model.z_cell_sizes)),
    y_cell_sizes = model.y_cell_sizes, z_cell_sizes = model.z_cell_sizes, receiver_positions = Float64.(receivers),
    frequencies = [1.0], n_air_cells = 0)

"""
    PlotInversion2D(run; true_model_path=nothing, maximum_depth_km=10.0,
                    resistivity_log10_range=(0.0, 4.0), background_resistivity=100.0) -> paths

Standard plots of a run returned by the six-file `Invert2D`, written to `run.run_dir/plots`:
the mesh, the start, final and (when given) true models, the data fit, and the
convergence for GN and NLCG.
"""
function PlotInversion2D(run; true_model_path::Union{Nothing, AbstractString} = nothing,
                         maximum_depth_km::Real = 10.0, resistivity_log10_range = (0.0, 4.0),
                         background_resistivity::Real = 100.0)
    dir = joinpath(run.run_dir, "plots")
    path(name) = joinpath(dir, name)
    core = (show_padding = false, maximum_depth_km, resistivity_log10_range)
    paths = String[
        plot_mt2d_mesh(run.mesh; output_path = path("Mesh.png"), region = :full, background_resistivity),
        plot_mt2d_mesh(run.mesh; output_path = path("MeshCore.png"), region = :core, background_resistivity),
        plot_mt2d_model(run.mesh, run.start; output_path = path("ModelStart.png"), core...),
        plot_mt2d_model(run.mesh, run.final; output_path = path("ModelFinal.png"), core...),
        plot_mt2d_model(run.mesh, run.final; output_path = path("ModelFinalFull.png"), resistivity_log10_range),
        plot_mt2d_data_fit(run.observed, run.predicted; output_path = path("DataFit.png")),
    ]
    run.history === nothing || push!(paths, plot_inv2d_convergence(run.history; output_path = path("Convergence.png"),
                                                                    target_rms = run.ctrl.target_rms))
    if true_model_path !== nothing
        truth = ReadModel2D(true_model_path)
        push!(paths, plot_mt2d_model(_mt2d_earth_mesh(truth; receivers = run.mesh.receiver_positions), truth.resistivity;
                                     output_path = path("ModelTrue.png"), core...))
    end
    paths
end

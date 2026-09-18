# 2D MT model plots
# Author: @pankajkmishra
# Resistivity sections, mesh layout with skin-depth core, and inversion convergence

using CairoMakie

"""
    _plot_core_indices_2d(cell_sizes; tol=0.20)

Inputs:
- Cell sizes and a core-detection tolerance.

Output:
- `UnitRange{Int}`: Detected core-cell range.

Description:
- Detects the central uniform-cell part of the 2D mesh for plotting.
"""
function _plot_core_indices_2d(cell_sizes::AbstractVector{<:Real}; tol::Real = 0.20)
    minimum_size = minimum(Float64.(cell_sizes))
    threshold = minimum_size * (1 + Float64(tol))
    indices = findall(size -> Float64(size) <= threshold + 1e-9, cell_sizes)
    isempty(indices) ? (1:length(cell_sizes)) : (first(indices):last(indices))
end

"""
    plot_mt2d_model(mesh, resistivity; output_path, show_air=false, show_grid=false, show_padding=true,
                    maximum_depth_km=Inf, resistivity_log10_range=(0.0, 4.0), title=nothing)

Inputs:
- 2D mesh, resistivity model, output image path, and plotting controls.

Output:
- `String`: Path to the written plot.

Description:
- Writes the standard 2D resistivity-model plot.
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
    title::Union{Nothing, AbstractString} = nothing,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    core_y = _plot_core_indices_2d(mesh.y_cell_sizes)
    column_range = show_padding ? (1:size(resistivity, 2)) : core_y
    y_edges = mesh.y_nodes[first(column_range):(last(column_range) + 1)] ./ 1000
    row_range = show_air ? (1:size(resistivity, 1)) : ((mesh.n_air_cells + 1):size(resistivity, 1))
    z_edges = show_air ? (mesh.z_nodes ./ 1000) : (mesh.z_nodes[(mesh.n_air_cells + 1):end] ./ 1000)
    rho_plot = log10.(resistivity[row_range, column_range])
    show_air && (rho_plot[1:mesh.n_air_cells, :] .= NaN)      # air drawn as nan_color

    figure = Figure(size = (1100, 650))
    axis = Axis(
        figure[1, 1],
        xlabel = "Offset (km)",
        ylabel = "Depth (km)",
        yreversed = true,
        title = something(title, show_air ? "2D resistivity model with air" : "2D resistivity model"),
    )
    heatmap = heatmap!(axis, y_edges, z_edges, rho_plot', colormap = :Spectral, colorrange = resistivity_log10_range,
                       nan_color = :aliceblue)
    Colorbar(figure[1, 2], heatmap, label = "log10(ρ)")
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

    save(output_path, figure)
    String(output_path)
end

"""
    PlotModel2D(model_path; output_path, show_grid=false, show_padding=true, maximum_depth_km=Inf, resistivity_log10_range=(0.0, 4.0))

Inputs:
- 2D model path, output image path, and plotting controls.

Output:
- `String`: Path to the written plot.

Description:
- Loads a saved 2D model and writes the standard model plot.
"""
function PlotModel2D(
    model_path::AbstractString;
    output_path::AbstractString,
    show_grid::Bool = false,
    show_padding::Bool = true,
    maximum_depth_km::Real = Inf,
    resistivity_log10_range::Tuple{Float64, Float64} = (0.0, 4.0),
)
    model = load_model2d(model_path)
    mesh = build_mesh_from_model2d(model; frequencies = [1.0], receiver_positions = Float64[])
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
    plot_mt2d_mesh(mesh; output_path, region=:full, background_resistivity=100.0)

Inputs:
- `mesh`: 2D MT mesh.
- `output_path`: Output image path.
- `region`: `:full` = whole mesh with air and padding; `:core` = uniform core only.
- `background_resistivity`: Resistivity for the skin-depth marker lines.

Output:
- `String`: Path to the written plot.

Description:
- Draws the cell edges, shades the air, outlines the uniform core, and marks the skin
  depth of the lowest and highest frequency.
"""
function plot_mt2d_mesh(
    mesh::MT2DMesh;
    output_path::AbstractString,
    region::Symbol = :full,
    background_resistivity::Real = 100.0,
)
    region in (:full, :core) || error("region must be :full or :core")
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    #---------- core extent ----------
    na = mesh.n_air_cells
    core_y = _plot_core_indices_2d(mesh.y_cell_sizes)
    core_z = na .+ _plot_core_indices_2d(mesh.z_cell_sizes[na+1:end])
    y_core = (mesh.y_nodes[first(core_y)], mesh.y_nodes[last(core_y)+1]) ./ 1000
    z_core = mesh.z_nodes[last(core_z)+1] / 1000
    δ_max = mt2d_skin_depth(background_resistivity, minimum(mesh.frequencies)) / 1000
    δ_min = mt2d_skin_depth(background_resistivity, maximum(mesh.frequencies)) / 1000

    #---------- visible edges ----------
    ys, zs = mesh.y_nodes ./ 1000, mesh.z_nodes ./ 1000
    if region == :core
        ys = ys[first(core_y):last(core_y)+1]
        zs = zs[na+1:last(core_z)+1]
    end

    ny, nz = length(mesh.y_cell_sizes), length(mesh.z_cell_sizes) - na
    title = region == :full ?
        @sprintf("2D mesh: %d × %d earth cells, %d air layers", ny, nz, na) :
        @sprintf("2D mesh core: dz = %.0f m to %.2f km, dy = %.0f m",
                 mesh.z_cell_sizes[na+1], z_core, mesh.y_cell_sizes[first(core_y)])
    figure = Figure(size = (1100, 650))
    axis = Axis(figure[1, 1], xlabel = "Offset (km)", ylabel = "Depth (km)", yreversed = true, title = title)

    region == :full && poly!(axis, Rect(ys[1], zs[1], ys[end] - ys[1], -zs[1]), color = :aliceblue)
    vlines!(axis, ys, color = (:black, 0.35), linewidth = 0.6)
    hlines!(axis, zs, color = (:black, 0.35), linewidth = 0.6)
    region == :full && lines!(axis, [y_core[1], y_core[2], y_core[2], y_core[1], y_core[1]],
                              [0, 0, z_core, z_core, 0], color = :firebrick, linewidth = 2, label = "uniform core")
    hlines!(axis, [δ_max], color = :darkorange, linestyle = :dash, linewidth = 2,
            label = @sprintf("δ(f_min) = %.1f km", δ_max))
    hlines!(axis, [δ_min], color = :teal, linestyle = :dot, linewidth = 2,
            label = @sprintf("δ(f_max) = %.2f km", δ_min))
    hlines!(axis, [0.0], color = :black, linewidth = 1.5)
    scatter!(axis, mesh.receiver_positions ./ 1000, zeros(length(mesh.receiver_positions));
             marker = :dtriangle, markersize = 12, color = :black)

    xlims!(axis, ys[1], ys[end])
    ylims!(axis, zs[end], min(zs[1], 0.0))
    Legend(figure[1, 2], axis, @sprintf("ρ_bg = %g Ω·m", background_resistivity), framevisible = false)
    save(output_path, figure)
    String(output_path)
end

#---------- inversion convergence ----------

"""
    plot_inv2d_convergence(history; output_path, target_rms=0.0)

Inputs:
- `history`: `Inv2DResult.history`.
- `output_path`: Output image path.
- `target_rms`: Target line on the rms panel; 0 = none.

Output:
- `String`: Path to the written plot.

Description:
- Plots rms and the objective terms (total, data, regularization) per iteration.
"""
function plot_inv2d_convergence(history::AbstractVector; output_path::AbstractString, target_rms::Real = 0.0)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    it = [h.iteration for h in history]
    figure = Figure(size = (1100, 450))
    ax1 = Axis(figure[1, 1], xlabel = "Iteration", ylabel = "RMS", title = "Data misfit")
    scatterlines!(ax1, it, [h.rms for h in history], color = :navy)
    target_rms > 0 && hlines!(ax1, [target_rms], color = :gray, linestyle = :dash)

    ax2 = Axis(figure[1, 2], xlabel = "Iteration", ylabel = "Value", yscale = log10, title = "Objective terms")
    scatterlines!(ax2, it, [h.objective for h in history], label = "objective")
    scatterlines!(ax2, it, [h.chi2 / 2 for h in history], label = "χ²/2")
    # regularization is exactly zero at the start model, skip it on a log axis
    reg = [(h.iteration, h.regularization) for h in history if h.regularization > 0]
    isempty(reg) || scatterlines!(ax2, first.(reg), last.(reg), label = "‖R(m - m_ref)‖²/2")
    axislegend(ax2, position = :rt)

    save(output_path, figure)
    String(output_path)
end

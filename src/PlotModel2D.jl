# 2D MT model plots
# Author: @pankajkmishra
# Resistivity sections and the mesh layout with its skin-depth core

using CairoMakie

"""
    plot_mt2d_model(mesh, resistivity; output_path, show_air=false, show_grid=false, show_padding=true,
                    maximum_depth_km=Inf, resistivity_log10_range=(0.0, 4.0), annotation=nothing,
                    water=nothing, log10_values=false, colormap=:Spectral,
                    colorbar_label="log10 ρ (Ω·m)") -> path

Resistivity section, depth below the model top, with the stations on the ground. Air,
topographic air included, is blank; with topography the ground is drawn as a line.
`water` is a model-shaped mask of water cells, shaded blue. `annotation` is a short
label drawn inside the section, top left, e.g. an iteration number; there is no title.
With `log10_values` the matrix is plotted as given, e.g. a log10 standard deviation.
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
    water::Union{Nothing, AbstractMatrix{Bool}} = nothing,
    log10_values::Bool = false,
    colormap = :Spectral,
    colorbar_label::AbstractString = "log10 ρ (Ω·m)",
)
    CairoMakie.activate!()

    core_y = _core_range(mesh.y_cell_sizes)
    column_range = show_padding ? (1:size(resistivity, 2)) : core_y
    y_edges = mesh.y_nodes[first(column_range):(last(column_range) + 1)] ./ 1000
    row_range = show_air ? (1:size(resistivity, 1)) : ((mesh.n_air_cells + 1):size(resistivity, 1))
    z_edges = show_air ? (mesh.z_nodes ./ 1000) : (mesh.z_nodes[(mesh.n_air_cells + 1):end] ./ 1000)
    rho_plot = log10_values ? Float64.(resistivity[row_range, column_range]) : log10.(resistivity[row_range, column_range])
    rho_plot[mt2d_air_mask(mesh)[row_range, column_range]] .= NaN      # air drawn as nan_color

    figure = Figure(size = (1100, 650))
    axis = _mt_axis(figure[1, 1]; xlabel = "Offset (km)", ylabel = "Depth (km)", yreversed = true)
    heatmap = heatmap!(axis, y_edges, z_edges, rho_plot', colormap = colormap, colorrange = resistivity_log10_range,
                       nan_color = :aliceblue)
    Colorbar(figure[1, 2], heatmap, label = colorbar_label, labelfont = :regular)
    xlims!(axis, minimum(y_edges), maximum(y_edges))
    if water !== nothing && any(water)
        wet = [water[iz, iy] ? 1.0 : NaN for iz in row_range, iy in column_range]
        heatmap!(axis, y_edges, z_edges, wet', colormap = [:lightskyblue, :lightskyblue], nan_color = :transparent)
    end
    _mt2d_ground_line!(axis, mesh, column_range)

    scatter!(
        axis,
        mesh.receiver_positions ./ 1000,
        mt2d_receiver_depths(mesh) ./ 1000;
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

# staircase ground surface in km, only when there is topography
function _mt2d_ground_line!(axis, mesh, columns = eachindex(mesh.y_cell_sizes))
    topo = mt2d_topo_air(mesh)
    any(>(0), topo) || return nothing
    z0 = mesh.z_nodes[mesh.n_air_cells+1]
    depth = [mesh.z_nodes[mesh.n_air_cells+1+topo[iy]] - z0 for iy in columns] ./ 1000
    stairs!(axis, mesh.y_nodes[first(columns):last(columns)+1] ./ 1000, vcat(depth, depth[end]);
            step = :post, color = :black, linewidth = 1.2)
    nothing
end

# water cells of an earth-cell mask filled blue with their bed (bathymetry) outlined; offset and depth below
# the model top in metres divided by `yunit` and `zunit` (1000 = km)
function _mt2d_water!(axis, mesh, water; yunit = 1000, zunit = 1000)
    (water === nothing || !any(water)) && return nothing
    na = mesh.n_air_cells
    y, z = mesh.y_nodes ./ yunit, (mesh.z_nodes .- mesh.z_nodes[na+1]) ./ zunit
    wet = findall(water)
    poly!(axis, [Rect2(y[i[2]], z[na+i[1]], y[i[2]+1] - y[i[2]], z[na+i[1]+1] - z[na+i[1]]) for i in wet];
          color = :lightskyblue)
    bed = [any(water[:, iy]) ? z[na+findlast(water[:, iy])+1] : NaN for iy in axes(water, 2)]
    seg = Point2f[]
    for iy in eachindex(bed)
        isnan(bed[iy]) && continue
        append!(seg, [Point2f(y[iy], bed[iy]), Point2f(y[iy+1], bed[iy])])
        iy < length(bed) && !isnan(bed[iy+1]) && append!(seg, [Point2f(y[iy+1], bed[iy]), Point2f(y[iy+1], bed[iy+1])])
    end
    linesegments!(axis, seg; color = :steelblue4, linewidth = 1.2)
    nothing
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
    plot_mt2d_mesh(mesh; output_path, region=:full, background_resistivity=100.0, water=nothing) -> path

Cell edges with the air shaded, the core outlined and the skin depths of the lowest and
highest frequency in `background_resistivity` marked. `water` is a model-shaped mask of
water cells (lakes, sea), filled blue with the bed line. The core is found as in 3D:
laterally `core_indices` on the cell widths, and down to the layer boundary nearest the
skin depth of the lowest frequency, as `z_indices_for_max_depth`; `region = :core` shows
it only.
"""
function plot_mt2d_mesh(
    mesh::MT2DMesh;
    output_path::AbstractString,
    region::Symbol = :full,
    background_resistivity::Real = 100.0,
    water::Union{Nothing, AbstractMatrix{Bool}} = nothing,
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
    _mt2d_water!(axis, mesh, water)
    vlines!(axis, ys, color = (:black, 0.35), linewidth = 0.6)
    hlines!(axis, zs, color = (:black, 0.35), linewidth = 0.6)
    region == :full && lines!(axis, [y_core[1], y_core[2], y_core[2], y_core[1], y_core[1]],
                              [0, 0, z_core, z_core, 0], color = :firebrick, linewidth = 2, label = "core")
    hlines!(axis, [δ_max], color = :darkorange, linestyle = :dash, linewidth = 2,
            label = @sprintf("δ(f_min) = %.1f km", δ_max))
    hlines!(axis, [δ_min], color = :teal, linestyle = :dot, linewidth = 2,
            label = @sprintf("δ(f_max) = %.2f km", δ_min))
    hlines!(axis, [0.0], color = :black, linewidth = 1.5)
    _mt2d_ground_line!(axis, mesh)
    scatter!(axis, mesh.receiver_positions ./ 1000, mt2d_receiver_depths(mesh) ./ 1000;
             marker = :dtriangle, markersize = 12, color = :black)

    xlims!(axis, ys[1], ys[end])
    ylims!(axis, max(zs[end], 1.02δ_max), min(zs[1], 0.0))       # the core may stop just above δ(f_min)
    Legend(figure[1, 2], axis, @sprintf("ρ_bg = %g Ω·m", background_resistivity), framevisible = false,
           titlefont = :regular, labelfont = :regular)
    _mt_save(output_path, figure)
end

#---------- model file grid ----------

# earth-only mesh for plotting a model file on its own grid, topographic air from the tags
function _mt2d_earth_mesh(model::ModelFile2D; receivers = Float64[])
    topo = _mt2d_model_topo_air(model.resistivity)
    MT2DMesh(y_nodes = model.origin[2] .+ vcat(0.0, cumsum(model.y_cell_sizes)), z_nodes = vcat(0.0, cumsum(model.z_cell_sizes)),
             y_cell_sizes = model.y_cell_sizes, z_cell_sizes = model.z_cell_sizes, receiver_positions = Float64.(receivers),
             frequencies = [1.0], n_air_cells = 0, topo_air = any(>(0), topo) ? topo : Int[])
end

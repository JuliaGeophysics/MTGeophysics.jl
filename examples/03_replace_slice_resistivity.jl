# Example: replace deep resistivity values below a cutoff depth.
# Author: @pankajkmishra
# This script applies bulk edits to deeper layers and opens the depth-slice viewer for review.
# Use it for quick scenario edits before saving a modified model.

using Pkg

using GLMakie
using Statistics
using Dates

include(joinpath(dirname(@__DIR__), "src", "Model.jl"))
include(joinpath(dirname(@__DIR__), "src", "PlotModel.jl"))
include(joinpath(@__DIR__, "02_plot_depth_slices.jl"))

model_file = joinpath(@__DIR__, "cascadiaInv", "cascad_half_inverse.ws")

cutoff_depth = 40000.0
target_resistivity = 10000.0
transition_layers = 1

log10_scale = true
colormap = Reverse(:turbo)
with_padding = true
max_depth = nothing
resistivity_range = (0.0, 4.0)
show_grid = true
grid_color = :black
grid_linewidth = 0.5
grid_alpha = 0.3
pad_tol = 0.5

function get_layer_depths(z_centers::AbstractVector)
    z_edges = edges_from_centers(z_centers)
    cumsum(diff(z_edges))
end

function find_cutoff_layer(layer_depths::AbstractVector, cutoff_depth::Real)
    idx = findfirst(>=(cutoff_depth), layer_depths)
    isnothing(idx) ? length(layer_depths) + 1 : idx
end

function modify_deep_resistivity!(A::Array{<:Real,3}, z_centers::AbstractVector;
                                   cutoff_depth::Real = 100000.0,
                                   target_resistivity::Real = 1000.0,
                                   transition_layers::Int = 3)
    layer_depths = get_layer_depths(z_centers)
    cutoff_layer = find_cutoff_layer(layer_depths, cutoff_depth)

    nz = size(A, 3)

    if cutoff_layer > nz
        println("Warning: cutoff_depth $(cutoff_depth) m is deeper than model. No modification made.")
        return A, cutoff_layer, layer_depths
    end

    println("Modifying resistivity:")
    println("  Cutoff depth: $(cutoff_depth) m")
    println("  Cutoff layer: $(cutoff_layer) / $(nz)")
    println("  Target resistivity: $(target_resistivity) ohm-m")
    println("  Transition layers: $(transition_layers)")

    if transition_layers > 0 && cutoff_layer > 1
        trans_start = max(1, cutoff_layer - transition_layers)
        trans_end = cutoff_layer - 1

        for (i, k) in enumerate(trans_start:trans_end)
            weight = i / (transition_layers + 1)
            A[:, :, k] .= (1 - weight) .* A[:, :, k] .+ weight * target_resistivity
        end
        println("  Transition applied to layers $(trans_start) - $(trans_end)")
    end

    for k in cutoff_layer:nz
        A[:, :, k] .= target_resistivity
    end
    println("  Replaced layers $(cutoff_layer) - $(nz) with target value")

    return A, cutoff_layer, layer_depths
end

function main()
    if !isfile(model_file)
        println("ERROR: Model file not found: $model_file")
        return nothing, nothing
    end

    println("Loading ModEM model: $model_file")
    M = load_model_modem(model_file)

    println("Model loaded:")
    println("  Grid size: $(M.nx) × $(M.ny) × $(M.nz)")

    A_modified = copy(M.A)

    A_modified, cutoff_layer, layer_depths = modify_deep_resistivity!(
        A_modified, M.cz;
        cutoff_depth = cutoff_depth,
        target_resistivity = target_resistivity,
        transition_layers = transition_layers
    )

    M_modified = (
        A = A_modified,
        x = M.x, y = M.y, z = M.z,
        cx = M.cx, cy = M.cy, cz = M.cz,
        nx = M.nx, ny = M.ny, nz = M.nz,
        npad = M.npad
    )

    model_name = splitext(basename(model_file))[1] * "_modified"

    println("\nCreating viewer for modified model...")

    fig, parts = depth_slice_viewer(M_modified;
        model_name = model_name,
        log10scale = log10_scale,
        cmap = colormap,
        withPadding = with_padding,
        max_depth = max_depth,
        pad_tol = pad_tol,
        resistivity_range = resistivity_range,
        show_grid = show_grid,
        grid_color = grid_color,
        grid_linewidth = grid_linewidth,
        grid_alpha = grid_alpha
    )

    save_grid = fig[4, 1:2] = GridLayout()
    btn_save = Button(save_grid[1, 1], label = "Save Model", fontsize = 14)
    save_label = Observable("")
    Label(save_grid[1, 2], save_label, fontsize = 12, color = :green)

    on(btn_save.clicks) do _
        outname = splitext(model_file)[1] * "_modified.rho"
        write_model_modem(outname, M.dx, M.dy, M.dz, A_modified, M.origin)
        save_label[] = "Saved: $(basename(outname))"
    end

    println("\nViewer ready!")
    println("  Cutoff at layer $(cutoff_layer), depth $(round(cutoff_depth/1000, digits=1)) km")
    println("  Click 'Save Model' to write the modified model to disk.")

    screen = display(fig)
    println("\nClose the figure window to exit...")
    wait(screen)

    return fig, parts, (M=M, A_modified=A_modified)
end

fig, parts, data = main()

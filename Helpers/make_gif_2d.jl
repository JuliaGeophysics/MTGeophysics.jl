#=
make_gif_2d.jl — Generate a convergence GIF from VFSA2DMT snapshot models.

Reads the averaged iteration snapshots (`avg_iter_NNNNN.rho`) written by
`run_mt2d_vfsa` in the `snapshots/` subdirectory of a run, renders each
frame as a 2-D resistivity heatmap, and assembles them into an animated GIF.

Usage (from the MTGeophysics.jl project root):

    julia --project=. Helpers/make_gif_2d.jl <run_dir> [output.gif] [--fps N] [--depth_km D] [--rho_range lo,hi]

Arguments:
    run_dir       Path to the VFSA run directory (must contain snapshots/).
    output.gif    Optional output path (default: <run_dir>/convergence.gif).

Options:
    --fps N            Frames per second (default: 4).
    --depth_km D       Maximum depth in km to display (default: Inf).
    --rho_range lo,hi  log10(ρ) colour range (default: 0.0,4.0).
=#

using Pkg
Pkg.activate(dirname(@__DIR__))
using MTGeophysics
using CairoMakie

function make_convergence_gif(
    run_dir::AbstractString;
    output_path::Union{Nothing, AbstractString} = nothing,
    fps::Int = 4,
    maximum_depth_km::Float64 = Inf,
    resistivity_log10_range::Tuple{Float64, Float64} = (0.0, 4.0),
)
    snapshot_dir = joinpath(run_dir, "snapshots")
    isdir(snapshot_dir) || error("snapshots/ directory not found in $run_dir")

    # Discover snapshot files sorted by iteration number.
    files = filter(f -> startswith(f, "avg_iter_") && endswith(f, ".rho"), readdir(snapshot_dir))
    sort!(files)
    isempty(files) && error("no avg_iter_*.rho files found in $snapshot_dir")

    gif_path = output_path === nothing ? joinpath(run_dir, "convergence.gif") : String(output_path)
    mkpath(dirname(gif_path))

    # Need an observed data file to reconstruct the mesh for plotting.
    obs_path = joinpath(run_dir, "observed.obs")
    if !isfile(obs_path)
        # Fall back to any .obs file in the run dir.
        obs_candidates = filter(f -> endswith(f, ".obs"), readdir(run_dir))
        isempty(obs_candidates) && error("no observed data (.obs) file found in $run_dir")
        obs_path = joinpath(run_dir, first(obs_candidates))
    end
    observed_data = load_data2d(obs_path)

    # Load the first snapshot to build the mesh.
    first_model = load_model2d(joinpath(snapshot_dir, first(files)))
    mesh = build_mesh_from_model2d(
        first_model;
        frequencies = observed_data.frequencies,
        receiver_positions = observed_data.receivers,
    )

    println("Rendering $(length(files)) frames …")
    CairoMakie.activate!()
    figure = Figure(size = (1100, 650))

    # Pre-compute mesh edges (always show full grid including padding).
    column_range = 1:size(first_model.resistivity, 2)
    y_edges = mesh.y_nodes[first(column_range):(last(column_range) + 1)] ./ 1000
    row_range = (mesh.n_air_cells + 1):size(first_model.resistivity, 1)
    z_edges = mesh.z_nodes[(mesh.n_air_cells + 1):end] ./ 1000
    depth_limit_km = isfinite(maximum_depth_km) ? min(maximum_depth_km, maximum(z_edges)) : maximum(z_edges)

    record(figure, gif_path, eachindex(files); framerate = fps) do frame_index
        empty!(figure)
        model = load_model2d(joinpath(snapshot_dir, files[frame_index]))
        rho_plot = log10.(model.resistivity[row_range, column_range])

        iter_str = replace(replace(files[frame_index], "avg_iter_" => ""), ".rho" => "")
        ax = Axis(
            figure[1, 1],
            xlabel = "Offset (km)",
            ylabel = "Depth (km)",
            yreversed = true,
            title = "Chain-mean best model — iteration $iter_str",
        )
        hm = heatmap!(ax, y_edges, z_edges, rho_plot', colormap = :Spectral, colorrange = resistivity_log10_range)
        Colorbar(figure[1, 2], hm, label = "log₁₀(ρ)")
        xlims!(ax, minimum(y_edges), maximum(y_edges))
        ylims!(ax, depth_limit_km, 0.0)
        scatter!(
            ax,
            mesh.receiver_positions ./ 1000,
            fill(0.0, length(mesh.receiver_positions));
            marker = :dtriangle,
            markersize = 12,
            color = :black,
        )
    end
    println("GIF written to $gif_path")
    gif_path
end

# --- CLI entry-point ---
if abspath(PROGRAM_FILE) == @__FILE__
    let args = copy(ARGS)
        isempty(args) && error("Usage: julia make_gif_2d.jl <run_dir> [output.gif] [--fps N] [--depth_km D] [--rho_range lo,hi]")

        run_dir = popfirst!(args)
        output = nothing
        fps = 4
        depth = Inf
        rho_lo, rho_hi = 0.0, 4.0

        while !isempty(args)
            arg = popfirst!(args)
            if arg == "--fps"
                fps = parse(Int, popfirst!(args))
            elseif arg == "--depth_km"
                depth = parse(Float64, popfirst!(args))
            elseif arg == "--rho_range"
                parts = split(popfirst!(args), ",")
                rho_lo = parse(Float64, parts[1])
                rho_hi = parse(Float64, parts[2])
            else
                output = arg
            end
        end

        make_convergence_gif(
            run_dir;
            output_path = output,
            fps = fps,
            maximum_depth_km = depth,
            resistivity_log10_range = (rho_lo, rho_hi),
        )
    end
end

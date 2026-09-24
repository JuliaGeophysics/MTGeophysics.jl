# Convergence GIF of a 2D VFSA run
# Author: @pankajkmishra
# Each frame is the cross-chain mean (log10) of the chains' best models at one snapshot iteration,
# read from run_dir/vfsa/chain_XX/best_iter_NNNNN.rho (set "Snapshot interval" in the VFSA control)
# Usage: julia --project=. helpers/make_gif_2D.jl <run_dir> [output.gif] [--fps N] [--depth_km D] [--rho_range lo,hi]

using MTGeophysics
using CairoMakie
using Printf

function make_convergence_gif(run_dir::AbstractString; output_path = nothing, fps::Int = 4,
                              maximum_depth_km::Float64 = Inf, resistivity_log10_range = (0.0, 4.0))
    vdir = joinpath(run_dir, "vfsa")
    chains = filter(d -> startswith(d, "chain_"), readdir(vdir))
    isempty(chains) && error("no chain_XX directories in $vdir")
    snaps = sort(unique(f for c in chains for f in readdir(joinpath(vdir, c)) if startswith(f, "best_iter_")))
    isempty(snaps) && error("no best_iter_*.rho snapshots in $vdir; set 'Snapshot interval' in the VFSA control")
    gif_path = output_path === nothing ? joinpath(run_dir, "plots", "Convergence.gif") : String(output_path)
    mkpath(dirname(gif_path))

    frames = mktempdir() do tmp
        paths = String[]
        for (k, snap) in enumerate(snaps)
            models = [ReadModel2D(joinpath(vdir, c, snap)) for c in chains if isfile(joinpath(vdir, c, snap))]
            ens = mt2d_ensemble([m.resistivity for m in models])
            ρ = 10.0 .^ ens.mean
            ρ[models[1].resistivity .> MTGeophysics.MT2D_AIR_THRESHOLD] .= MTGeophysics.MT2D_AIR_TAG
            model = ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = models[1].y_cell_sizes,
                                z_cell_sizes = models[1].z_cell_sizes, resistivity = ρ, n_air_cells = 0,
                                origin = models[1].origin, rotation = 0.0, format = "LOGE")
            it = parse(Int, match(r"(\d+)", snap).captures[1])
            push!(paths, plot_mt2d_model(MTGeophysics._mt2d_earth_mesh(model), ρ; output_path = joinpath(tmp, @sprintf("f%05d.png", k)),
                                         show_padding = false, maximum_depth_km, resistivity_log10_range,
                                         annotation = "iteration $it, $(length(models)) chains"))
        end
        [load(p) for p in paths]
    end
    record(Figure(size = reverse(size(frames[1])) .÷ 3), gif_path, eachindex(frames); framerate = fps) do k
        empty!(current_figure())
        image!(Axis(current_figure()[1, 1]; aspect = DataAspect()), rotr90(frames[k]))
        hidedecorations!(current_axis()); hidespines!(current_axis())
    end
    gif_path
end

function main(args = ARGS)
    isempty(args) && error("usage: julia --project=. helpers/make_gif_2D.jl <run_dir> [output.gif] [--fps N] [--depth_km D] [--rho_range lo,hi]")
    opts = Dict{String, String}()
    positional = String[]
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            opts[args[i]] = args[i+1]; i += 2
        else
            push!(positional, args[i]); i += 1
        end
    end
    lohi = parse.(Float64, split(get(opts, "--rho_range", "0,4"), ','))
    path = make_convergence_gif(positional[1]; output_path = get(positional, 2, nothing), fps = parse(Int, get(opts, "--fps", "4")),
                                maximum_depth_km = parse(Float64, get(opts, "--depth_km", "Inf")),
                                resistivity_log10_range = (lohi[1], lohi[2]))
    println("GIF = ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# 2D mesh tool window
# Author: @pankajkmishra
# Front end of MakeMesh2D: one slider per setting, a live mesh preview with the topography, water, stations and
# skin depths, a core/full view toggle as in MakeMesh3D, the mesh summary and advice, and Save to write the inputs. The window is built on any Makie
# backend, so CairoMakie can drive it headless; _make_mesh2D_gui shows it with GLMakie

# the window and its widgets; refresh! rebuilds the preview, the sliders and Save are wired up
function _makemesh2d_window(data, topo, water, p0, out_dir, ctrls)
    specs = [
        (:cell_width_frac, "Cell width / station spacing", 0.1:0.05:1.0),
        (:core_margin_cells, "Core margin (cells)", 0:1:20),
        (:n_pad, "Padding cells", 4:1:30),
        (:pad_factor, "Padding growth", 1.1:0.05:2.0),
        (:first_layer_div, "δ(f_max) / first layer", 2.0:0.5:20.0),
        (:vertical_factor, "Layer growth", 1.02:0.01:1.4),
        (:depth_mult, "Depth × δ(f_min)", 1.0:0.25:6.0),
        (:air_layers, "Air layers", 4:1:20),
        (:cov_smoothing, "Covariance smoothing", 0.0:0.05:0.9),
    ]
    fig = Figure(size = (1500, 900))
    controls = fig[1:2, 1] = GridLayout()
    grid = SliderGrid(controls[1, 1], [(label = l, range = r, startvalue = getproperty(p0, k)) for (k, l, r) in specs]...;
                      width = 380)
    buttons = controls[2, 1] = GridLayout()
    save_button = Button(buttons[1, 1], label = "Save inputs", tellwidth = false)
    view_button = Button(buttons[1, 2], label = "Show full", tellwidth = false)
    show_core = Ref(true)
    info = Label(controls[3, 1], "", tellwidth = false, justification = :left, word_wrap = true, width = 380)
    axis = Axis(fig[1, 2]; xlabel = "Offset (km)", ylabel = "Depth (km)", yreversed = true,
                xgridvisible = false, ygridvisible = false)
    # the ground, water and stations near the surface, where the whole-mesh view is too coarse
    surface = Axis(fig[2, 2]; xlabel = "Offset (km)", ylabel = "Depth (m)", yreversed = true,
                   xgridvisible = false, ygridvisible = false)
    colsize!(fig.layout, 1, Fixed(400))
    rowsize!(fig.layout, 2, Relative(0.35))

    params() = merge(p0, NamedTuple{Tuple(first.(specs))}(Tuple(s.value[] for s in grid.sliders)))
    built = Ref{Any}(nothing)
    function refresh!()
        b = try
            _makemesh2d_build(data, topo, water, params())
        catch err
            info.text[] = "✗ " * sprint(showerror, err)
            return
        end
        built[] = b
        empty!(axis)
        m = b.mesh
        wet = b.mask .== MT2D_MASK_WATER          # not `water`, which the closure shares with the build
        ys, zs = m.y_nodes ./ 1000, (m.z_nodes .- m.z_nodes[m.n_air_cells+1]) ./ 1000
        _mt2d_water!(axis, m, wet)
        vlines!(axis, ys, color = (:black, 0.25), linewidth = 0.5)
        hlines!(axis, zs[m.n_air_cells+1:end], color = (:black, 0.25), linewidth = 0.5)
        _mt2d_ground_line!(axis, m)
        hlines!(axis, [b.δmax / 1000], color = :darkorange, linestyle = :dash)
        scatter!(axis, m.receiver_positions ./ 1000, mt2d_receiver_depths(m) ./ 1000; marker = :dtriangle,
                 color = :black, markersize = 10)
        core = _core_range(m.y_cell_sizes)
        if show_core[]
            xlims!(axis, ys[first(core)], ys[last(core)+1])
            ylims!(axis, 1.02 * b.δmax / 1000, -0.02 * b.δmax / 1000)
        else
            xlims!(axis, ys[1], ys[end])
            ylims!(axis, zs[end], -0.02 * zs[end])
        end

        empty!(surface)
        na, zm = m.n_air_cells, 1000 .* zs
        _mt2d_water!(surface, m, wet; zunit = 1)
        vlines!(surface, ys, color = (:black, 0.25), linewidth = 0.5)
        hlines!(surface, zm[na+1:end], color = (:black, 0.25), linewidth = 0.5)
        ground = zm[na .+ 1 .+ mt2d_topo_air(m)]
        stairs!(surface, ys, vcat(ground, ground[end]); step = :post, color = :black, linewidth = 1.2)
        scatter!(surface, m.receiver_positions ./ 1000, mt2d_receiver_depths(m); marker = :dtriangle,
                 color = :black, markersize = 10)
        xlims!(surface, ys[first(core)], ys[last(core)+1])
        deepest = maximum(ground[core])
        ylims!(surface, deepest + 3 * (zm[na+2] - zm[na+1]) + 1, -1)
        info.text[] = b.summary * "\n" * (isempty(b.notes) ? "✓ mesh looks well sized" : join("• " .* b.notes, "\n"))
    end
    for s in grid.sliders
        on(_ -> refresh!(), s.value)
    end
    on(view_button.clicks) do _
        show_core[] = !show_core[]
        view_button.label[] = show_core[] ? "Show full" : "Show core"
        refresh!()
    end
    on(save_button.clicks) do _
        built[] === nothing && return
        paths = _makemesh2d_write(built[], out_dir, params(), ctrls, topo !== nothing)
        info.text[] = "Saved to $(out_dir)\n" * built[].summary
        foreach(println, values(paths))
    end
    refresh!()
    (; fig, grid, save_button, view_button, info, built, params, axis, surface)
end

function _make_mesh2D_gui(data, topo, water, p0, out_dir, ctrls)
    GLMakie.activate!()
    w = _makemesh2d_window(data, topo, water, p0, out_dir, ctrls)
    display(w.fig)
    w.fig
end

# 2D MT mesh and model
# Author: @pankajkmishra
# Profile mesh geometry, skin-depth vertical core, model builders, and model file I/O

using LinearAlgebra
using Printf

"""
    MT2DMesh

Inputs:
- Horizontal and vertical mesh nodes, cell sizes, receiver positions, frequencies, and air-cell count.

Output:
- `MT2DMesh`: Container for a 2D MT profile mesh.

Description:
- Stores the profile mesh geometry and survey definition used by the 2D workflows.
"""
Base.@kwdef struct MT2DMesh
    y_nodes::Vector{Float64}
    z_nodes::Vector{Float64}
    y_cell_sizes::Vector{Float64}
    z_cell_sizes::Vector{Float64}
    receiver_positions::Vector{Float64}
    frequencies::Vector{Float64}
    n_air_cells::Int
end

"""
    ModelFile2D

Inputs:
- Model metadata, mesh spacings, resistivity values, air-cell count, origin, and rotation.

Output:
- `ModelFile2D`: Parsed 2D model file container.

Description:
- Stores a 2D model file exactly as needed for round-tripping between disk and the forward solver.
"""
Base.@kwdef struct ModelFile2D
    title::String
    x_cell_sizes::Vector{Float64}
    y_cell_sizes::Vector{Float64}
    z_cell_sizes::Vector{Float64}
    resistivity::Matrix{Float64}
    n_air_cells::Int
    origin::Vector{Float64}
    rotation::Float64
    format::String
    path::String = ""
end

const μ₀_2D = 4π * 1e-7

"""
    mt2d_y_centers(mesh)

Inputs:
- `mesh`: 2D MT mesh.

Output:
- `Vector{Float64}`: Horizontal cell-center positions.

Description:
- Computes the horizontal cell centers of the 2D mesh.
"""
mt2d_y_centers(mesh::MT2DMesh) = 0.5 .* (mesh.y_nodes[1:end-1] .+ mesh.y_nodes[2:end])

"""
    mt2d_z_centers(mesh)

Inputs:
- `mesh`: 2D MT mesh.

Output:
- `Vector{Float64}`: Vertical cell-center positions.

Description:
- Computes the vertical cell centers of the 2D mesh.
"""
mt2d_z_centers(mesh::MT2DMesh) = 0.5 .* (mesh.z_nodes[1:end-1] .+ mesh.z_nodes[2:end])

"""
    mt2d_center_station(mesh)

Inputs:
- `mesh`: 2D MT mesh.

Output:
- `Int`: Index of the middle receiver location.

Description:
- Returns the central survey station used by plotting and smoke tests.
"""
mt2d_center_station(mesh::MT2DMesh) = cld(length(mesh.receiver_positions), 2)

"""
    mt2d_skin_depth(resistivity, frequency)

Inputs:
- `resistivity`: Halfspace resistivity in ohm metres.
- `frequency`: Frequency in Hz.

Output:
- `Float64`: Skin depth `δ = sqrt(2ρ/(ωμ₀)) ≈ 503·sqrt(ρ/f)` in metres.
"""
mt2d_skin_depth(resistivity::Real, frequency::Real) = sqrt(2 * resistivity / (2π * frequency * μ₀_2D))

"""
    mt2d_skin_depth_layers(frequencies; background_resistivity=100.0, z_core_cell=nothing,
                           z_core_skin_depths=1.0, z_bottom_skin_depths=4.0,
                           max_core_layers=80, pad_factor=1.25)

Inputs:
- `frequencies`: Survey frequencies in Hz.
- `background_resistivity`: Reference halfspace resistivity for the skin depths.
- `z_core_cell`: Uniform core thickness; `nothing` = `max(δ_min/3, δ_max/max_core_layers)`.
- `z_core_skin_depths`: Depth of the uniform core in skin depths of the lowest frequency.
- `z_bottom_skin_depths`: Total mesh depth in skin depths of the lowest frequency.
- `pad_factor`: Geometric growth of the padding layers below the core.

Output:
- `Vector{Float64}`: Ground layer thicknesses, surface downward.

Description:
- The core is regular (constant dz) from the surface down to at least `z_core_skin_depths`
  skin depths of the longest period; below it the layers grow geometrically until the
  mesh bottom reaches `z_bottom_skin_depths` skin depths, far enough for the 1D
  boundary fields to have decayed.
"""
function mt2d_skin_depth_layers(
    frequencies::AbstractVector{<:Real};
    background_resistivity::Real = 100.0,
    z_core_cell::Union{Nothing, Real} = nothing,
    z_core_skin_depths::Real = 1.0,
    z_bottom_skin_depths::Real = 4.0,
    max_core_layers::Integer = 80,
    pad_factor::Real = 1.25,
)
    background_resistivity > 0 || error("background_resistivity must be positive")
    z_core_skin_depths > 0 || error("z_core_skin_depths must be positive")
    z_bottom_skin_depths >= z_core_skin_depths || error("z_bottom_skin_depths must be >= z_core_skin_depths")
    pad_factor > 1 || error("pad_factor must be greater than 1")

    δ_min = mt2d_skin_depth(background_resistivity, maximum(frequencies))
    δ_max = mt2d_skin_depth(background_resistivity, minimum(frequencies))
    dz = something(z_core_cell, max(δ_min / 3, δ_max / max_core_layers))
    dz > 0 || error("z_core_cell must be positive")

    # uniform core, rounded up to a whole number of cells
    n_core = ceil(Int, z_core_skin_depths * δ_max / dz - 1e-9)
    layers = fill(Float64(dz), n_core)

    # geometric padding to the bottom
    z_bottom = z_bottom_skin_depths * δ_max
    depth, Δz = n_core * dz, Float64(dz)
    while depth < z_bottom
        Δz *= pad_factor
        depth += Δz
        push!(layers, Δz)
    end
    layers
end

"""
    build_mt2d_mesh(; frequencies=..., y_core_range=(-6000.0, 6000.0), y_core_cell=300.0,
                    y_padding=9000.0, pad_factor=1.25, air_top=-12000.0, air_cells=8,
                    ground_layers=nothing, background_resistivity=100.0, z_core_cell=nothing,
                    z_core_skin_depths=1.0, z_bottom_skin_depths=4.0, max_core_layers=80,
                    receiver_stride=2, receiver_positions=nothing)

Inputs:
- Frequency axis and horizontal/vertical mesh controls.
- `ground_layers`: Explicit ground layer thicknesses; `nothing` = skin-depth design
  from `mt2d_skin_depth_layers` with the `background_resistivity`, `z_core_*`,
  `z_bottom_skin_depths`, `max_core_layers` and `pad_factor` keywords.

Output:
- `MT2DMesh`: Survey mesh and receiver geometry.

Description:
- Builds the padded 2D MT profile mesh used by the forward and inversion workflows.
  By default the vertical core is regular down to one skin depth of the lowest
  frequency in the background resistivity.
"""
function build_mt2d_mesh(;
    frequencies::AbstractVector{<:Real} = collect(10 .^ range(-2, 2, length = 10)),
    y_core_range::Tuple{<:Real, <:Real} = (-6000.0, 6000.0),
    y_core_cell::Real = 300.0,
    y_padding::Real = 9000.0,
    pad_factor::Real = 1.25,
    air_top::Real = -12000.0,
    air_cells::Integer = 8,
    ground_layers::Union{Nothing, AbstractVector{<:Real}} = nothing,
    background_resistivity::Real = 100.0,
    z_core_cell::Union{Nothing, Real} = nothing,
    z_core_skin_depths::Real = 1.0,
    z_bottom_skin_depths::Real = 4.0,
    max_core_layers::Integer = 80,
    receiver_stride::Integer = 2,
    receiver_positions::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    y_core_cell > 0 || error("y_core_cell must be positive")
    y_padding > 0 || error("y_padding must be positive")
    pad_factor > 1 || error("pad_factor must be greater than 1")
    air_cells > 0 || error("air_cells must be positive")
    receiver_stride > 0 || error("receiver_stride must be positive")
    if ground_layers === nothing
        ground_layers = mt2d_skin_depth_layers(frequencies; background_resistivity, z_core_cell,
            z_core_skin_depths, z_bottom_skin_depths, max_core_layers, pad_factor)
    end
    any(Δz -> Δz <= 0, ground_layers) && error("all ground layer thicknesses must be positive")

    y1, y2 = Float64.(y_core_range)
    y_core_nodes = collect(y1:y_core_cell:y2)
    abs(y_core_nodes[end] - y2) > 1e-9 && push!(y_core_nodes, y2)

    left_nodes = Float64[]
    Δy = Float64(y_core_cell)
    y = y1
    while y - Δy > y1 - y_padding - 1e-9
        Δy *= pad_factor
        y -= Δy
        push!(left_nodes, y)
    end
    reverse!(left_nodes)

    right_nodes = Float64[]
    Δy = Float64(y_core_cell)
    y = y2
    while y + Δy < y2 + y_padding + 1e-9
        Δy *= pad_factor
        y += Δy
        push!(right_nodes, y)
    end

    y_nodes = vcat(left_nodes, y_core_nodes, right_nodes)
    y_cell_sizes = diff(y_nodes)

    z_air = collect(range(Float64(air_top), 0.0, length = air_cells + 1))
    z_ground = vcat(0.0, cumsum(Float64.(ground_layers)))
    z_nodes = vcat(z_air[1:end-1], z_ground)
    z_cell_sizes = diff(z_nodes)

    receivers = if receiver_positions === nothing
        y_receivers = 0.5 .* (y_core_nodes[1:end-1] .+ y_core_nodes[2:end])
        collect(y_receivers[1:receiver_stride:end])
    else
        Float64.(receiver_positions)
    end

    MT2DMesh(
        y_nodes = y_nodes,
        z_nodes = z_nodes,
        y_cell_sizes = y_cell_sizes,
        z_cell_sizes = z_cell_sizes,
        receiver_positions = receivers,
        frequencies = Float64.(frequencies),
        n_air_cells = air_cells,
    )
end

"""
    BuildMesh2D(; kwargs...)

Inputs:
- Keyword arguments accepted by `build_mt2d_mesh`.

Output:
- `MT2DMesh`: Survey mesh and receiver geometry.

Description:
- Public alias for `build_mt2d_mesh`.
"""
BuildMesh2D(; kwargs...) = build_mt2d_mesh(; kwargs...)

"""
    build_mt2d_block_model(mesh; background_resistivity=100.0, blocks=NamedTuple[])

Inputs:
- `mesh`: 2D MT mesh.
- `background_resistivity`: Host resistivity in ohm metres.
- `blocks`: Rectangular anomaly definitions.

Output:
- `Matrix{Float64}`: Cell resistivity model.

Description:
- Builds a 2D block model on the provided mesh.
"""
function build_mt2d_block_model(
    mesh::MT2DMesh;
    background_resistivity::Real = 100.0,
    blocks::AbstractVector = NamedTuple[],
)
    n_z = length(mesh.z_cell_sizes)
    n_y = length(mesh.y_cell_sizes)
    ρ = fill(Float64(background_resistivity), n_z, n_y)
    ρ[1:mesh.n_air_cells, :] .= 1e9

    y_centers = mt2d_y_centers(mesh)
    z_centers = mt2d_z_centers(mesh)
    for block in blocks
        y1, y2 = Float64.(block.y_range)
        z1, z2 = Float64.(block.z_range)
        ρblock = Float64(block.resistivity)
        for iy in eachindex(y_centers), iz in eachindex(z_centers)
            if y1 <= y_centers[iy] <= y2 && z1 <= z_centers[iz] <= z2
                ρ[iz, iy] = ρblock
            end
        end
    end

    ρ
end

"""
    build_mt2d_halfspace_model(mesh; background_resistivity=100.0)

Inputs:
- `mesh`: 2D MT mesh.
- `background_resistivity`: Half-space resistivity in ohm metres.

Output:
- `Matrix{Float64}`: Half-space resistivity model.

Description:
- Builds a homogeneous half-space on the given mesh.
"""
build_mt2d_halfspace_model(mesh::MT2DMesh; background_resistivity::Real = 100.0) =
    build_mt2d_block_model(mesh; background_resistivity = background_resistivity, blocks = NamedTuple[])

"""
    build_mt2d_layered_model(mesh; layer_resistivities, interface_depths)

Inputs:
- `mesh`: 2D MT mesh.
- `layer_resistivities`: Layer resistivities including the basement.
- `interface_depths`: Interface depths in metres.

Output:
- `Matrix{Float64}`: Layered resistivity model.

Description:
- Builds a laterally uniform layered model on the 2D mesh.
"""
function build_mt2d_layered_model(
    mesh::MT2DMesh;
    layer_resistivities::AbstractVector{<:Real},
    interface_depths::AbstractVector{<:Real},
)
    length(layer_resistivities) == length(interface_depths) + 1 || error("layer_resistivities must contain one more value than interface_depths")

    n_z = length(mesh.z_cell_sizes)
    n_y = length(mesh.y_cell_sizes)
    ρ = fill(Float64(layer_resistivities[end]), n_z, n_y)
    ρ[1:mesh.n_air_cells, :] .= 1e9

    z_centers = mt2d_z_centers(mesh)
    z_interfaces = Float64.(interface_depths)
    ρlayers = Float64.(layer_resistivities)
    for iz in (mesh.n_air_cells + 1):n_z
        layer = searchsortedfirst(z_interfaces, z_centers[iz])
        ρ[iz, :] .= ρlayers[clamp(layer, 1, length(ρlayers))]
    end

    ρ
end

"""
    add_mt2d_rect!(resistivity, mesh; y_range, z_range, resistivity_value)

Inputs:
- `resistivity`: Existing 2D resistivity model.
- `mesh`: 2D MT mesh.
- `y_range`, `z_range`, `resistivity_value`: Rectangle geometry and value.

Output:
- `AbstractMatrix`: Updated resistivity model.

Description:
- Overwrites a rectangular region of a 2D model with a new resistivity value.
"""
function add_mt2d_rect!(
    resistivity::AbstractMatrix{<:Real},
    mesh::MT2DMesh;
    y_range::Tuple{<:Real, <:Real},
    z_range::Tuple{<:Real, <:Real},
    resistivity_value::Real,
)
    y1, y2 = Float64.(y_range)
    z1, z2 = Float64.(z_range)
    ρrect = Float64(resistivity_value)
    y_centers = mt2d_y_centers(mesh)
    z_centers = mt2d_z_centers(mesh)

    for iy in eachindex(y_centers), iz in (mesh.n_air_cells + 1):length(z_centers)
        if y1 <= y_centers[iy] <= y2 && z1 <= z_centers[iz] <= z2
            resistivity[iz, iy] = ρrect
        end
    end

    resistivity
end

"""
    build_mt2d_comemi_models(mesh)

Inputs:
- `mesh`: 2D MT mesh.

Output:
- `Vector`: Named tuples containing the benchmark names, labels, and resistivity models.

Description:
- Builds the three package COMEMI-style 2D benchmark models.
"""
function build_mt2d_comemi_models(mesh::MT2DMesh)
    case1 = build_mt2d_layered_model(mesh; layer_resistivities = [100.0, 500.0], interface_depths = [2000.0])
    add_mt2d_rect!(case1, mesh; y_range = (-1200.0, 1200.0), z_range = (200.0, 3500.0), resistivity_value = 5.0)

    case2 = build_mt2d_layered_model(mesh; layer_resistivities = [30.0, 100.0], interface_depths = [1500.0])
    add_mt2d_rect!(case2, mesh; y_range = (-5000.0, -500.0), z_range = (600.0, 2800.0), resistivity_value = 800.0)
    add_mt2d_rect!(case2, mesh; y_range = (2000.0, 6500.0), z_range = (1200.0, 4500.0), resistivity_value = 400.0)

    case3 = build_mt2d_layered_model(mesh; layer_resistivities = [80.0, 20.0, 300.0], interface_depths = [800.0, 3500.0])
    add_mt2d_rect!(case3, mesh; y_range = (-7000.0, -2000.0), z_range = (300.0, 1800.0), resistivity_value = 3.0)
    add_mt2d_rect!(case3, mesh; y_range = (1500.0, 6000.0), z_range = (2200.0, 6500.0), resistivity_value = 1000.0)

    [
        (name = "comemi2d_case1_dyke", label = "COMEMI2D 1", resistivity = case1),
        (name = "comemi2d_case2_resistive_blocks", label = "COMEMI2D 2", resistivity = case2),
        (name = "comemi2d_case3_mixed", label = "COMEMI2D 3", resistivity = case3),
    ]
end

"""
    mt2d_resistivity_at(mesh, resistivity, y, z)

Inputs:
- `mesh`: 2D MT mesh.
- `resistivity`: Cell resistivity model.
- `y`, `z`: Query coordinates in metres.

Output:
- `Float64`: Resistivity at the requested cell.

Description:
- Returns the resistivity of the cell containing the requested coordinate.
"""
function mt2d_resistivity_at(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real}, y::Real, z::Real)
    iy = findfirst(i -> mesh.y_nodes[i] <= y < mesh.y_nodes[i + 1], 1:length(mesh.y_cell_sizes))
    iz = findfirst(i -> mesh.z_nodes[i] <= z < mesh.z_nodes[i + 1], 1:length(mesh.z_cell_sizes))
    iy === nothing && error("requested y=$y m is outside mesh bounds")
    iz === nothing && error("requested z=$z m is outside mesh bounds")
    Float64(resistivity[iz, iy])
end

"""
    validate_mt2d_comemi_models(mesh, models)

Inputs:
- `mesh`: 2D MT mesh.
- `models`: COMEMI benchmark models.

Output:
- `Vector{String}`: Descriptions of the passed geometry checks.

Description:
- Verifies that the benchmark models contain the expected anomalies and host values.
"""
function validate_mt2d_comemi_models(mesh::MT2DMesh, models)
    lookup = Dict(model.name => model.resistivity for model in models)
    required = [
        "comemi2d_case1_dyke",
        "comemi2d_case2_resistive_blocks",
        "comemi2d_case3_mixed",
    ]
    for name in required
        haskey(lookup, name) || error("missing COMEMI benchmark model: $name")
    end

    checks = [
        ("comemi2d_case1_dyke", 0.0, 1000.0, 5.0, "conductive dyke core"),
        ("comemi2d_case1_dyke", 4000.0, 1000.0, 100.0, "upper host away from dyke"),
        ("comemi2d_case1_dyke", 4000.0, 3000.0, 500.0, "lower host away from dyke"),
        ("comemi2d_case2_resistive_blocks", -2000.0, 1200.0, 800.0, "left resistive block"),
        ("comemi2d_case2_resistive_blocks", 4000.0, 3000.0, 400.0, "right resistive block"),
        ("comemi2d_case2_resistive_blocks", 0.0, 1000.0, 30.0, "upper conductive host"),
        ("comemi2d_case2_resistive_blocks", 0.0, 3000.0, 100.0, "lower host"),
        ("comemi2d_case3_mixed", -4000.0, 1000.0, 3.0, "left conductive anomaly"),
        ("comemi2d_case3_mixed", 3000.0, 5000.0, 1000.0, "right deep resistive anomaly"),
        ("comemi2d_case3_mixed", 0.0, 500.0, 80.0, "shallow background"),
        ("comemi2d_case3_mixed", 0.0, 2000.0, 20.0, "middle background"),
        ("comemi2d_case3_mixed", 0.0, 8000.0, 300.0, "deep background"),
    ]

    results = String[]
    for (name, y, z, expected, label) in checks
        value = mt2d_resistivity_at(mesh, lookup[name], y, z)
        isapprox(value, expected; rtol = 1e-8, atol = 1e-8) || error("$name failed check '$label': got $value Ω·m, expected $expected Ω·m")
        push!(results, "$(name) | $(label) | rho=$(round(value; digits = 3)) ohm.m")
    end

    results
end

"""
    build_default_mt2d_mesh()

Inputs:
- None.

Output:
- `MT2DMesh`: Default COMEMI benchmark mesh.

Description:
- Builds the standard 2D benchmark mesh used by the package examples and tests.
"""
function build_default_mt2d_mesh()
    build_mt2d_mesh(
        frequencies = collect(10 .^ range(-2, 2, length = 7)),
        y_core_range = (-9000.0, 9000.0),
        y_core_cell = 400.0,
        y_padding = 12_000.0,
        pad_factor = 1.25,
        air_top = -15_000.0,
        air_cells = 6,
        ground_layers = vcat(fill(200.0, 8), fill(400.0, 10), fill(800.0, 10)),
        receiver_positions = collect(-8000.0:1600.0:8000.0),
    )
end

"""
    MakeMesh2D(; output_dir=...)

Inputs:
- `output_dir`: Directory where the benchmark models are written.

Output:
- Named tuple with `mesh` and `model_paths`.

Description:
- Builds the standard 2D benchmark mesh and writes the COMEMI-style models to disk.
"""
function MakeMesh2D(;
    output_dir::AbstractString = joinpath(dirname(@__DIR__), "Models"),
)
    mesh = build_default_mt2d_mesh()
    models = build_mt2d_comemi_models(mesh)
    paths = Dict{String, String}()
    filename_map = Dict(
        "comemi2d_case1_dyke" => "Comemi2D1.true",
        "comemi2d_case2_resistive_blocks" => "Comemi2D2.true",
        "comemi2d_case3_mixed" => "Comemi2D3.true",
    )

    for model in models
        filename = get(
            filename_map,
            model.name,
            replace(join(uppercasefirst.(split(model.name, "_")), ""), "1d" => "1D", "2d" => "2D", "3d" => "3D") * ".true",
        )
        paths[model.name] = write_model2d(joinpath(output_dir, filename), mesh, model.resistivity; title = model.label)
    end

    (mesh = mesh, model_paths = paths)
end

"""
    _write_vector_lines(io, values; per_line=12)

Inputs:
- Output stream, numeric values, and the number of values per line.

Output:
- Nothing.

Description:
- Writes a numeric vector to an open text stream using fixed-width scientific notation.
"""
function _write_vector_lines(io, values::AbstractVector{<:Real}; per_line::Int = 12)
    for first_index in 1:per_line:length(values)
        last_index = min(first_index + per_line - 1, length(values))
        println(io, join([@sprintf("%.8e", Float64(values[idx])) for idx in first_index:last_index], " "))
    end
end

"""
    write_model2d(path, mesh, resistivity; title="MTGeophysics.jl 2D profile model", use_loge=true)

Inputs:
- Output path, 2D mesh, resistivity model, and file-format options.

Output:
- `String`: Path to the written model file.

Description:
- Writes a 2D profile model file that can be reloaded by the package.
"""
function write_model2d(
    path::AbstractString,
    mesh::MT2DMesh,
    resistivity::AbstractMatrix{<:Real};
    title::AbstractString = "MTGeophysics.jl 2D profile model",
    use_loge::Bool = true,
)
    n_z = length(mesh.z_cell_sizes)
    n_y = length(mesh.y_cell_sizes)
    size(resistivity) == (n_z, n_y) || error("resistivity must be size ($(n_z), $(n_y))")

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# $(title)")
        println(io, "# NZA=$(mesh.n_air_cells)")
        println(io, "1 $(n_y) $(n_z) 0 $(use_loge ? "LOGE" : "LINEAR")")
        _write_vector_lines(io, [1.0])
        _write_vector_lines(io, mesh.y_cell_sizes)
        _write_vector_lines(io, mesh.z_cell_sizes)

        values = Float64[]
        sizehint!(values, n_y * n_z)
        for iz in 1:n_z, iy in 1:n_y
            ρ = Float64(resistivity[iz, iy])
            push!(values, use_loge ? log(ρ) : ρ)
        end
        _write_vector_lines(io, values)

        println(io, @sprintf("%.8e %.8e %.8e", 0.0, first(mesh.y_nodes), first(mesh.z_nodes)))
        println(io, "0.0")
    end

    String(path)
end

"""
    load_model2d(path)

Inputs:
- `path`: Path to a 2D model file.

Output:
- `ModelFile2D`: Parsed model file.

Description:
- Reads a 2D model file written by the package and reconstructs the stored metadata.
"""
function load_model2d(path::AbstractString)
    isfile(path) || error("model file not found: $path")
    lines = readlines(path)

    title = "MTGeophysics.jl 2D profile model"
    n_air_cells = 0
    dims_index = nothing
    for (line_index, line) in enumerate(lines)
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "#")
            title == "MTGeophysics.jl 2D profile model" && (title = strip(replace(stripped, "#" => "")))
            if occursin("NZA=", stripped)
                parts = split(stripped, "NZA=")
                length(parts) > 1 && (n_air_cells = something(tryparse(Int, strip(parts[2])), 0))
            end
            continue
        end
        tokens = split(stripped)
        if length(tokens) >= 5 && all(token -> occursin(r"^-?\d+$", token), tokens[1:4])
            dims_index = line_index
            break
        end
    end
    dims_index === nothing && error("could not locate model dimension line in $path")

    dims_tokens = split(strip(lines[dims_index]))
    n_x = parse(Int, dims_tokens[1])
    n_y = parse(Int, dims_tokens[2])
    n_z = parse(Int, dims_tokens[3])
    format = uppercase(dims_tokens[5])
    n_x == 1 || error("this 2D profile file expects nx=1, got $n_x")

    numeric_tokens = String[]
    for line in lines[(dims_index + 1):end]
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        append!(numeric_tokens, split(stripped))
    end
    values = parse.(Float64, numeric_tokens)
    required_values = n_x + n_y + n_z + n_x * n_y * n_z + 4
    length(values) >= required_values || error("model file is incomplete: expected at least $required_values numeric values, got $(length(values))")

    index = 1
    x_cell_sizes = values[index:(index + n_x - 1)]
    index += n_x
    y_cell_sizes = values[index:(index + n_y - 1)]
    index += n_y
    z_cell_sizes = values[index:(index + n_z - 1)]
    index += n_z

    block_values = values[index:(index + n_x * n_y * n_z - 1)]
    index += n_x * n_y * n_z
    ρvalues = format == "LOGE" ? exp.(block_values) : block_values

    resistivity = Matrix{Float64}(undef, n_z, n_y)
    value_index = 1
    for iz in 1:n_z, iy in 1:n_y
        resistivity[iz, iy] = ρvalues[value_index]
        value_index += 1
    end

    origin = Float64[values[index], values[index + 1], values[index + 2]]
    index += 3
    rotation = values[index]

    ModelFile2D(
        title = title,
        x_cell_sizes = Float64.(x_cell_sizes),
        y_cell_sizes = Float64.(y_cell_sizes),
        z_cell_sizes = Float64.(z_cell_sizes),
        resistivity = resistivity,
        n_air_cells = n_air_cells,
        origin = origin,
        rotation = rotation,
        format = format,
        path = String(path),
    )
end

"""
    build_mesh_from_model2d(model; frequencies, receiver_positions)

Inputs:
- `model`: Parsed 2D model file.
- `frequencies`: Survey frequencies in hertz.
- `receiver_positions`: Receiver offsets in metres.

Output:
- `MT2DMesh`: Solver mesh matching the stored model.

Description:
- Reconstructs a 2D solver mesh from a saved model file and survey definition.
"""
function build_mesh_from_model2d(
    model::ModelFile2D;
    frequencies::AbstractVector{<:Real},
    receiver_positions::AbstractVector{<:Real},
)
    y_nodes = model.origin[2] .+ vcat(0.0, cumsum(model.y_cell_sizes))
    z_nodes = model.origin[3] .+ vcat(0.0, cumsum(model.z_cell_sizes))
    MT2DMesh(
        y_nodes = Float64.(y_nodes),
        z_nodes = Float64.(z_nodes),
        y_cell_sizes = Float64.(model.y_cell_sizes),
        z_cell_sizes = Float64.(model.z_cell_sizes),
        receiver_positions = Float64.(receiver_positions),
        frequencies = Float64.(frequencies),
        n_air_cells = model.n_air_cells,
    )
end


# 2D MT mesh and model
# Author: @pankajkmishra
# Profile mesh geometry, skin-depth vertical core, model builders, and model file I/O

using LinearAlgebra
using Printf

"""
    MT2DMesh

Profile mesh and survey: y (along the profile) and z (down) nodes and cell sizes, air
rows first, receiver positions and frequencies. `topo_air` holds the topographic air
cells at the top of each earth column (empty = flat); every receiver sits on the ground
surface of its column; next to a topographic step TM Ey is averaged over the columns
within half of `dipole_length` metres each side.
"""
Base.@kwdef struct MT2DMesh
    y_nodes::Vector{Float64}
    z_nodes::Vector{Float64}
    y_cell_sizes::Vector{Float64}
    z_cell_sizes::Vector{Float64}
    receiver_positions::Vector{Float64}
    frequencies::Vector{Float64}
    n_air_cells::Int
    air_resistivity::Float64 = 1e9
    topo_air::Vector{Int} = Int[]
    dipole_length::Float64 = 100.0
end

"""
    ModelFile2D

A 2D model file on disk: cell sizes, resistivity `(nz, ny)`, air row count and origin.
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

#---------- topography ----------
#
# as in 3D (MakeMesh3D, load_ws3d_model): model files hold topographic air as 1e17 ohm m
# and any cell above 1e15 ohm m reads as air; the cov mask gives air 0 and water 9

const MT2D_AIR_TAG = 1e17
const MT2D_AIR_THRESHOLD = 1e15
const MT2D_MASK_WATER = 9

# topographic air cells per earth column, 0 everywhere when flat
mt2d_topo_air(mesh::MT2DMesh) = isempty(mesh.topo_air) ? zeros(Int, length(mesh.y_cell_sizes)) : mesh.topo_air

"""
    mt2d_air_mask(mesh) -> BitMatrix

Air cells of the full mesh `(nz, ny)`: the air layers plus the topographic air.
"""
function mt2d_air_mask(mesh::MT2DMesh)
    air = falses(length(mesh.z_cell_sizes), length(mesh.y_cell_sizes))
    for (iy, k) in enumerate(mt2d_topo_air(mesh))
        air[1:mesh.n_air_cells+k, iy] .= true
    end
    air
end

# column of each receiver, and the depth of its ground surface below the model top
mt2d_receiver_columns(mesh::MT2DMesh) = [searchsortedlast(mesh.y_nodes, y) for y in mesh.receiver_positions]
mt2d_receiver_depths(mesh::MT2DMesh) =
    [mesh.z_nodes[mesh.n_air_cells+1+mt2d_topo_air(mesh)[iy]] - mesh.z_nodes[mesh.n_air_cells+1]
     for iy in mt2d_receiver_columns(mesh)]

# leading tagged cells of each column of an earth model; tagged cells under the ground are an error
function _mt2d_model_topo_air(ρ::AbstractMatrix{<:Real})
    tagged = ρ .> MT2D_AIR_THRESHOLD
    topo = [something(findfirst(!, tagged[:, iy]), size(ρ, 1) + 1) - 1 for iy in axes(ρ, 2)]
    all(<(size(ρ, 1)), topo) || throw(ArgumentError("a model column holds air only"))
    for iy in axes(ρ, 2)
        stray = findfirst(tagged[topo[iy]+1:end, iy])
        stray === nothing || throw(ArgumentError("air-tagged cell under the ground in column $iy, row $(topo[iy] + stray)"))
    end
    topo
end

"""
    mt2d_station_offsets(mesh, data) -> Vector{NamedTuple}

Snap of each station to the ground surface of its column: `site`, `z` (data depth
below the model top), `surface` (the mesh's ground depth there), `offset = z - surface`
and `tolerance`, half the thickness of the first earth cell of that column.
"""
function mt2d_station_offsets(mesh::MT2DMesh, data)
    columns, depths = mt2d_receiver_columns(mesh), mt2d_receiver_depths(mesh)
    topo = mt2d_topo_air(mesh)
    map(eachindex(columns)) do i
        dz = mesh.z_cell_sizes[mesh.n_air_cells+topo[columns[i]]+1]
        (site = data.site_names[i], z = data.z_positions[i], surface = depths[i],
         offset = data.z_positions[i] - depths[i], tolerance = dz / 2)
    end
end

mt2d_y_centers(mesh::MT2DMesh) = 0.5 .* (mesh.y_nodes[1:end-1] .+ mesh.y_nodes[2:end])

mt2d_z_centers(mesh::MT2DMesh) = 0.5 .* (mesh.z_nodes[1:end-1] .+ mesh.z_nodes[2:end])

# δ = sqrt(2ρ/(ωμ₀)) ≈ 503 sqrt(ρ/f) metres
mt2d_skin_depth(resistivity::Real, frequency::Real) = sqrt(2 * resistivity / (2π * frequency * μ₀_2D))

"""
    mt2d_skin_depth_layers(frequencies; background_resistivity=100.0, z_core_cell=nothing,
                           z_core_skin_depths=1.0, z_bottom_skin_depths=4.0,
                           max_core_layers=80, pad_factor=1.25) -> thicknesses

Earth layers from the surface down: a uniform core of `z_core_cell` (default
`max(δ_min/3, δ_max/max_core_layers)`) down to `z_core_skin_depths` skin depths of the
lowest frequency, then layers growing by `pad_factor` down to `z_bottom_skin_depths`.
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
    mt2d_geometric_layers(frequencies; background_resistivity=100.0, first_layer_div=5.0,
                          vertical_factor=1.1, depth_mult=4.0) -> thicknesses

Earth layers from the surface down, as in `MakeMesh3D`: the first layer is the skin
depth of the highest frequency over `first_layer_div`, each next layer is
`vertical_factor` thicker, down to `depth_mult` skin depths of the lowest frequency.
"""
function mt2d_geometric_layers(frequencies::AbstractVector{<:Real}; background_resistivity::Real = 100.0,
                               first_layer_div::Real = 5.0, vertical_factor::Real = 1.1, depth_mult::Real = 4.0)
    first_layer_div > 0 && vertical_factor >= 1 && depth_mult > 0 ||
        throw(ArgumentError("need first_layer_div > 0, vertical_factor >= 1, depth_mult > 0"))
    layers = [mt2d_skin_depth(background_resistivity, maximum(frequencies)) / first_layer_div]
    bottom = depth_mult * mt2d_skin_depth(background_resistivity, minimum(frequencies))
    while sum(layers) < bottom
        push!(layers, layers[end] * vertical_factor)
    end
    layers
end

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
    BuildMesh2D(; frequencies, y_core_range=(-6000.0, 6000.0), y_core_cell=300.0,
                y_padding=9000.0, pad_factor=1.25, air_top=-12000.0, air_cells=8,
                ground_layers=nothing, receiver_stride=2, receiver_positions=nothing,
                background_resistivity, z_core_cell, z_core_skin_depths,
                z_bottom_skin_depths, max_core_layers) -> MT2DMesh

Padded profile mesh: a uniform core of `y_core_cell` over `y_core_range`, padding
growing by `pad_factor` over `y_padding` on each side, `air_cells` uniform air layers up
to `air_top`, and `ground_layers` below the surface (default `mt2d_skin_depth_layers`).
Receivers default to every `receiver_stride`-th core cell centre.
"""
BuildMesh2D(; kwargs...) = build_mt2d_mesh(; kwargs...)

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

# homogeneous halfspace under 1e9 ohm m air
build_mt2d_halfspace_model(mesh::MT2DMesh; background_resistivity::Real = 100.0) =
    build_mt2d_block_model(mesh; background_resistivity = background_resistivity, blocks = NamedTuple[])

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

# the three COMEMI-style benchmark models of helpers/benchmarks_2D.jl
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

# small benchmark mesh used by the tests
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

# fixed-width scientific notation, per_line values per line
function _write_vector_lines(io, values::AbstractVector{<:Real}; per_line::Int = 12)
    for first_index in 1:per_line:length(values)
        last_index = min(first_index + per_line - 1, length(values))
        println(io, join([@sprintf("%.8e", Float64(values[idx])) for idx in first_index:last_index], " "))
    end
end

# the older MTGeophysics model layout with air rows; ReadModel2D still reads it, the tests write it
function _write_model2d_legacy(
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

function _load_model2d_legacy(path::AbstractString)
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

# solver mesh of an older-layout model, air rows included
function _mesh_from_legacy_model2d(
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

#---------- ModEM-layout model files ----------
#
# layout, earth cells only, after one '#' description line as in
# ModEM 3D model files:
#   ny nz LOGE            (or LINEAR)
#   ny cell widths (m)
#   nz layer thicknesses (m), top to bottom
#   0
#   nz rows of ny values, top row first, ln(ρ) for LOGE, ρ for LINEAR
# the grid is centred on y = 0, the same origin as the data file's local y

"""
    ReadModel2D(path) -> ModelFile2D

Read a ModEM-layout 2D model, earth cells only, centred on y = 0. The older
MTGeophysics layout with air rows is also accepted; its air rows are dropped.
"""
function ReadModel2D(path::AbstractString)
    isfile(path) || error("model file not found: $path")
    lines = filter(l -> !isempty(strip(l)) && !startswith(strip(l), "#"), readlines(path))
    isempty(lines) && error("$path: empty model file")
    head = split(strip(lines[1]))
    if length(head) >= 5 && all(t -> occursin(r"^-?\d+$", t), head[1:4])
        legacy = _load_model2d_legacy(path)
        na = legacy.n_air_cells
        return ModelFile2D(title = legacy.title, x_cell_sizes = [1.0], y_cell_sizes = legacy.y_cell_sizes,
                           z_cell_sizes = legacy.z_cell_sizes[na+1:end], resistivity = legacy.resistivity[na+1:end, :],
                           n_air_cells = 0, origin = [0.0, legacy.origin[2], 0.0], rotation = legacy.rotation,
                           format = legacy.format, path = String(path))
    end
    length(head) >= 2 || error("$path: first line must be 'ny nz LOGE' or 'ny nz LINEAR'")
    ny, nz = parse(Int, head[1]), parse(Int, head[2])
    format = occursin("LOGE", uppercase(strip(lines[1]))) ? "LOGE" : "LINEAR"
    tokens = parse.(Float64, reduce(vcat, split.(strip.(lines[2:end]))))
    length(tokens) == ny + nz + 1 + ny * nz ||
        error("$path: expected $(ny + nz + 1 + ny * nz) numbers after the header, got $(length(tokens))")
    dy = tokens[1:ny]
    dz = tokens[ny+1:ny+nz]
    values = permutedims(reshape(tokens[ny+nz+2:end], ny, nz))
    ρ = format == "LOGE" ? exp.(values) : values
    all(>(0), dy) && all(>(0), dz) && all(x -> isfinite(x) && x > 0, ρ) ||
        error("$path: cell sizes and resistivities must be positive")
    ModelFile2D(title = "", x_cell_sizes = [1.0], y_cell_sizes = dy, z_cell_sizes = dz, resistivity = ρ,
                n_air_cells = 0, origin = [0.0, -sum(dy) / 2, 0.0], rotation = 0.0, format = format,
                path = String(path))
end

"""
    WriteModel2D(path, y_cell_sizes, z_cell_sizes, ρ; loge=true) -> path
    WriteModel2D(path, mesh, ρ; loge=true) -> path

Write a ModEM-layout 2D model. With a mesh, `ρ` holds every mesh row and the air
rows are left out.
"""
function WriteModel2D(path::AbstractString, dy::AbstractVector{<:Real}, dz::AbstractVector{<:Real},
                      ρ::AbstractMatrix{<:Real}; loge::Bool = true)
    size(ρ) == (length(dz), length(dy)) || throw(DimensionMismatch("model must be (nz, ny) = ($(length(dz)), $(length(dy)))"))
    all(x -> isfinite(x) && x > 0, ρ) || throw(ArgumentError("resistivity must be positive and finite"))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# 2D MT model written by MTGeophysics.jl in ModEM format")
        @printf(io, "%5d%5d %s\n", length(dy), length(dz), loge ? "LOGE" : "LINEAR")
        _write_model_rows(io, "%12.3f", dy)
        _write_model_rows(io, "%12.3f", dz)
        println(io, "0")
        for iz in axes(ρ, 1)
            println(io)
            _write_model_rows(io, "%13.5E", loge ? log.(ρ[iz, :]) : ρ[iz, :])
        end
    end
    String(path)
end

# ten fixed-width values per line
function _write_model_rows(io, fmt::String, values)
    f = Printf.Format(fmt)
    for row in Iterators.partition(values, 10)
        foreach(v -> Printf.format(io, f, v), row)
        println(io)
    end
end

function WriteModel2D(path::AbstractString, mesh::MT2DMesh, ρ::AbstractMatrix{<:Real}; loge::Bool = true)
    # the layout has no origin, readers centre the grid on y = 0
    abs(mesh.y_nodes[1] + mesh.y_nodes[end]) <= 1e-6 * (mesh.y_nodes[end] - mesh.y_nodes[1]) ||
        throw(ArgumentError("mesh must be centred on y = 0, it spans $(mesh.y_nodes[1]) to $(mesh.y_nodes[end]) m"))
    na = mesh.n_air_cells
    earth = Matrix{Float64}(ρ[na+1:end, :])
    earth[mt2d_air_mask(mesh)[na+1:end, :]] .= MT2D_AIR_TAG
    WriteModel2D(path, mesh.y_cell_sizes, mesh.z_cell_sizes[na+1:end], earth; loge)
end

#---------- mesh from model, data and fwd.ctrl ----------

"""
    mt2d_air_layers(n, thickness, growth) -> Vector{Float64}

Air layer thicknesses, top to bottom: `n` layers adding up to `thickness`, each
`growth` times thicker than the one below it.
"""
function mt2d_air_layers(n::Integer, thickness::Real, growth::Real)
    n > 0 && thickness > 0 && growth >= 1 || throw(ArgumentError("need n > 0, thickness > 0, growth >= 1"))
    t1 = growth == 1 ? thickness / n : thickness * (growth - 1) / (growth^n - 1)
    reverse(t1 .* growth .^ (0:n-1))
end

"""
    Mesh2DFromInputs(model::ModelFile2D, data::DataFile2D, fwd::FwdCtrl2D; warn=true) -> (mesh, ρ)

Mesh and full resistivity (air rows first) for a ModEM-layout model, the survey in
`data`, and the air in `fwd`. The model is centred on the data's local y = 0.
Topographic air (cells above 1e15 ohm m) takes the `fwd` air resistivity. Each station
sits on the ground of its column; a data depth Z more than half a surface cell away
from it is warned about (see `mt2d_station_offsets`).
"""
function Mesh2DFromInputs(model::ModelFile2D, data, fwd::FwdCtrl2D; warn::Bool = true)
    model.n_air_cells == 0 || throw(ArgumentError("model must hold earth cells only, read it with ReadModel2D"))
    topo = _mt2d_model_topo_air(model.resistivity)
    air = mt2d_air_layers(fwd.air_layers, fwd.air_thickness, fwd.air_growth)
    dz = vcat(air, model.z_cell_sizes)
    y_nodes = model.origin[2] .+ vcat(0.0, cumsum(model.y_cell_sizes))
    z_nodes = vcat(0.0, cumsum(dz)) .- sum(air)
    all(y -> y_nodes[1] <= y < y_nodes[end], data.receivers) ||
        throw(ArgumentError("stations lie outside the model, which spans $(y_nodes[1]) to $(y_nodes[end]) m"))
    mesh = MT2DMesh(y_nodes = y_nodes, z_nodes = z_nodes, y_cell_sizes = Float64.(model.y_cell_sizes),
                    z_cell_sizes = dz, receiver_positions = Float64.(data.receivers),
                    frequencies = Float64.(data.frequencies), n_air_cells = fwd.air_layers,
                    air_resistivity = fwd.air_resistivity, topo_air = any(>(0), topo) ? topo : Int[],
                    dipole_length = fwd.dipole_length)
    ρ = vcat(fill(fwd.air_resistivity, fwd.air_layers, length(model.y_cell_sizes)), model.resistivity)
    ρ[mt2d_air_mask(mesh)] .= fwd.air_resistivity
    if warn
        off = filter(o -> abs(o.offset) > o.tolerance, mt2d_station_offsets(mesh, data))
        isempty(off) || @warn "$(length(off)) station(s) moved to the ground of their column by more than half a surface cell" *
            join([@sprintf("\n  %s: Z %.1f m, ground %.1f m", o.site, o.z, o.surface) for o in off])
    end
    mesh, ρ
end

using CairoMakie
using LinearAlgebra
using Printf
using Random
using SparseArrays

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
    MT2DResponse

Inputs:
- Frequencies, periods, receiver positions, and TE/TM impedance-derived responses.

Output:
- `MT2DResponse`: Container for a 2D MT forward response.

Description:
- Stores the TE and TM apparent resistivity, phase, and impedance responses for a profile survey.
"""
Base.@kwdef struct MT2DResponse
    frequencies::Vector{Float64}
    periods::Vector{Float64}
    receivers::Vector{Float64}
    rho_xy::Matrix{Float64}
    phase_xy::Matrix{Float64}
    z_xy::Matrix{ComplexF64}
    rho_yx::Matrix{Float64}
    phase_yx::Matrix{Float64}
    z_yx::Matrix{ComplexF64}
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

"""
    DataFile2D

Inputs:
- Survey metadata, site geometry, impedance tensors, errors, and derived resistivity/phase values.

Output:
- `DataFile2D`: Parsed 2D MT data container.

Description:
- Stores the impedance observations and uncertainties used by the 2D file-based workflows.
"""
Base.@kwdef struct DataFile2D
    title::String
    periods::Vector{Float64}
    frequencies::Vector{Float64}
    site_names::Vector{String}
    receivers::Vector{Float64}
    x_positions::Vector{Float64}
    z_positions::Vector{Float64}
    z_xy::Matrix{ComplexF64}
    z_xy_error::Matrix{Float64}
    z_yx::Matrix{ComplexF64}
    z_yx_error::Matrix{Float64}
    z_xx::Matrix{ComplexF64}
    z_xx_error::Matrix{Float64}
    z_yy::Matrix{ComplexF64}
    z_yy_error::Matrix{Float64}
    rho_xy::Matrix{Float64}
    phase_xy::Matrix{Float64}
    rho_yx::Matrix{Float64}
    phase_yx::Matrix{Float64}
    path::String = ""
end

"""
    FitSummary2D

Inputs:
- Chi-square, RMS, and sample count values.

Output:
- `FitSummary2D`: Simple data-fit summary.

Description:
- Stores the scalar fit metrics used by the 2D forward and inversion workflows.
"""
Base.@kwdef struct FitSummary2D
    chi2::Float64
    rms::Float64
    count::Int
end

"""
    TensorMesh2D

Inputs:
- Cell sizes, tensor-mesh dimensions, origin, conductivity, and sparse operators.

Output:
- `TensorMesh2D`: Solver-ready tensor mesh.

Description:
- Stores the sparse-operator form of the 2D tensor mesh used by the finite-volume solver.
"""
mutable struct TensorMesh2D
    y_lengths::Vector{Float64}
    z_lengths::Vector{Float64}
    grid_size::Vector{Int}
    origin::Vector{Float64}
    conductivity::Vector{Float64}
    face::SparseMatrixCSC{Float64, Int}
    gradient::SparseMatrixCSC{Float64, Int}
    average_cell_to_node::SparseMatrixCSC{Float64, Int}
    average_cell_to_face::SparseMatrixCSC{Float64, Int}
    setup::Bool
end

"""
    CoeffMat

Inputs:
- Real and imaginary interior/interior and interior/exterior sparse blocks.

Output:
- `CoeffMat`: Sparse coefficient blocks for a TE or TM system.

Description:
- Stores the block matrices used to assemble and solve the 2D TE and TM linear systems.
"""
mutable struct CoeffMat{T<:Float64}
    real_ii::SparseMatrixCSC{T, Int}
    imag_ii::SparseMatrixCSC{T, Int}
    real_io::SparseMatrixCSC{T, Int}
    imag_io::SparseMatrixCSC{T, Int}
end

const μ₀_2D = 4π * 1e-7
const MU0_2D = μ₀_2D

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
    build_mt2d_mesh(; frequencies=..., y_core_range=(-6000.0, 6000.0), y_core_cell=300.0, y_padding=9000.0, pad_factor=1.25, air_top=-12000.0, air_cells=8, ground_layers=..., receiver_stride=2, receiver_positions=nothing)

Inputs:
- Frequency axis and horizontal/vertical mesh controls.

Output:
- `MT2DMesh`: Survey mesh and receiver geometry.

Description:
- Builds the padded 2D MT profile mesh used by the forward and inversion workflows.
"""
function build_mt2d_mesh(;
    frequencies::AbstractVector{<:Real} = collect(10 .^ range(-2, 2, length = 10)),
    y_core_range::Tuple{<:Real, <:Real} = (-6000.0, 6000.0),
    y_core_cell::Real = 300.0,
    y_padding::Real = 9000.0,
    pad_factor::Real = 1.25,
    air_top::Real = -12000.0,
    air_cells::Integer = 8,
    ground_layers::AbstractVector{<:Real} = vcat(fill(100.0, 8), fill(250.0, 10), fill(500.0, 10)),
    receiver_stride::Integer = 2,
    receiver_positions::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    y_core_cell > 0 || error("y_core_cell must be positive")
    y_padding > 0 || error("y_padding must be positive")
    pad_factor > 1 || error("pad_factor must be greater than 1")
    air_cells > 0 || error("air_cells must be positive")
    receiver_stride > 0 || error("receiver_stride must be positive")
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
    spunit(n)

Inputs:
- `n`: Matrix dimension.

Output:
- Sparse identity matrix.

Description:
- Builds an `n × n` sparse identity matrix.
"""
spunit(n::Integer) = sparse(1.0I, n, n)

"""
    spdiag((x1, x2), (d1, d2), m, n)

Inputs:
- Two diagonal vectors, their offsets, and matrix dimensions.

Output:
- Sparse matrix.

Description:
- Convenience helper for constructing a sparse two-diagonal matrix.
"""
spdiag((x1, x2), (d1, d2), m, n) = spdiagm(m, n, d1 => x1, d2 => x2)

"""
    ddx(n)

Inputs:
- `n`: Number of cells.

Output:
- Sparse first-difference operator.

Description:
- Builds the nodal first-difference operator used by the tensor-mesh discretization.
"""
ddx(n::Integer) = spdiag((-ones(n), ones(n)), (0, 1), n, n + 1)

"""
    av(n)

Inputs:
- `n`: Number of cells.

Output:
- Sparse averaging operator.

Description:
- Builds the simple one-dimensional node-to-cell averaging operator.
"""
av(n::Integer) = spdiag((0.5 .* ones(n), 0.5 .* ones(n)), (0, 1), n, n + 1)

"""
    avcn(n)

Inputs:
- `n`: Number of cells.

Output:
- Sparse averaging operator with boundary preservation.

Description:
- Builds the cell-to-node averaging operator used by the tensor-mesh assembly.
"""
function avcn(n::Integer)
    A = spdiagm(n + 1, n, -1 => 0.5 .* ones(n), 0 => 0.5 .* ones(n))
    A[1, 1] = 1.0
    A[end, end] = 1.0
    A
end

"""
    sdiag(values)

Inputs:
- `values`: Diagonal entries.

Output:
- Sparse diagonal matrix.

Description:
- Builds a sparse diagonal matrix from a vector.
"""
sdiag(values::AbstractVector) = spdiagm(0 => values)

"""
    mesh_geo_face_2d(d1, d2)

Inputs:
- `d1`, `d2`: Cell sizes along the two tensor-mesh axes.

Output:
- Sparse diagonal face-geometry matrix.

Description:
- Builds the face geometry matrix used in the 2D tensor-mesh discretization.
"""
mesh_geo_face_2d(d1::Vector, d2::Vector) = kron(sdiag(d2), sdiag(d1))

"""
    mesh_geo_edge_inv_2d(d1, d2)

Inputs:
- `d1`, `d2`: Cell sizes along the two tensor-mesh axes.

Output:
- Sparse inverse edge-geometry matrix.

Description:
- Builds the inverse edge geometry matrix for the 2D tensor-mesh discretization.
"""
function mesh_geo_edge_inv_2d(d1::Vector, d2::Vector)
    n1 = length(d1)
    n2 = length(d2)
    left = kron(spunit(n2 + 1), sdiag(1.0 ./ d1))
    right = kron(sdiag(1.0 ./ d2), spunit(n1 + 1))
    blockdiag(left, right)
end

"""
    nodal_gradient_2d(d1, d2)

Inputs:
- `d1`, `d2`: Cell sizes along the two tensor-mesh axes.

Output:
- Sparse nodal-gradient matrix.

Description:
- Builds the 2D tensor-mesh nodal-gradient operator.
"""
function nodal_gradient_2d(d1::Vector, d2::Vector)
    n1 = length(d1)
    n2 = length(d2)
    g1 = kron(spunit(n2 + 1), ddx(n1))
    g2 = kron(ddx(n2), spunit(n1 + 1))
    mesh_geo_edge_inv_2d(d1, d2) * [g1; g2]
end

"""
    average_cell_to_node_2d(grid_size)

Inputs:
- `grid_size`: Two-element vector `[n_y, n_z]`.

Output:
- Sparse averaging matrix.

Description:
- Builds the 2D tensor-mesh cell-to-node averaging operator.
"""
average_cell_to_node_2d(grid_size::Vector{Int}) = kron(avcn(grid_size[2]), avcn(grid_size[1]))

"""
    average_cell_to_face_2d(grid_size)

Inputs:
- `grid_size`: Two-element vector `[n_y, n_z]`.

Output:
- Sparse averaging matrix.

Description:
- Builds the 2D tensor-mesh cell-to-face averaging operator.
"""
function average_cell_to_face_2d(grid_size::Vector{Int})
    face_y = kron(spunit(grid_size[2]), avcn(grid_size[1]))
    face_z = kron(avcn(grid_size[2]), spunit(grid_size[1]))
    [face_z; face_y]
end

"""
    TensorMesh2D(y_lengths, z_lengths; origin=[0.0, 0.0], conductivity=nothing)

Inputs:
- Cell sizes, origin, and optional conductivity vector.

Output:
- `TensorMesh2D`: Unassembled solver mesh.

Description:
- Builds the tensor-mesh container used by the 2D finite-volume solver.
"""
function TensorMesh2D(
    y_lengths::Vector{Float64},
    z_lengths::Vector{Float64};
    origin = [0.0, 0.0],
    conductivity = nothing,
)
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    σ = conductivity === nothing ? fill(1 / 100.0, n_y * n_z) : conductivity
    TensorMesh2D(
        y_lengths,
        z_lengths,
        [n_y, n_z],
        Float64.(origin),
        σ,
        spzeros(Float64, n_y * n_z, n_y * n_z),
        spzeros(Float64, (n_y + 1 + n_z) * (n_y + 1), (n_y + 1) * (n_z + 1)),
        spzeros(Float64, (n_y + 1) * (n_z + 1), n_y * n_z),
        spzeros(Float64, (n_y * n_z + (n_y + 1) * n_z), n_y * n_z),
        false,
    )
end

"""
    setup_tensor_mesh_2d!(mesh)

Inputs:
- `mesh`: Tensor mesh to assemble.

Output:
- `TensorMesh2D`: Updated mesh with sparse operators.

Description:
- Assembles the sparse geometry and averaging operators required by the 2D solver.
"""
function setup_tensor_mesh_2d!(mesh::TensorMesh2D)
    mesh.face = mesh_geo_face_2d(mesh.y_lengths, mesh.z_lengths)
    mesh.gradient = nodal_gradient_2d(mesh.y_lengths, mesh.z_lengths)
    mesh.average_cell_to_node = average_cell_to_node_2d(mesh.grid_size)
    mesh.average_cell_to_face = average_cell_to_face_2d(mesh.grid_size)
    mesh.setup = true
    mesh
end

"""
    mt1d_boundary_field(frequency, conductivity, z_nodes; return_magnetic=false)

Inputs:
- `frequency`: Frequency in hertz.
- `conductivity`: Layer conductivities.
- `z_nodes`: Layer nodes in metres.
- `return_magnetic`: Whether to also return the magnetic field.

Output:
- Electric field vector, or electric and magnetic field vectors.

Description:
- Solves the 1D layered boundary problem used to impose TE and TM boundary conditions.
"""
function mt1d_boundary_field(
    frequency::Float64,
    conductivity::Vector{Float64},
    z_nodes::Vector{Float64};
    return_magnetic::Bool = false,
)
    length(conductivity) == length(z_nodes) - 1 || error("conductivity layers and z-nodes do not match")

    Etop = 1.0 + 0im
    ε₀ = 8.85e-12
    ω = 2π * frequency
    σext = vcat(conductivity, conductivity[end])
    n_layers = length(z_nodes)
    Δz = diff(z_nodes)

    k = sqrt(μ₀_2D * ε₀ * ω^2 - μ₀_2D * σext[end] * ω * 1im)
    Zs = ω * μ₀_2D / k
    for i in (n_layers - 1):-1:1
        k = sqrt(μ₀_2D * ε₀ * ω^2 - μ₀_2D * σext[i] * ω * 1im)
        Zi = ω * μ₀_2D / k
        q = tanh(k * Δz[i] * 1im)
        Zs = Zi * (Zs + Zi * q) / (Zi + Zs * q)
    end

    layers = zeros(ComplexF64, 2, n_layers)
    layers[1, 1] = 0.5 * Etop * (1 - ω * μ₀_2D / (Zs * k))
    layers[2, 1] = 0.5 * Etop * (1 + ω * μ₀_2D / (Zs * k))
    kall = sqrt.(μ₀_2D * ε₀ * ω^2 .- μ₀_2D .* σext .* ω .* 1im)

    for i in 1:(n_layers - 1)
        ratio = kall[i] / kall[i + 1]
        T = 0.5 .* [1 + ratio 1 - ratio; 1 - ratio 1 + ratio]
        P = [
            exp(kall[i] * Δz[i] * 1im) 0
            0 exp(-kall[i] * Δz[i] * 1im)
        ]
        layers[:, i + 1] = T * P * layers[:, i]
        downstream = abs(layers[1, i + 1] + layers[2, i + 1])
        upstream = abs(layers[1, i] + layers[2, i])
        if downstream > upstream || isnan(downstream)
            layers[:, i + 1:end] .= 0.0
            break
        end
    end

    E = transpose(sum(layers, dims = 1))
    if !return_magnetic
        return E
    end

    Hlayers = [
        layers[1:1, :] * sparse(Diagonal(-kall)) / (ω * μ₀_2D * 1im)
        layers[2:2, :] * sparse(Diagonal(kall)) / (ω * μ₀_2D * 1im)
    ]
    H = transpose(sum(Hlayers, dims = 1))
    E, H
end

"""
    get_boundary_index(n_y, n_z)

Inputs:
- `n_y`, `n_z`: Number of cells in the horizontal and vertical directions.

Output:
- Tuple of interior and exterior node indices.

Description:
- Splits the tensor-mesh nodes into interior and boundary sets.
"""
function get_boundary_index(n_y::Int, n_z::Int)
    n_nodes = (n_y + 1) * (n_z + 1)
    index_grid = reshape(collect(1:n_nodes), n_y + 1, n_z + 1)'
    inside = reshape((index_grid[2:end-1, 2:end-1])', (n_y - 1) * (n_z - 1))
    top = reshape(index_grid[1, :], n_y + 1)
    left = index_grid[2:end, 1]
    right = index_grid[2:end, end]
    bottom = reshape(index_grid[end, 2:end-1], n_y - 1)
    outside = [top; left; right; bottom]
    inside, outside
end

"""
    get_boundary_mt2d_te(frequency, y_lengths, z_lengths, conductivity)

Inputs:
- Frequency, tensor-mesh cell sizes, and conductivity vector.

Output:
- `Vector{ComplexF64}`: TE boundary field values.

Description:
- Builds the TE boundary conditions from stitched 1D edge solutions.
"""
function get_boundary_mt2d_te(
    frequency::Float64,
    y_lengths::Vector{Float64},
    z_lengths::Vector{Float64},
    conductivity::Vector{Float64},
)
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    z_nodes = [0.0; cumsum(z_lengths)]
    σ2d = reshape(conductivity, n_y, n_z)'
    boundary = zeros(ComplexF64, 2 * (n_y + n_z))

    boundary[1:n_y+1] .= 1.0 + 0.0im
    σ1d = σ2d[:, 1]
    E = mt1d_boundary_field(frequency, σ1d, z_nodes)
    E ./= E[1]
    boundary[n_y+2:n_y+n_z+1] = E[2:end]

    σ1d = σ2d[:, end]
    E = mt1d_boundary_field(frequency, σ1d, z_nodes)
    E ./= E[1]
    boundary[n_y+n_z+2:n_y+2*n_z+1] = E[2:end]

    for iy in 2:n_y
        σ1d = (σ2d[:, iy - 1] * y_lengths[iy - 1] + σ2d[:, iy] * y_lengths[iy]) / (y_lengths[iy - 1] + y_lengths[iy])
        E = mt1d_boundary_field(frequency, σ1d, z_nodes)
        boundary[n_y + 2 * n_z + iy] = E[end] / E[1]
    end

    boundary
end

"""
    get_boundary_mt2d_tm(frequency, y_lengths, z_lengths, conductivity)

Inputs:
- Frequency, tensor-mesh cell sizes, and conductivity vector.

Output:
- `Vector{ComplexF64}`: TM boundary field values.

Description:
- Builds the TM boundary conditions from stitched 1D edge solutions.
"""
function get_boundary_mt2d_tm(
    frequency::Float64,
    y_lengths::Vector{Float64},
    z_lengths::Vector{Float64},
    conductivity::Vector{Float64},
)
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    z_nodes = [0.0; cumsum(z_lengths)]
    σ2d = reshape(conductivity, n_y, n_z)'
    boundary = zeros(ComplexF64, 2 * (n_y + n_z))

    boundary[1:n_y+1] .= 1.0 + 0.0im
    σ1d = σ2d[:, 1]
    _, H = mt1d_boundary_field(frequency, σ1d, z_nodes; return_magnetic = true)
    H ./= H[1]
    boundary[n_y+2:n_y+n_z+1] = H[2:end]

    σ1d = σ2d[:, end]
    _, H = mt1d_boundary_field(frequency, σ1d, z_nodes; return_magnetic = true)
    H ./= H[1]
    boundary[n_y+n_z+2:n_y+2*n_z+1] = H[2:end]

    for iy in 2:n_y
        σ1d = (σ2d[:, iy - 1] * y_lengths[iy - 1] + σ2d[:, iy] * y_lengths[iy]) / (y_lengths[iy - 1] + y_lengths[iy])
        _, H = mt1d_boundary_field(frequency, σ1d, z_nodes; return_magnetic = true)
        boundary[n_y + 2 * n_z + iy] = H[end] / H[1]
    end

    boundary
end

"""
    compute_fields_at_receivers_te(ω, receiver_locations, y_nodes, first_cell_thickness, sigma_row, electric_pair)

Inputs:
- Angular frequency, receiver coordinates, surface nodes, first-cell thickness, conductivity row, and electric fields.

Output:
- Tuple of electric and magnetic receiver fields.

Description:
- Interpolates the TE electric and magnetic fields from the solved tensor mesh to the receivers.
"""
function compute_fields_at_receivers_te(
    ω::Float64,
    receiver_locations::Matrix{Float64},
    y_nodes::Vector{Float64},
    first_cell_thickness::Float64,
    sigma_row::Vector{Float64},
    electric_pair::Matrix{ComplexF64},
)
    y_lengths = diff(y_nodes)
    n_y = length(y_lengths)
    μ = μ₀_2D .* ones(n_y)
    n_receivers = size(receiver_locations, 1)

    E = electric_pair[:, 1]
    Hz0 = (ddx(n_y) * electric_pair[:, 1]) ./ y_lengths ./ (1im * ω)
    Hz1 = (ddx(n_y) * electric_pair[:, 2]) ./ y_lengths ./ (1im * ω)
    Hzquarter = (0.75 .* Hz0 .+ 0.25 .* Hz1) ./ μ
    Hyhalf = -(electric_pair[2:end-1, 2] .- electric_pair[2:end-1, 1]) ./ first_cell_thickness ./ (1im * ω * μ₀_2D)
    Equarter = 0.75 .* electric_pair[2:end-1, 1] .+ 0.25 .* electric_pair[2:end-1, 2]
    σavg = (av(n_y - 1) * (sigma_row .* y_lengths)) ./ (av(n_y - 1) * y_lengths)
    ∂Hz∂y = (ddx(n_y - 1) * Hzquarter) ./ (av(n_y - 1) * y_lengths)

    Hysurface = zeros(ComplexF64, n_y + 1)
    Hysurface[2:end-1] = Hyhalf .- (∂Hz∂y .- σavg .* Equarter) .* (0.5 * first_cell_thickness)
    Hysurface[1] = Hysurface[2]
    Hysurface[end] = Hysurface[end - 1]

    Erec = zeros(ComplexF64, n_receivers)
    Hrec = zeros(ComplexF64, n_receivers)
    for i in 1:n_receivers
        y = receiver_locations[i, 1]
        node = findfirst(v -> v > y, y_nodes)
        Δy1 = y - y_nodes[node - 1]
        Δy2 = y_nodes[node] - y
        Erec[i] = E[node - 1] * Δy2 + E[node] * Δy1
        Hrec[i] = Hysurface[node - 1] * Δy2 + Hysurface[node] * Δy1
    end

    Erec, Hrec
end

"""
    compute_fields_at_receivers_tm(ω, receiver_locations, y_nodes, first_cell_thickness, sigma_row, magnetic_pair)

Inputs:
- Angular frequency, receiver coordinates, surface nodes, first-cell thickness, conductivity row, and magnetic fields.

Output:
- Tuple of electric and magnetic receiver fields.

Description:
- Interpolates the TM electric and magnetic fields from the solved tensor mesh to the receivers.
"""
function compute_fields_at_receivers_tm(
    ω::Float64,
    receiver_locations::Matrix{Float64},
    y_nodes::Vector{Float64},
    first_cell_thickness::Float64,
    sigma_row::Vector{Float64},
    magnetic_pair::Matrix{ComplexF64},
)
    y_lengths = diff(y_nodes)
    n_y = length(y_lengths)
    n_receivers = size(receiver_locations, 1)

    H = magnetic_pair[:, 1]
    Jz0 = -(ddx(n_y) * magnetic_pair[:, 1]) ./ y_lengths
    Jz1 = -(ddx(n_y) * magnetic_pair[:, 2]) ./ y_lengths
    Ezquarter = (0.75 .* Jz0 .+ 0.25 .* Jz1) ./ sigma_row
    Jyhalf = (magnetic_pair[2:end-1, 2] .- magnetic_pair[2:end-1, 1]) ./ first_cell_thickness
    ρavg = (av(n_y - 1) * ((1.0 ./ sigma_row) .* y_lengths)) ./ (av(n_y - 1) * y_lengths)
    Eyhalf = Jyhalf .* ρavg
    Hquarter = 0.75 .* magnetic_pair[2:end-1, 1] .+ 0.25 .* magnetic_pair[2:end-1, 2]
    ∂Ez∂y = (ddx(n_y - 1) * Ezquarter) ./ (av(n_y - 1) * y_lengths)

    Eysurface = zeros(ComplexF64, n_y + 1)
    Eysurface[2:end-1] = Eyhalf .- (∂Ez∂y .+ 1im * ω * μ₀_2D .* Hquarter) .* (0.5 * first_cell_thickness)
    Eysurface[1] = Eysurface[2]
    Eysurface[end] = Eysurface[end - 1]

    Erec = zeros(ComplexF64, n_receivers)
    Hrec = zeros(ComplexF64, n_receivers)
    for i in 1:n_receivers
        y = receiver_locations[i, 1]
        node = findfirst(v -> v > y, y_nodes)
        Δy1 = y - y_nodes[node - 1]
        Δy2 = y_nodes[node] - y
        Erec[i] = Eysurface[node - 1] * Δy2 + Eysurface[node] * Δy1
        Hrec[i] = H[node - 1] * Δy2 + H[node] * Δy1
    end

    Erec, Hrec
end

"""
    _phase_fold_to_0_90(phases)

Inputs:
- `phases`: Phase values in degrees.

Output:
- Folded phase values in the range `[0, 90]`.

Description:
- Folds TM phases into the convention used by the package plots and comparisons.
"""
_phase_fold_to_0_90(phases::AbstractArray) = map(ϕ -> begin
    folded = ϕ < 0 ? ϕ + 180 : ϕ
    folded > 90 ? 180 - folded : folded
end, phases)

"""
    compute_mt_response_te(ω, electric, magnetic, data_type)

Inputs:
- Angular frequency, receiver electric fields, receiver magnetic fields, and output type.

Output:
- Matrix of TE apparent resistivity/phase or TE impedance values.

Description:
- Converts TE electric and magnetic fields into MT observables.
"""
function compute_mt_response_te(
    ω::Float64,
    electric::Vector{ComplexF64},
    magnetic::Vector{ComplexF64},
    data_type::String,
)
    Z = electric ./ magnetic
    if occursin("Impedance", data_type)
        return [real(Z) imag(Z)]
    end
    ρ = abs.(Z) .^ 2 ./ (ω * μ₀_2D)
    ϕ = rad2deg.(atan.(imag.(Z), real.(Z)))
    [ρ ϕ]
end

"""
    compute_mt_response_tm(ω, electric, magnetic, data_type; fold_phase=true)

Inputs:
- Angular frequency, receiver electric fields, receiver magnetic fields, output type, and phase-fold flag.

Output:
- Matrix of TM apparent resistivity/phase or TM impedance values.

Description:
- Converts TM electric and magnetic fields into MT observables.
"""
function compute_mt_response_tm(
    ω::Float64,
    electric::Vector{ComplexF64},
    magnetic::Vector{ComplexF64},
    data_type::String;
    fold_phase::Bool = true,
)
    Z = electric ./ magnetic
    if occursin("Impedance", data_type)
        return [real(Z) imag(Z)]
    end
    ρ = abs.(Z) .^ 2 ./ (ω * μ₀_2D)
    ϕ = rad2deg.(atan.(imag.(Z), real.(Z)))
    ϕ = fold_phase ? _phase_fold_to_0_90(ϕ) : ϕ
    [ρ ϕ]
end

"""
    solve_mt2d_te(frequency, mesh, coefficients, receiver_locations, data_type)

Inputs:
- Frequency, tensor mesh, TE coefficient blocks, receiver locations, and requested output type.

Output:
- Tuple with TE observables and the flattened field solution.

Description:
- Solves the TE system for one frequency and samples the receiver responses.
"""
function solve_mt2d_te(
    frequency::Float64,
    mesh::TensorMesh2D,
    coefficients::CoeffMat,
    receiver_locations::Matrix{Float64},
    data_type::String,
)
    y_lengths = mesh.y_lengths
    z_lengths = mesh.z_lengths
    conductivity = mesh.conductivity
    y_nodes = [0.0; cumsum(y_lengths)] .- mesh.origin[1]
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    ω = 2π * frequency

    Aii = coefficients.real_ii + 1im * ω * coefficients.imag_ii
    Aio = coefficients.real_io + 1im * ω * coefficients.imag_io
    boundary = get_boundary_mt2d_te(frequency, y_lengths, z_lengths, conductivity)
    rhs = -Aio * boundary
    interior = lu(Aii) \ rhs

    field = zeros(ComplexF64, n_z + 1, n_y + 1)
    field[1, :] = boundary[1:n_y+1]
    field[2:end, 1] = boundary[n_y+2:n_y+n_z+1]
    field[2:end, end] = boundary[n_y+n_z+2:n_y+2*n_z+1]
    field[end, 2:end-1] = boundary[n_y+2*n_z+2:end]
    field[2:end-1, 2:end-1] = copy(transpose(reshape(interior, n_y - 1, n_z - 1)))

    surface_depth = receiver_locations[1, 2]
    z_index = findfirst(v -> v < 1e-9, abs.(([0.0; cumsum(z_lengths)] .- mesh.origin[2]) .- surface_depth))
    Epair = copy(transpose(field[z_index:z_index+1, :]))
    σrow = conductivity[(z_index - 1) * n_y + 1:z_index * n_y]
    Δz1 = z_lengths[z_index]
    Erec, Hrec = compute_fields_at_receivers_te(ω, receiver_locations, y_nodes, Δz1, σrow, Epair)

    response = compute_mt_response_te(ω, Erec, Hrec, data_type)
    response, vec(copy(transpose(field)))
end

"""
    solve_mt2d_tm(frequency, mesh, coefficients, receiver_locations, data_type; fold_phase=true)

Inputs:
- Frequency, tensor mesh, TM coefficient blocks, receiver locations, output type, and phase-fold flag.

Output:
- Tuple with TM observables and the flattened field solution.

Description:
- Solves the TM system for one frequency and samples the receiver responses.
"""
function solve_mt2d_tm(
    frequency::Float64,
    mesh::TensorMesh2D,
    coefficients::CoeffMat,
    receiver_locations::Matrix{Float64},
    data_type::String;
    fold_phase::Bool = true,
)
    y_lengths = mesh.y_lengths
    z_lengths = mesh.z_lengths
    conductivity = mesh.conductivity
    y_nodes = [0.0; cumsum(y_lengths)] .- mesh.origin[1]
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    ω = 2π * frequency

    Aii = coefficients.real_ii + 1im * ω * coefficients.imag_ii
    Aio = coefficients.real_io + 1im * ω * coefficients.imag_io
    boundary = get_boundary_mt2d_tm(frequency, y_lengths, z_lengths, conductivity)
    rhs = -Aio * boundary
    interior = lu(Aii) \ rhs

    field = zeros(ComplexF64, n_z + 1, n_y + 1)
    field[1, :] = boundary[1:n_y+1]
    field[2:end, 1] = boundary[n_y+2:n_y+n_z+1]
    field[2:end, end] = boundary[n_y+n_z+2:n_y+2*n_z+1]
    field[end, 2:end-1] = boundary[n_y+2*n_z+2:end]
    field[2:end-1, 2:end-1] = copy(transpose(reshape(interior, n_y - 1, n_z - 1)))

    surface_depth = receiver_locations[1, 2]
    z_index = findfirst(v -> v < 1e-9, abs.(([0.0; cumsum(z_lengths)] .- mesh.origin[2]) .- surface_depth))
    Hpair = copy(transpose(field[z_index:z_index+1, :]))
    σrow = conductivity[(z_index - 1) * n_y + 1:z_index * n_y]
    Δz1 = z_lengths[z_index]
    Erec, Hrec = compute_fields_at_receivers_tm(ω, receiver_locations, y_nodes, Δz1, σrow, Hpair)

    response = compute_mt_response_tm(ω, Erec, Hrec, data_type; fold_phase = fold_phase)
    response, vec(copy(transpose(field)))
end

"""
    _assemble_mt2d_system(mesh, resistivity)

Inputs:
- `mesh`: 2D MT mesh.
- `resistivity`: Cell resistivity model.

Output:
- Tuple with the tensor mesh, TE coefficients, TM coefficients, and receiver locations.

Description:
- Assembles the sparse finite-volume operators needed for the 2D forward solve.
"""
function _assemble_mt2d_system(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real})
    n_z = length(mesh.z_cell_sizes)
    n_y = length(mesh.y_cell_sizes)
    size(resistivity) == (n_z, n_y) || error("resistivity must be size ($(n_z), $(n_y))")

    σ = 1.0 ./ Matrix{Float64}(resistivity)
    σ[1:mesh.n_air_cells, :] .= 1e-9

    tensor_mesh = TensorMesh2D(
        mesh.y_cell_sizes,
        mesh.z_cell_sizes;
        origin = [0.0, 0.0],
        conductivity = vec(copy(transpose(σ))),
    )
    setup_tensor_mesh_2d!(tensor_mesh)

    receiver_y = mesh.receiver_positions .- first(mesh.y_nodes)
    receiver_z = abs(first(mesh.z_nodes))
    receiver_locations = hcat(receiver_y, fill(receiver_z, length(mesh.receiver_positions)))

    n_y_cells, n_z_cells = tensor_mesh.grid_size
    interior, exterior = get_boundary_index(n_y_cells, n_z_cells)

    σcell_to_node = tensor_mesh.average_cell_to_node * (tensor_mesh.face * tensor_mesh.conductivity) |> values -> sparse(Diagonal(values))
    μface = tensor_mesh.average_cell_to_face * (tensor_mesh.face * (1.0 ./ (μ₀_2D .* ones(n_y_cells * n_z_cells)))) |> values -> sparse(Diagonal(values))
    grad_te = tensor_mesh.gradient' * μface * tensor_mesh.gradient
    coeffs_te = CoeffMat(
        grad_te[interior, interior],
        σcell_to_node[interior, interior],
        grad_te[interior, exterior],
        σcell_to_node[interior, exterior],
    )

    μcell_to_node = tensor_mesh.average_cell_to_node * (tensor_mesh.face * (μ₀_2D .* ones(n_y_cells * n_z_cells))) |> values -> sparse(Diagonal(values))
    σface = tensor_mesh.average_cell_to_face * (tensor_mesh.face * (1.0 ./ tensor_mesh.conductivity)) |> values -> sparse(Diagonal(values))
    grad_tm = tensor_mesh.gradient' * σface * tensor_mesh.gradient
    coeffs_tm = CoeffMat(
        grad_tm[interior, interior],
        μcell_to_node[interior, interior],
        grad_tm[interior, exterior],
        μcell_to_node[interior, exterior],
    )

    tensor_mesh, coeffs_te, coeffs_tm, receiver_locations
end

"""
    run_mt2d_forward(mesh, resistivity; mode=:TETM)

Inputs:
- `mesh`: 2D MT mesh.
- `resistivity`: Cell resistivity model.
- `mode`: `:TE`, `:TM`, or `:TETM`.

Output:
- `MT2DResponse`: TE/TM MT response.

Description:
- Runs the 2D forward solver across all survey frequencies.
"""
function run_mt2d_forward(
    mesh::MT2DMesh,
    resistivity::AbstractMatrix{<:Real};
    mode::Symbol = :TETM,
)
    tensor_mesh, coeffs_te, coeffs_tm, receiver_locations = _assemble_mt2d_system(mesh, resistivity)

    n_f = length(mesh.frequencies)
    n_r = length(mesh.receiver_positions)
    rho_xy = zeros(n_f, n_r)
    phase_xy = zeros(n_f, n_r)
    z_xy = zeros(ComplexF64, n_f, n_r)
    rho_yx = zeros(n_f, n_r)
    phase_yx = zeros(n_f, n_r)
    z_yx = zeros(ComplexF64, n_f, n_r)

    for (i, f) in enumerate(mesh.frequencies)
        if mode in (:TE, :TETM)
            response_te, _ = solve_mt2d_te(f, tensor_mesh, coeffs_te, receiver_locations, "Rho_Pha")
            rho_xy[i, :] .= response_te[:, 1]
            phase_xy[i, :] .= response_te[:, 2]
            impedance_te, _ = solve_mt2d_te(f, tensor_mesh, coeffs_te, receiver_locations, "Impedance")
            z_xy[i, :] .= complex.(impedance_te[:, 1], impedance_te[:, 2])
        end
        if mode in (:TM, :TETM)
            response_tm, _ = solve_mt2d_tm(f, tensor_mesh, coeffs_tm, receiver_locations, "Rho_Pha"; fold_phase = false)
            rho_yx[i, :] .= response_tm[:, 1]
            phase_yx[i, :] .= response_tm[:, 2]
            impedance_tm, _ = solve_mt2d_tm(f, tensor_mesh, coeffs_tm, receiver_locations, "Impedance"; fold_phase = false)
            z_yx[i, :] .= complex.(impedance_tm[:, 1], impedance_tm[:, 2])
        end
    end

    MT2DResponse(
        frequencies = mesh.frequencies,
        periods = 1.0 ./ mesh.frequencies,
        receivers = mesh.receiver_positions,
        rho_xy = rho_xy,
        phase_xy = phase_xy,
        z_xy = z_xy,
        rho_yx = rho_yx,
        phase_yx = phase_yx,
        z_yx = z_yx,
    )
end

"""
    write_mt2d_response_csv(path, response)

Inputs:
- `path`: Output CSV path.
- `response`: 2D MT response.

Output:
- `String`: Path to the written CSV file.

Description:
- Writes TE and TM response curves for all stations and frequencies to a CSV file.
"""
function write_mt2d_response_csv(path::AbstractString, response::MT2DResponse)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "frequency_hz,period_s,receiver_m,rho_xy_ohm_m,phase_xy_deg,rho_yx_ohm_m,phase_yx_deg,re_zxy,im_zxy,re_zyx,im_zyx")
        for ir in eachindex(response.receivers)
            for ifreq in eachindex(response.frequencies)
                @printf(
                    io,
                    "%.8e,%.8e,%.3f,%.8e,%.6f,%.8e,%.6f,%.8e,%.8e,%.8e,%.8e\n",
                    response.frequencies[ifreq],
                    response.periods[ifreq],
                    response.receivers[ir],
                    response.rho_xy[ifreq, ir],
                    response.phase_xy[ifreq, ir],
                    response.rho_yx[ifreq, ir],
                    response.phase_yx[ifreq, ir],
                    real(response.z_xy[ifreq, ir]),
                    imag(response.z_xy[ifreq, ir]),
                    real(response.z_yx[ifreq, ir]),
                    imag(response.z_yx[ifreq, ir]),
                )
            end
        end
    end
    String(path)
end

"""
    Forward2D(mesh, resistivity)

Inputs:
- `mesh`: 2D MT mesh.
- `resistivity`: Cell resistivity model.

Output:
- `MT2DResponse`: TE/TM response.

Description:
- Public alias for `run_mt2d_forward`.
"""
Forward2D(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real}) = run_mt2d_forward(mesh, resistivity)

"""
    RunForward2D(mesh, resistivity)

Inputs:
- `mesh`: 2D MT mesh.
- `resistivity`: Cell resistivity model.

Output:
- `MT2DResponse`: TE/TM response.

Description:
- Compatibility alias for `run_mt2d_forward`.
"""
RunForward2D(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real}) = run_mt2d_forward(mesh, resistivity)

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

"""
    _impedance_to_rho_phase(impedance, frequency)

Inputs:
- `impedance`: Complex impedance.
- `frequency`: Frequency in hertz.

Output:
- Named tuple with `rho` and `phase`.

Description:
- Converts a complex impedance into apparent resistivity and phase.
"""
function _impedance_to_rho_phase(impedance::ComplexF64, frequency::Float64)
    if !isfinite(real(impedance)) || !isfinite(imag(impedance))
        return (rho = NaN, phase = NaN)
    end
    ω = 2π * frequency
    ρ = abs2(impedance) / (μ₀_2D * ω)
    ϕ = rad2deg(angle(impedance))
    (rho = ρ, phase = ϕ)
end

"""
    data_from_response2d(response; z_xy_error=nothing, z_yx_error=nothing, z_xx_error=nothing, z_yy_error=nothing, impedance_error_fraction=0.05, title="MTGeophysics.jl 2D profile data", site_names=nothing, x_positions=nothing, z_positions=nothing)

Inputs:
- `response`: 2D MT response.
- Optional error floors, metadata, and station coordinates.

Output:
- `DataFile2D`: 2D data container.

Description:
- Converts a forward response into the package data-file representation.
"""
function data_from_response2d(
    response::MT2DResponse;
    z_xy_error::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    z_yx_error::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    z_xx_error::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    z_yy_error::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
    impedance_error_fraction::Real = 0.05,
    title::AbstractString = "MTGeophysics.jl 2D profile data",
    site_names::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    x_positions::Union{Nothing, AbstractVector{<:Real}} = nothing,
    z_positions::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    n_f, n_r = size(response.z_xy)
    sites = site_names === nothing ? [@sprintf("Site%03d", i) for i in 1:n_r] : String.(site_names)
    xvals = x_positions === nothing ? zeros(n_r) : Float64.(x_positions)
    zvals = z_positions === nothing ? zeros(n_r) : Float64.(z_positions)

    default_xy_error = impedance_error_fraction .* abs.(response.z_xy)
    default_yx_error = impedance_error_fraction .* abs.(response.z_yx)
    default_xx_error = fill(max(impedance_error_fraction, 1e-6), n_f, n_r)
    default_yy_error = fill(max(impedance_error_fraction, 1e-6), n_f, n_r)

    data = DataFile2D(
        title = String(title),
        periods = Float64.(response.periods),
        frequencies = Float64.(response.frequencies),
        site_names = sites,
        receivers = Float64.(response.receivers),
        x_positions = xvals,
        z_positions = zvals,
        z_xy = ComplexF64.(response.z_xy),
        z_xy_error = z_xy_error === nothing ? default_xy_error : Float64.(z_xy_error),
        z_yx = ComplexF64.(response.z_yx),
        z_yx_error = z_yx_error === nothing ? default_yx_error : Float64.(z_yx_error),
        z_xx = zeros(ComplexF64, n_f, n_r),
        z_xx_error = z_xx_error === nothing ? default_xx_error : Float64.(z_xx_error),
        z_yy = zeros(ComplexF64, n_f, n_r),
        z_yy_error = z_yy_error === nothing ? default_yy_error : Float64.(z_yy_error),
        rho_xy = Float64.(response.rho_xy),
        phase_xy = Float64.(response.phase_xy),
        rho_yx = Float64.(response.rho_yx),
        phase_yx = Float64.(response.phase_yx),
    )
    size(data.z_xy_error) == size(data.z_xy) || error("z_xy_error does not match z_xy size")
    size(data.z_yx_error) == size(data.z_yx) || error("z_yx_error does not match z_yx size")
    data
end

function _resolve_forwardsolve2d_errors(template::DataFile2D, response::MT2DResponse, data_path::AbstractString)
    is_reference_template = lowercase(splitext(data_path)[2]) == ".ref" ||
        (all(iszero, template.z_xy) && all(iszero, template.z_yx))
    if is_reference_template
        return (
            z_xy_error = max.(Float64.(template.z_xy_error) .* abs.(response.z_xy), 1e-6),
            z_yx_error = max.(Float64.(template.z_yx_error) .* abs.(response.z_yx), 1e-6),
            z_xx_error = Float64.(template.z_xx_error),
            z_yy_error = Float64.(template.z_yy_error),
        )
    end

    (
        z_xy_error = Float64.(template.z_xy_error),
        z_yx_error = Float64.(template.z_yx_error),
        z_xx_error = Float64.(template.z_xx_error),
        z_yy_error = Float64.(template.z_yy_error),
    )
end

"""
    write_data2d(path, response_or_data; impedance_error_fraction=0.05, title="MTGeophysics.jl 2D profile data")

Inputs:
- Output path and either a `MT2DResponse` or `DataFile2D`.

Output:
- `String`: Path to the written data file.

Description:
- Writes a 2D MT data file, creating one from a response if needed.
"""
function write_data2d(
    path::AbstractString,
    response_or_data;
    impedance_error_fraction::Real = 0.05,
    title::AbstractString = "MTGeophysics.jl 2D profile data",
)
    data = response_or_data isa DataFile2D ?
        response_or_data :
        data_from_response2d(response_or_data; impedance_error_fraction = impedance_error_fraction, title = title)

    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# $(data.title)")
        println(io, "# Generated by MTGeophysics.jl/src/MTGeophysics2D.jl")
        println(io, "# Period(s) Site Lat Lon X(m) Y(m) Z(m) Component Real Imag Error")
        println(io, "> Full_Impedance")
        println(io, "> exp(+i\\omega t)")
        println(io, "> [V/m]/[T]")
        println(io, "> 0.00")
        println(io, "> 0.000 0.000")
        println(io, @sprintf("> %d %d", length(data.receivers), length(data.periods) * length(data.receivers) * 4))

        for ir in eachindex(data.receivers)
            site = data.site_names[ir]
            x = data.x_positions[ir]
            y = data.receivers[ir]
            z = data.z_positions[ir]
            for ifreq in eachindex(data.periods)
                period = data.periods[ifreq]
                entries = (
                    ("ZXX", data.z_xx[ifreq, ir], data.z_xx_error[ifreq, ir]),
                    ("ZXY", data.z_xy[ifreq, ir], data.z_xy_error[ifreq, ir]),
                    ("ZYX", data.z_yx[ifreq, ir], data.z_yx_error[ifreq, ir]),
                    ("ZYY", data.z_yy[ifreq, ir], data.z_yy_error[ifreq, ir]),
                )
                for (component, impedance, σ) in entries
                    @printf(
                        io,
                        "%.8e %s %.6f %.6f %.3f %.3f %.3f %s %.8e %.8e %.8e\n",
                        period,
                        site,
                        0.0,
                        0.0,
                        x,
                        y,
                        z,
                        component,
                        real(impedance),
                        imag(impedance),
                        Float64(σ),
                    )
                end
            end
        end
    end

    String(path)
end

"""
    load_data2d(path)

Inputs:
- `path`: Path to a 2D MT data file.

Output:
- `DataFile2D`: Parsed data file.

Description:
- Reads a 2D MT impedance-data file written by the package.
"""
function load_data2d(path::AbstractString)
    isfile(path) || error("data file not found: $path")
    lines = readlines(path)

    title = "MTGeophysics.jl 2D profile data"
    rows = NamedTuple[]
    for line in lines
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, "#")
            title == "MTGeophysics.jl 2D profile data" && (title = strip(replace(stripped, "#" => "")))
            continue
        end
        startswith(stripped, ">") && continue

        tokens = split(stripped)
        length(tokens) < 11 && continue
        push!(
            rows,
            (
                period = parse(Float64, tokens[1]),
                site = String(tokens[2]),
                x = parse(Float64, tokens[5]),
                y = parse(Float64, tokens[6]),
                z = parse(Float64, tokens[7]),
                component = uppercase(tokens[8]),
                real = parse(Float64, tokens[9]),
                imag = parse(Float64, tokens[10]),
                sigma = parse(Float64, tokens[11]),
            ),
        )
    end
    isempty(rows) && error("no impedance rows were found in $path")

    periods = unique([row.period for row in rows])
    frequencies = 1.0 ./ periods

    site_names = String[]
    receiver_lookup = Dict{String, Tuple{Float64, Float64, Float64}}()
    for row in rows
        if !haskey(receiver_lookup, row.site)
            push!(site_names, row.site)
            receiver_lookup[row.site] = (row.x, row.y, row.z)
        end
    end
    n_periods = length(periods)
    n_receivers = length(site_names)

    z_xy = fill(ComplexF64(NaN, NaN), n_periods, n_receivers)
    z_xy_error = fill(NaN, n_periods, n_receivers)
    z_yx = fill(ComplexF64(NaN, NaN), n_periods, n_receivers)
    z_yx_error = fill(NaN, n_periods, n_receivers)
    z_xx = fill(ComplexF64(0.0, 0.0), n_periods, n_receivers)
    z_xx_error = fill(NaN, n_periods, n_receivers)
    z_yy = fill(ComplexF64(0.0, 0.0), n_periods, n_receivers)
    z_yy_error = fill(NaN, n_periods, n_receivers)

    period_index = Dict(period => index for (index, period) in enumerate(periods))
    site_index = Dict(site => index for (index, site) in enumerate(site_names))
    for row in rows
        ip = period_index[row.period]
        is = site_index[row.site]
        impedance = ComplexF64(row.real, row.imag)
        if row.component == "ZXY"
            z_xy[ip, is] = impedance
            z_xy_error[ip, is] = row.sigma
        elseif row.component == "ZYX"
            z_yx[ip, is] = impedance
            z_yx_error[ip, is] = row.sigma
        elseif row.component == "ZXX"
            z_xx[ip, is] = impedance
            z_xx_error[ip, is] = row.sigma
        elseif row.component == "ZYY"
            z_yy[ip, is] = impedance
            z_yy_error[ip, is] = row.sigma
        end
    end

    rho_xy = Matrix{Float64}(undef, n_periods, n_receivers)
    phase_xy = similar(rho_xy)
    rho_yx = similar(rho_xy)
    phase_yx = similar(rho_xy)
    for is in 1:n_receivers, ip in 1:n_periods
        values_xy = _impedance_to_rho_phase(z_xy[ip, is], frequencies[ip])
        rho_xy[ip, is] = values_xy.rho
        phase_xy[ip, is] = values_xy.phase
        values_yx = _impedance_to_rho_phase(z_yx[ip, is], frequencies[ip])
        rho_yx[ip, is] = values_yx.rho
        phase_yx[ip, is] = values_yx.phase
    end

    coordinates = [receiver_lookup[site] for site in site_names]
    x_positions = [value[1] for value in coordinates]
    receivers = [value[2] for value in coordinates]
    z_positions = [value[3] for value in coordinates]

    DataFile2D(
        title = title,
        periods = Float64.(periods),
        frequencies = Float64.(frequencies),
        site_names = site_names,
        receivers = receivers,
        x_positions = x_positions,
        z_positions = z_positions,
        z_xy = z_xy,
        z_xy_error = z_xy_error,
        z_yx = z_yx,
        z_yx_error = z_yx_error,
        z_xx = z_xx,
        z_xx_error = z_xx_error,
        z_yy = z_yy,
        z_yy_error = z_yy_error,
        rho_xy = rho_xy,
        phase_xy = phase_xy,
        rho_yx = rho_yx,
        phase_yx = phase_yx,
        path = String(path),
    )
end

"""
    data_to_response2d(data)

Inputs:
- `data`: Parsed 2D data file.

Output:
- `MT2DResponse`: Response container derived from the data file.

Description:
- Converts a `DataFile2D` object into the forward-response container used by plotting code.
"""
function data_to_response2d(data::DataFile2D)
    MT2DResponse(
        frequencies = Float64.(data.frequencies),
        periods = Float64.(data.periods),
        receivers = Float64.(data.receivers),
        rho_xy = Float64.(data.rho_xy),
        phase_xy = Float64.(data.phase_xy),
        z_xy = ComplexF64.(data.z_xy),
        rho_yx = Float64.(data.rho_yx),
        phase_yx = Float64.(data.phase_yx),
        z_yx = ComplexF64.(data.z_yx),
    )
end

"""
    chi2_rms2d(observed, predicted; components=["ZXY", "ZYX"])

Inputs:
- Observed and predicted 2D data objects.
- Optional list of tensor components to evaluate.

Output:
- `FitSummary2D`: Chi-square, RMS, and sample count.

Description:
- Computes the weighted complex-impedance misfit between observed and predicted 2D data.
"""
function chi2_rms2d(
    observed::DataFile2D,
    predicted::DataFile2D;
    components::AbstractVector{<:AbstractString} = ["ZXY", "ZYX"],
)
    size(observed.z_xy) == size(predicted.z_xy) || error("observed and predicted surveys do not match")
    χ² = 0.0
    count = 0

    for component in uppercase.(String.(components))
        if component == "ZXY"
            for index in eachindex(observed.z_xy)
                zo = observed.z_xy[index]
                zp = predicted.z_xy[index]
                σ = observed.z_xy_error[index]
                if isfinite(real(zo)) && isfinite(imag(zo)) && isfinite(real(zp)) && isfinite(imag(zp)) && isfinite(σ) && σ > 0
                    residual_real = (real(zp) - real(zo)) / σ
                    residual_imag = (imag(zp) - imag(zo)) / σ
                    χ² += residual_real^2 + residual_imag^2
                    count += 2
                end
            end
        elseif component == "ZYX"
            for index in eachindex(observed.z_yx)
                zo = observed.z_yx[index]
                zp = predicted.z_yx[index]
                σ = observed.z_yx_error[index]
                if isfinite(real(zo)) && isfinite(imag(zo)) && isfinite(real(zp)) && isfinite(imag(zp)) && isfinite(σ) && σ > 0
                    residual_real = (real(zp) - real(zo)) / σ
                    residual_imag = (imag(zp) - imag(zo)) / σ
                    χ² += residual_real^2 + residual_imag^2
                    count += 2
                end
            end
        end
    end

    FitSummary2D(
        chi2 = χ²,
        rms = count > 0 ? sqrt(χ² / count) : NaN,
        count = count,
    )
end

"""
    chi2_rms2d(observed_path, predicted_path; components=["ZXY", "ZYX"])

Inputs:
- Paths to observed and predicted 2D data files.

Output:
- `FitSummary2D`: Chi-square, RMS, and sample count.

Description:
- Convenience method that loads data files before calling `chi2_rms2d`.
"""
chi2_rms2d(observed_path::AbstractString, predicted_path::AbstractString; components::AbstractVector{<:AbstractString} = ["ZXY", "ZYX"]) =
    chi2_rms2d(load_data2d(observed_path), load_data2d(predicted_path); components = components)

"""
    build_mt2d_data_template(mesh; impedance_error_fraction=0.05, title="2D data template")

Inputs:
- `mesh`: 2D MT mesh.
- Impedance error floor and optional title.

Output:
- `DataFile2D`: Empty data template with relative impedance-error fractions.

Description:
- Builds the package 2D data-template object from a survey mesh.
- The stored impedance-error values are interpreted as relative fractions when the
- template is later consumed through a `.ref` workflow.
"""
function build_mt2d_data_template(
    mesh::MT2DMesh;
    impedance_error_fraction::Real = 0.05,
    title::AbstractString = "2D data template",
)
    n_f = length(mesh.frequencies)
    n_r = length(mesh.receiver_positions)
    error_fraction = fill(max(Float64(impedance_error_fraction), 1e-6), n_f, n_r)

    DataFile2D(
        title = String(title),
        periods = 1.0 ./ mesh.frequencies,
        frequencies = Float64.(mesh.frequencies),
        site_names = [@sprintf("Site%03d", index) for index in 1:n_r],
        receivers = Float64.(mesh.receiver_positions),
        x_positions = zeros(n_r),
        z_positions = zeros(n_r),
        z_xy = zeros(ComplexF64, n_f, n_r),
        z_xy_error = copy(error_fraction),
        z_yx = zeros(ComplexF64, n_f, n_r),
        z_yx_error = copy(error_fraction),
        z_xx = zeros(ComplexF64, n_f, n_r),
        z_xx_error = copy(error_fraction),
        z_yy = zeros(ComplexF64, n_f, n_r),
        z_yy_error = copy(error_fraction),
        rho_xy = fill(NaN, n_f, n_r),
        phase_xy = fill(NaN, n_f, n_r),
        rho_yx = fill(NaN, n_f, n_r),
        phase_yx = fill(NaN, n_f, n_r),
    )
end

"""
    write_mt2d_data_template(path, mesh; impedance_error_fraction=0.05)

Inputs:
- Output path, 2D MT mesh, and impedance error floor.

Output:
- `String`: Path to the written template.

Description:
- Writes a 2D data template to disk.
"""
function write_mt2d_data_template(
    path::AbstractString,
    mesh::MT2DMesh;
    impedance_error_fraction::Real = 0.05,
)
    write_data2d(path, build_mt2d_data_template(mesh; impedance_error_fraction = impedance_error_fraction))
end

"""
    _apply_mt2d_noise(data; rng_seed)

Inputs:
- `data`: Noise-free 2D data object.
- `rng_seed`: Random-number seed.

Output:
- `DataFile2D`: Noisy data object.

Description:
- Perturbs the complex impedance entries using the stored impedance error floors.
"""
function _apply_mt2d_noise(data::DataFile2D; rng_seed::Integer)
    noisy = deepcopy(data)
    rng = MersenneTwister(rng_seed)

    for ir in axes(noisy.z_xy, 2), ifreq in axes(noisy.z_xy, 1)
        xy_error = noisy.z_xy_error[ifreq, ir] / sqrt(2)
        yx_error = noisy.z_yx_error[ifreq, ir] / sqrt(2)
        noisy.z_xy[ifreq, ir] += xy_error * (randn(rng) + im * randn(rng))
        noisy.z_yx[ifreq, ir] += yx_error * (randn(rng) + im * randn(rng))

        xy = _impedance_to_rho_phase(noisy.z_xy[ifreq, ir], noisy.frequencies[ifreq])
        yx = _impedance_to_rho_phase(noisy.z_yx[ifreq, ir], noisy.frequencies[ifreq])
        noisy.rho_xy[ifreq, ir] = xy.rho
        noisy.phase_xy[ifreq, ir] = xy.phase
        noisy.rho_yx[ifreq, ir] = yx.rho
        noisy.phase_yx[ifreq, ir] = yx.phase
    end

    noisy
end

"""
    ForwardSolve2D(model_path, data_path; add_noise=false, output_path=nothing, rng_seed=20260308)

Inputs:
- `model_path`: Saved 2D model file.
- `data_path`: 2D data template or observed-data file.
- `add_noise`, `output_path`, `rng_seed`: File-driven forward options.

Output:
- `String`: Path to the written predicted or observed data file.

Description:
- Loads a saved 2D model, runs the forward solver, optionally adds noise, and writes the result to disk.
"""
function ForwardSolve2D(
    model_path::AbstractString,
    data_path::AbstractString;
    add_noise::Bool = false,
    output_path::Union{Nothing, AbstractString} = nothing,
    rng_seed::Integer = 20260308,
)
    model = load_model2d(model_path)
    template = load_data2d(data_path)
    mesh = build_mesh_from_model2d(
        model;
        frequencies = template.frequencies,
        receiver_positions = template.receivers,
    )

    response = run_mt2d_forward(mesh, model.resistivity)
    resolved_errors = _resolve_forwardsolve2d_errors(template, response, data_path)
    predicted = data_from_response2d(
        response;
        z_xy_error = resolved_errors.z_xy_error,
        z_yx_error = resolved_errors.z_yx_error,
        z_xx_error = resolved_errors.z_xx_error,
        z_yy_error = resolved_errors.z_yy_error,
        site_names = template.site_names,
        x_positions = template.x_positions,
        z_positions = template.z_positions,
        title = add_noise ? "2D observed data" : "2D predicted data",
    )

    observed = add_noise ? _apply_mt2d_noise(predicted; rng_seed = rng_seed) : predicted
    destination = something(output_path, joinpath(dirname(data_path), "Data.obs"))
    write_data2d(destination, observed)
end

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
    plot_mt2d_model(mesh, resistivity; output_path, show_air=false, show_grid=false, show_padding=true, maximum_depth_km=Inf, resistivity_log10_range=(0.0, 4.0))

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
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    core_y = _plot_core_indices_2d(mesh.y_cell_sizes)
    column_range = show_padding ? (1:size(resistivity, 2)) : core_y
    y_edges = mesh.y_nodes[first(column_range):(last(column_range) + 1)] ./ 1000
    row_range = show_air ? (1:size(resistivity, 1)) : ((mesh.n_air_cells + 1):size(resistivity, 1))
    z_edges = show_air ? (mesh.z_nodes ./ 1000) : (mesh.z_nodes[(mesh.n_air_cells + 1):end] ./ 1000)
    rho_plot = log10.(resistivity[row_range, column_range])

    figure = Figure(size = (1100, 650))
    axis = Axis(
        figure[1, 1],
        xlabel = "Offset (km)",
        ylabel = "Depth (km)",
        yreversed = !show_air,
        title = show_air ? "2D resistivity model with air" : "2D resistivity model",
    )
    heatmap = heatmap!(axis, y_edges, z_edges, rho_plot', colormap = :Spectral, colorrange = resistivity_log10_range)
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

"""
    phase180(values)

Inputs:
- Phase values in degrees.

Output:
- Phase values folded to `[-180, 180]`.

Description:
- Wraps phase angles for the 2D map plot.
"""
phase180(values) = ((values .+ 180) .% 360) .- 180

"""
    plot_mt2d_data_maps(response; output_path)

Inputs:
- 2D MT response and output image path.

Output:
- `String`: Path to the written plot.

Description:
- Writes the standard 2D TE/TM response maps.
"""
function plot_mt2d_data_maps(
    response::MT2DResponse;
    output_path::AbstractString,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    rx_km = response.receivers ./ 1000
    periods = response.periods
    log_periods = log10.(periods)
    tm_phase = _phase_fold_to_0_90(response.phase_yx)
    figure = Figure(size = (1400, 900))

    # Clamp negative apparent resistivities (numerical artifacts) to eps() before log10
    rho_xy_safe = max.(response.rho_xy, eps())
    rho_yx_safe = max.(response.rho_yx, eps())

    ax1 = Axis(figure[1, 1], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TE log10(ρxy)")
    hm1 = heatmap!(ax1, rx_km, log_periods, log10.(rho_xy_safe)', colormap = :Spectral)
    Colorbar(figure[1, 2], hm1)

    ax2 = Axis(figure[1, 3], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TM log10(ρyx)")
    hm2 = heatmap!(ax2, rx_km, log_periods, log10.(rho_yx_safe)', colormap = :Spectral)
    Colorbar(figure[1, 4], hm2)

    ax3 = Axis(figure[2, 1], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TE phase")
    hm3 = heatmap!(ax3, rx_km, log_periods, phase180(response.phase_xy)', colormap = :Spectral, colorrange = (-180, 180))
    Colorbar(figure[2, 2], hm3, ticks = [-180, -90, 0, 90, 180])

    ax4 = Axis(figure[2, 3], xlabel = "Position (km)", ylabel = "log10 Period (s)", title = "TM phase")
    hm4 = heatmap!(ax4, rx_km, log_periods, tm_phase', colormap = :Spectral, colorrange = (0, 90))
    Colorbar(figure[2, 4], hm4, ticks = [0, 30, 60, 90])

    save(output_path, figure)
    String(output_path)
end

"""
    plot_mt2d_site_curves(response; station_index=cld(length(response.receivers), 2), output_path)

Inputs:
- 2D MT response, station index, and output image path.

Output:
- `String`: Path to the written plot.

Description:
- Writes the representative-station TE/TM response curves.
"""
function plot_mt2d_site_curves(
    response::MT2DResponse;
    station_index::Integer = cld(length(response.receivers), 2),
    output_path::AbstractString,
)
    CairoMakie.activate!()
    mkpath(dirname(output_path))

    tm_phase = _phase_fold_to_0_90(response.phase_yx[:, station_index])
    figure = Figure(size = (900, 700))
    rho_axis = Axis(
        figure[1, 1],
        xlabel = "Period (s)",
        ylabel = "Apparent resistivity (Ω·m)",
        xscale = log10,
        yscale = log10,
        title = "Station $(station_index) response",
    )
    phase_axis = Axis(
        figure[2, 1],
        xlabel = "Period (s)",
        ylabel = "Phase (deg)",
        xscale = log10,
        title = "Phase",
    )

    # Clamp negative apparent resistivities (numerical artifacts) to eps() before log-scale plot
    rho_xy_site = max.(response.rho_xy[:, station_index], eps())
    rho_yx_site = max.(response.rho_yx[:, station_index], eps())
    lines!(rho_axis, response.periods, rho_xy_site, color = :navy, linewidth = 3, label = "TE")
    lines!(rho_axis, response.periods, rho_yx_site, color = :darkorange, linewidth = 3, label = "TM")
    lines!(phase_axis, response.periods, response.phase_xy[:, station_index], color = :navy, linewidth = 3, label = "TE")
    lines!(phase_axis, response.periods, tm_phase, color = :darkorange, linewidth = 3, label = "TM")
    axislegend(rho_axis, position = :rb)
    axislegend(phase_axis, position = :rb)

    save(output_path, figure)
    String(output_path)
end

"""
    PlotData2D(data_path; maps_output_path, curves_output_path, station_index=nothing)

Inputs:
- 2D data path, map image path, curve image path, and optional station index.

Output:
- Named tuple with the written plot paths.

Description:
- Loads a 2D data file and writes both the response maps and a representative-station plot.
"""
function PlotData2D(
    data_path::AbstractString;
    maps_output_path::AbstractString,
    curves_output_path::AbstractString,
    station_index::Union{Nothing, Int} = nothing,
)
    data = load_data2d(data_path)
    response = data_to_response2d(data)
    plot_mt2d_data_maps(response; output_path = maps_output_path)
    plot_mt2d_site_curves(
        response;
        station_index = something(station_index, cld(length(response.receivers), 2)),
        output_path = curves_output_path,
    )
    (maps_output_path = maps_output_path, curves_output_path = curves_output_path)
end

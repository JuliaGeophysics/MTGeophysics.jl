# 2D MT forward modelling
# Author: @pankajkmishra
# Finite-difference (box integration) TE/TM solver, data file I/O, misfit, and Fréchet derivatives
# Fréchet derivatives in Tarantola (2005) notation, d = g(m), G = ∂g/∂m: FrechetDerivative2D (G),
# ApplyFrechet2D (δd = G δm), ApplyFrechetTranspose2D (δm̂ = Gᵗ δd̂); the discrete system is
# differentiated implicitly, no finite differences

#***********************************************************************
# Description: 2D MT forward problem
#
#   convention   exp(+iωt), quasi-static, μ = μ₀
#                strike along x, y along the profile, z positive down
#                from the top of the air layer
#   maxwell      ∇×E = -iωμ H,   ∇×H = σ E
#   domain       Ω = [y₀, y₁] × [0, z_max], air included with σ_air = 1e-9 S/m
#
#   TE (E-polarisation), E = Ex(y,z) x̂
#     ∇·(μ⁻¹ ∇Ex) - iωσ Ex = 0                                   in Ω
#     Hy = -(1/iωμ) ∂z Ex,   Hz = (1/iωμ) ∂y Ex
#
#   TM (H-polarisation), H = Hx(y,z) x̂,  ρ = 1/σ
#     ∇·(ρ ∇Hx) - iωμ Hx = 0                                     in Ω
#     Ey = ρ ∂z Hx,   Ez = -ρ ∂y Hx
#
#   boundary conditions, dirichlet on ∂Ω, u = Ex (TE) or Hx (TM)
#     top      u(y, 0) = 1
#     sides    u(y₀, z) = u¹ᴰ[σ(y₀, ·)](z),   u(y₁, z) = u¹ᴰ[σ(y₁, ·)](z)
#     bottom   u(y, z_max) = u¹ᴰ[σ(y, ·)](z_max)
#     u¹ᴰ[σ] = 1D field of the column σ(z) at that y, normalised to u¹ᴰ(0) = 1
#
#   1D problem for a column σ(z) over a halfspace continuing its deepest value
#     d²E/dz² + k² E = 0,   k² = ω²μ₀ε₀ - iωμ₀σ,   Im k < 0
#     E and dE/dz continuous, only the decaying exp(-ikz) wave in the halfspace
#     TE takes u¹ᴰ = E
#     TM takes u¹ᴰ = H = -(1/iωμ) dE/dz, which solves d/dz(ρ dH/dz) - iωμ H = 0
#     the ω²μ₀ε₀ term is kept only in this 1D problem
#
#   data, fields at each receiver on the air/earth interface
#     Zxy = Ex / Hy (TE),   Zyx = Ey / Hx (TM)
#     ρa = |Z|² / (ωμ₀),   φ = atan(Im Z, Re Z)
#     Zxy lies in the first quadrant, Zyx in the third
#***********************************************************************

using LinearAlgebra
using Printf
using Random
using SparseArrays
using ForwardDiff

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
- Stores the sparse-operator form of the 2D tensor mesh used by the finite-difference (box integration) solver.
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

# sparse identity
spunit(n::Integer) = sparse(1.0I, n, n)

# sparse two-diagonal matrix
spdiag((x1, x2), (d1, d2), m, n) = spdiagm(m, n, d1 => x1, d2 => x2)

# first difference, n cells to n+1 nodes
ddx(n::Integer) = spdiag((-ones(n), ones(n)), (0, 1), n, n + 1)

# node-to-cell average
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

# sparse diagonal
sdiag(values::AbstractVector) = spdiagm(0 => values)

# cell areas as a sparse diagonal
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

# 2D cell-to-node average
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
- Builds the tensor-mesh container used by the 2D finite-difference (box integration) solver.
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
    conductivity::AbstractVector{<:Real},
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

    layers = zeros(Complex{eltype(conductivity)}, 2, n_layers)
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
    conductivity::AbstractVector{<:Real},
)
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    z_nodes = [0.0; cumsum(z_lengths)]
    σ2d = reshape(conductivity, n_y, n_z)'
    boundary = zeros(Complex{eltype(conductivity)}, 2 * (n_y + n_z))

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
    conductivity::AbstractVector{<:Real},
)
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    z_nodes = [0.0; cumsum(z_lengths)]
    σ2d = reshape(conductivity, n_y, n_z)'
    boundary = zeros(Complex{eltype(conductivity)}, 2 * (n_y + n_z))

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
    sigma_row::AbstractVector{<:Real},
    electric_pair::AbstractMatrix{<:Complex},
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

    CT = promote_type(Complex{eltype(sigma_row)}, eltype(electric_pair))
    Hysurface = zeros(CT, n_y + 1)
    Hysurface[2:end-1] = Hyhalf .- (∂Hz∂y .- σavg .* Equarter) .* (0.5 * first_cell_thickness)
    Hysurface[1] = Hysurface[2]
    Hysurface[end] = Hysurface[end - 1]

    Erec = zeros(CT, n_receivers)
    Hrec = zeros(CT, n_receivers)
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
    sigma_row::AbstractVector{<:Real},
    magnetic_pair::AbstractMatrix{<:Complex},
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

    CT = promote_type(Complex{eltype(sigma_row)}, eltype(magnetic_pair))
    Eysurface = zeros(CT, n_y + 1)
    Eysurface[2:end-1] = Eyhalf .- (∂Ez∂y .+ 1im * ω * μ₀_2D .* Hquarter) .* (0.5 * first_cell_thickness)
    Eysurface[1] = Eysurface[2]
    Eysurface[end] = Eysurface[end - 1]

    Erec = zeros(CT, n_receivers)
    Hrec = zeros(CT, n_receivers)
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

# fold phases into [0, 90] degrees
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
    data_type::String;
    return_state::Bool = false,
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
    factor = lu(Aii)
    interior = factor \ rhs

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
    u = vec(copy(transpose(field)))
    return_state && return response, u, (; factor, Aio, u, z_index, frequency)
    response, u
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
    return_state::Bool = false,
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
    factor = lu(Aii)
    interior = factor \ rhs

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
    u = vec(copy(transpose(field)))
    return_state && return response, u, (; factor, Aio, u, z_index, frequency)
    response, u
end

"""
    _assemble_mt2d_system(mesh, resistivity)

Inputs:
- `mesh`: 2D MT mesh.
- `resistivity`: Cell resistivity model.

Output:
- Tuple with the tensor mesh, TE coefficients, TM coefficients, and receiver locations.

Description:
- Assembles the sparse finite-difference (box integration) operators needed for the 2D forward solve.
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

function _mt2d_forward_cache(mesh, resistivity; mode = :TETM, cache_fields::Bool = true)
    mode in (:TE, :TM, :TETM) || throw(ArgumentError("mode must be :TE, :TM, or :TETM"))
    size(resistivity) == (length(mesh.z_cell_sizes), length(mesh.y_cell_sizes)) ||
        throw(DimensionMismatch("resistivity must have shape (nz, ny)"))
    0 <= mesh.n_air_cells < size(resistivity, 1) || throw(ArgumentError("invalid air cell count"))
    all(x -> isfinite(x) && x > 0, resistivity[mesh.n_air_cells+1:end, :]) ||
        throw(ArgumentError("earth resistivities must be finite and positive"))
    all(x -> isfinite(x) && x > 0, mesh.frequencies) || throw(ArgumentError("frequencies must be finite and positive"))
    all(y -> first(mesh.y_nodes) <= y < last(mesh.y_nodes), mesh.receiver_positions) ||
        throw(ArgumentError("receivers must lie in [first(y_nodes), last(y_nodes))"))
    t, te, tm, locations = _assemble_mt2d_system(mesh, resistivity)
    inside, outside = get_boundary_index(t.grid_size...)
    nf, nr = length(mesh.frequencies), length(mesh.receiver_positions)
    arrays = (rho_xy = zeros(nf, nr), phase_xy = zeros(nf, nr), z_xy = zeros(ComplexF64, nf, nr),
              rho_yx = zeros(nf, nr), phase_yx = zeros(nf, nr), z_yx = zeros(ComplexF64, nf, nr))
    states = []
    if nr > 0
        for (i, f) in enumerate(mesh.frequencies), pol in (:TE, :TM)
            mode in (pol, :TETM) || continue
            solver, coeff = pol == :TE ? (solve_mt2d_te, te) : (solve_mt2d_tm, tm)
            data, _, state = solver(f, t, coeff, locations, "Impedance"; return_state = true)
            Z = complex.(data[:, 1], data[:, 2])
            rho, phase, z = pol == :TE ? (arrays.rho_xy, arrays.phase_xy, arrays.z_xy) :
                                       (arrays.rho_yx, arrays.phase_yx, arrays.z_yx)
            z[i, :] = Z
            rho[i, :] = abs2.(Z) ./ (2π * f * μ₀_2D)
            phase[i, :] = rad2deg.(angle.(Z))
            cache_fields && push!(states, (; state..., pol, i))
        end
    end
    response = MT2DResponse(; frequencies = mesh.frequencies, periods = 1 ./ mesh.frequencies,
                            receivers = mesh.receiver_positions, arrays...)
    response, (; mesh, t, locations, inside, outside, states, response)
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
    response, _ = _mt2d_forward_cache(mesh, resistivity; mode, cache_fields = false)
    response
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
Forward2D(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real}; kwargs...) = run_mt2d_forward(mesh, resistivity; kwargs...)

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
        println(io, "# Generated by MTGeophysics.jl/src/Fwd2D.jl")
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


#---------- fréchet derivatives ----------


const _MT2D_FIELDS = (:rho_xy, :phase_xy, :z_xy, :rho_yx, :phase_yx, :z_yx)

# forwarddiff only on the small model-dependent boundary and receiver kernels
# the sparse lu is differentiated implicitly, A δu = -δA u with δu = δb on the boundary
_mt2d_seed(x::Real, δx::Real) = ForwardDiff.Dual{Nothing}(x, δx)
_mt2d_seed(x::Complex, δx::Complex) = complex(_mt2d_seed(real(x), real(δx)), _mt2d_seed(imag(x), imag(δx)))
_mt2d_partial(x::ForwardDiff.Dual) = ForwardDiff.partials(x)[1]
_mt2d_partial(x::Complex) = complex(_mt2d_partial(real(x)), _mt2d_partial(imag(x)))

function _mt2d_parameter_scale(mesh, rho, parameterization)
    parameterization in (:resistivity, :log_resistivity, :log10_resistivity) ||
        throw(ArgumentError("parameterization must be :resistivity, :log_resistivity, or :log10_resistivity"))
    scale = parameterization == :resistivity ? ones(size(rho)) :
            parameterization == :log_resistivity ? Float64.(rho) : log(10.0) .* rho
    scale[1:mesh.n_air_cells, :] .= 0
    scale
end

function _mt2d_surface_input(c, s)
    ny = c.t.grid_size[1]
    row = (s.z_index-1)*ny+1:s.z_index*ny
    nodes = (s.z_index-1)*(ny+1)+1:(s.z_index+1)*(ny+1)
    vcat(c.t.conductivity[row], real.(s.u[nodes]), imag.(s.u[nodes])), row, nodes
end

function _mt2d_sample(c, s, x)
    ny = c.t.grid_size[1]
    n = 2*(ny+1)
    σ = x[1:ny]
    pair = reshape(complex.(x[ny+1:ny+n], x[ny+n+1:ny+2n]), ny+1, 2)
    sampler = s.pol == :TE ? compute_fields_at_receivers_te : compute_fields_at_receivers_tm
    E, H = sampler(2π*s.frequency, c.locations, [0.0; cumsum(c.t.y_lengths)],
                   c.t.z_lengths[s.z_index], σ, pair)
    E ./ H
end

function _mt2d_frechet(c, δρ)
    σ = c.t.conductivity
    δσ = -σ.^2 .* vec(permutedims(δρ))
    δσ[1:c.mesh.n_air_cells*c.t.grid_size[1]] .= 0
    Mn = c.t.average_cell_to_node * c.t.face
    Mf = c.t.average_cell_to_face * c.t.face
    D = c.t.gradient
    out = NamedTuple{_MT2D_FIELDS}(Tuple(zeros(eltype(getproperty(c.response, k)), size(getproperty(c.response, k))) for k in _MT2D_FIELDS))
    for s in c.states
        boundary = s.pol == :TE ? get_boundary_mt2d_te : get_boundary_mt2d_tm
        δb = _mt2d_partial.(boundary(s.frequency, c.t.y_lengths, c.t.z_lengths, _mt2d_seed.(σ, δσ)))
        δAu = s.pol == :TE ? 1im*2π*s.frequency .* (Mn * δσ) .* s.u :
                             D' * ((Mf * (-δσ ./ σ.^2)) .* (D*s.u))
        δu = zeros(ComplexF64, length(s.u))
        δu[c.outside] = δb
        δu[c.inside] = s.factor \ (-δAu[c.inside] - s.Aio * δb)
        x, row, nodes = _mt2d_surface_input(c, s)
        δx = vcat(δσ[row], real.(δu[nodes]), imag.(δu[nodes]))
        δz = _mt2d_partial.(_mt2d_sample(c, s, _mt2d_seed.(x, δx)))
        rk, pk, zk = s.pol == :TE ? (:rho_xy, :phase_xy, :z_xy) : (:rho_yx, :phase_yx, :z_yx)
        z = getproperty(c.response, zk)[s.i, :]
        getproperty(out, zk)[s.i, :] = δz
        getproperty(out, rk)[s.i, :] = 2 .* real.(conj.(z) .* δz) ./ (2π*s.frequency*μ₀_2D)
        getproperty(out, pk)[s.i, :] = rad2deg.(imag.(δz ./ z))
    end
    out
end

# Gᵗ through the 1D boundary columns, one column at a time, so no dense boundary
# fréchet block is formed and the reverse pass stays local to each profile
function _mt2d_boundary_frechet_transpose(c, s, b̂)
    ny, nz = c.t.grid_size
    σ = reshape(c.t.conductivity, ny, nz)
    g = zeros(ny, nz)
    zn = [0.0; cumsum(c.t.z_lengths)]
    function profile(x)
        if s.pol == :TE
            E = mt1d_boundary_field(s.frequency, x, zn)
            return vec(E ./ E[1])
        end
        _, H = mt1d_boundary_field(s.frequency, x, zn; return_magnetic = true)
        vec(H ./ H[1])
    end
    for (iy, inds) in ((1, ny+2:ny+nz+1), (ny, ny+nz+2:ny+2nz+1))
        weights = b̂[inds]
        g[iy, :] .+= ForwardDiff.gradient(x -> real(dot(weights, profile(x)[2:end])), σ[iy, :])
    end
    for iy in 2:ny
        w = c.t.y_lengths[iy-1] / (c.t.y_lengths[iy-1] + c.t.y_lengths[iy])
        x = w .* σ[iy-1, :] .+ (1-w) .* σ[iy, :]
        weight = b̂[ny+2nz+iy]
        x̂ = ForwardDiff.gradient(x -> real(conj(weight)*profile(x)[end]), x)
        g[iy-1, :] .+= w .* x̂
        g[iy, :] .+= (1-w) .* x̂
    end
    vec(g)
end

function _mt2d_dual(δd̂, key, template)
    !hasproperty(δd̂, key) && return zero(template)
    value = getproperty(δd̂, key)
    size(value) == size(template) || throw(DimensionMismatch("δd̂ $key has the wrong shape"))
    value
end

_mt2d_ops(c) = (Mn = c.t.average_cell_to_node * c.t.face, Mf = c.t.average_cell_to_face * c.t.face, D = c.t.gradient)

# Gᵗ for one frequency and mode: adds ∂ real(dot(ẑ, z)) / ∂σ into σ̂,
# ẑ = impedance part of δd̂ at the receivers, λ = adjoint field
function _mt2d_state_frechet_transpose!(σ̂, c, s, ẑ, ops)
    σ = c.t.conductivity
    x, row, nodes = _mt2d_surface_input(c, s)
    x̂ = ForwardDiff.gradient(x -> real(dot(ẑ, _mt2d_sample(c, s, x))), x)
    ny = c.t.grid_size[1]
    n = length(nodes)
    û = zeros(ComplexF64, length(s.u))
    û[nodes] = complex.(x̂[ny+1:ny+n], x̂[ny+n+1:ny+2n])
    σ̂[row] .+= x̂[1:ny]
    λ = zeros(ComplexF64, length(s.u))
    λ[c.inside] = s.factor' \ û[c.inside]
    b̂ = û[c.outside] - s.Aio' * λ[c.inside]
    if s.pol == :TE
        σ̂ .+= ops.Mn' * real.(-1im*2π*s.frequency .* conj.(λ) .* s.u)
    else
        σ̂ .+= (ops.Mf' * real.(conj.(ops.D*λ) .* (ops.D*s.u))) ./ σ.^2
    end
    σ̂ .+= _mt2d_boundary_frechet_transpose(c, s, b̂)
    σ̂
end

# σ̂ to model-shaped ρ̂ = -σ² σ̂, air zeroed
function _mt2d_dual_to_rho(c, σ̂)
    ρ̂ = permutedims(reshape(-c.t.conductivity.^2 .* σ̂, c.t.grid_size...))
    ρ̂[1:c.mesh.n_air_cells, :] .= 0
    ρ̂
end

function _mt2d_frechet_transpose(c, δd̂)
    σ̂ = zeros(length(c.t.conductivity))
    ops = _mt2d_ops(c)
    d̂ = NamedTuple{_MT2D_FIELDS}(Tuple(_mt2d_dual(δd̂, k, getproperty(c.response, k)) for k in _MT2D_FIELDS))
    for s in c.states
        rk, pk, zk = s.pol == :TE ? (:rho_xy, :phase_xy, :z_xy) : (:rho_yx, :phase_yx, :z_yx)
        z = getproperty(c.response, zk)[s.i, :]
        ẑ = getproperty(d̂, zk)[s.i, :] .+
               (2/(2π*s.frequency*μ₀_2D)) .* getproperty(d̂, rk)[s.i, :] .* z .+
               (180/π) .* getproperty(d̂, pk)[s.i, :] .* (1im ./ conj.(z))
        _mt2d_state_frechet_transpose!(σ̂, c, s, ẑ, ops)
    end
    _mt2d_dual_to_rho(c, σ̂)
end

# rows of G by transpose solves, each (key, index, weight) row gives
# ∂ real(conj(weight) z[key][index]) / ∂ρ over the whole model, one adjoint
# solve per row in its own frequency and mode, rows × model cells, (z, y) order
function _mt2d_frechet_rows(c, rows)
    ops = _mt2d_ops(c)
    nr = size(c.response.z_xy, 2)
    state = Dict((s.pol == :TE ? :z_xy : :z_yx, s.i) => s for s in c.states)
    out = zeros(length(rows), prod(c.t.grid_size))
    σ̂ = zeros(length(c.t.conductivity))
    ẑ = zeros(ComplexF64, nr)
    for (k, row) in enumerate(rows)
        f, r = Tuple(CartesianIndices(size(c.response.z_xy))[row.index])
        fill!(σ̂, 0); fill!(ẑ, 0)
        ẑ[r] = row.weight
        _mt2d_state_frechet_transpose!(σ̂, c, state[(row.key, f)], ẑ, ops)
        out[k, :] = vec(_mt2d_dual_to_rho(c, σ̂))
    end
    out
end

"""
    ApplyFrechet2D(mesh, ρ, δm; mode=:TETM, parameterization=:resistivity)

Tangent linear application δd = G δm without forming G, with G = ∂g/∂m the Fréchet
derivative of the forward relation d = g(m) (Tarantola, 2005).
- `ρ`: model in ohm metres, air included
- `δm`: model-shaped perturbation in `parameterization` coordinates,
  `:resistivity`, `:log_resistivity` (natural log), or `:log10_resistivity`
- returns δd as a named tuple of `rho_xy`, `phase_xy`, `z_xy`, `rho_yx`, `phase_yx`, `z_yx`,
  each `(frequency, receiver)`, phases in degrees; air cells are fixed
"""
function ApplyFrechet2D(mesh::MT2DMesh, ρ::AbstractMatrix{<:Real}, δm::AbstractMatrix{<:Real};
                  mode::Symbol = :TETM, parameterization::Symbol = :resistivity)
    size(δm) == size(ρ) || throw(DimensionMismatch("δm must match resistivity"))
    _, c = _mt2d_forward_cache(mesh, ρ; mode)
    _mt2d_frechet(c, δm .* _mt2d_parameter_scale(mesh, ρ, parameterization))
end

"""
    ApplyFrechetTranspose2D(mesh, ρ, δd̂; mode=:TETM, parameterization=:resistivity)

Transpose application δm̂ = Gᵗ δd̂ (Tarantola, 2005), one adjoint solve per frequency and mode.
- `δd̂`: dual data vector, a named tuple with any subset of the six response fields,
  missing fields are zero; impedances pair as `real(dot(δd̂, δz))`
- returns δm̂, model-shaped `(z, y)`, in `parameterization` coordinates, so that
  `sum(δm̂ .* δm)` equals the pairing of δd̂ with `ApplyFrechet2D(mesh, ρ, δm)`
"""
function ApplyFrechetTranspose2D(mesh::MT2DMesh, ρ::AbstractMatrix{<:Real}, δd̂;
                  mode::Symbol = :TETM, parameterization::Symbol = :resistivity)
    _, c = _mt2d_forward_cache(mesh, ρ; mode)
    _mt2d_frechet_transpose(c, δd̂) .* _mt2d_parameter_scale(mesh, ρ, parameterization)
end

"""
    FrechetDerivative2D(mesh, ρ; mode=:TETM, parameterization=:resistivity, active_cells=nothing)

Explicit Fréchet derivative G, Gⁱ_α = ∂gⁱ/∂mᵅ (Tarantola, 2005), one tangent linear solve per
column with the forward lu factors reused.
- returns `(response, cells, parameterization, rho_xy, phase_xy, z_xy, rho_yx, phase_yx, z_yx)`,
  each block `(nf * nr, length(cells))`; rows follow `vec(response.field)`, frequency fastest;
  columns follow `cells`, Cartesian indices in model `(z, y)` order
- `active_cells`: `nothing` for all earth cells, a model-shaped Bool mask, or Cartesian/linear
  indices; requested air columns are zero
- for large models use `ApplyFrechet2D` and `ApplyFrechetTranspose2D` instead
"""
function FrechetDerivative2D(mesh::MT2DMesh, ρ::AbstractMatrix{<:Real};
                           mode::Symbol = :TETM, parameterization::Symbol = :resistivity,
                           active_cells = nothing)
    response, c = _mt2d_forward_cache(mesh, ρ; mode)
    scale = _mt2d_parameter_scale(mesh, ρ, parameterization)
    allcells = CartesianIndices(ρ)
    cells = if active_cells === nothing
        [i for i in allcells if i[1] > mesh.n_air_cells]
    elseif active_cells isa AbstractArray{Bool}
        size(active_cells) == size(ρ) || throw(DimensionMismatch("active mask must match resistivity"))
        findall(active_cells)
    else
        [i isa CartesianIndex{2} ? i : allcells[i] for i in active_cells]
    end
    all(i -> checkbounds(Bool, ρ, i), cells) || throw(BoundsError(ρ, cells))
    G = NamedTuple{_MT2D_FIELDS}(Tuple(zeros(eltype(getproperty(response, k)), length(response.z_xy), length(cells)) for k in _MT2D_FIELDS))
    δm = zeros(size(ρ))
    for (j, cell) in enumerate(cells)
        δm[cell] = scale[cell]
        δd = _mt2d_frechet(c, δm)
        for k in _MT2D_FIELDS
            getproperty(G, k)[:, j] = vec(getproperty(δd, k))
        end
        δm[cell] = 0
    end
    (; response, cells, parameterization, G...)
end


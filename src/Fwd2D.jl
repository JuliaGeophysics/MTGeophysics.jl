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
#   domain       Ω = [y₀, y₁] × [0, z_max], air included with σ_air = 1/ρ_air
#   topography   σ = σ_air above the ground surface z = s(y), which may vary with y
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
#   data, fields at each receiver on the ground surface, z = s(y)
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

TE and TM response of a profile, each field `(frequency, receiver)`: apparent
resistivity (ohm m), phase (degrees) and impedance (ohm), exp(+iωt).
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

2D impedance data, each field `(frequency, site)`: TE = ZXY and TM = ZYX with their
errors (ohm, exp(+iωt)), derived apparent resistivity and phase, site names, local
positions (`receivers` along the profile, `x_positions`, `z_positions`), WGS84
`latitudes`/`longitudes` and the survey `origin` (lat, lon).
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
    latitudes::Vector{Float64} = Float64[]
    longitudes::Vector{Float64} = Float64[]
    origin::Vector{Float64} = [0.0, 0.0]
end

# chi2 over real data (Re and Im count separately), rms = sqrt(chi2 / count)
Base.@kwdef struct FitSummary2D
    chi2::Float64
    rms::Float64
    count::Int
end

#---------- tensor mesh operators ----------

# sparse operators of the box-integration scheme; cells y fastest
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

# interior/interior and interior/boundary blocks, A = real + iω imag
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

# cell-to-node average, end nodes take their one cell
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

# inverse edge lengths, y edges then z edges
function mesh_geo_edge_inv_2d(d1::Vector, d2::Vector)
    n1 = length(d1)
    n2 = length(d2)
    left = kron(spunit(n2 + 1), sdiag(1.0 ./ d1))
    right = kron(sdiag(1.0 ./ d2), spunit(n1 + 1))
    blockdiag(left, right)
end

# nodes to edges
function nodal_gradient_2d(d1::Vector, d2::Vector)
    n1 = length(d1)
    n2 = length(d2)
    g1 = kron(spunit(n2 + 1), ddx(n1))
    g2 = kron(ddx(n2), spunit(n1 + 1))
    mesh_geo_edge_inv_2d(d1, d2) * [g1; g2]
end

# 2D cell-to-node average
average_cell_to_node_2d(grid_size::Vector{Int}) = kron(avcn(grid_size[2]), avcn(grid_size[1]))

# cells to edges, z-directed faces then y-directed
function average_cell_to_face_2d(grid_size::Vector{Int})
    face_y = kron(spunit(grid_size[2]), avcn(grid_size[1]))
    face_z = kron(avcn(grid_size[2]), spunit(grid_size[1]))
    [face_z; face_y]
end

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

function setup_tensor_mesh_2d!(mesh::TensorMesh2D)
    mesh.face = mesh_geo_face_2d(mesh.y_lengths, mesh.z_lengths)
    mesh.gradient = nodal_gradient_2d(mesh.y_lengths, mesh.z_lengths)
    mesh.average_cell_to_node = average_cell_to_node_2d(mesh.grid_size)
    mesh.average_cell_to_face = average_cell_to_face_2d(mesh.grid_size)
    mesh.setup = true
    mesh
end

#---------- boundary conditions ----------

# 1D layered-earth E (and H) at the nodes z_nodes, E = 1 at the top
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

# interior nodes, and boundary nodes ordered top, left, right, bottom
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

# TE boundary values from the 1D E field of the side columns and of each bottom node
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

# TM boundary values, the same from the 1D H field
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

#---------- receiver fields ----------

# TE: Ex and Hy on the nodes of a surface row, from its two node rows
function _mt2d_surface_fields_te(
    ω::Float64,
    y_nodes::Vector{Float64},
    first_cell_thickness::Float64,
    sigma_row::AbstractVector{<:Real},
    electric_pair::AbstractMatrix{<:Complex},
)
    y_lengths = diff(y_nodes)
    n_y = length(y_lengths)
    μ = μ₀_2D .* ones(n_y)

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
    E, Hysurface
end

# linear interpolation of nodal surface fields to the receivers, unnormalised; only E / H is used
function _mt2d_interpolate_receivers(receiver_locations, y_nodes, E, H)
    n_receivers = size(receiver_locations, 1)
    Erec = zeros(eltype(E), n_receivers)
    Hrec = zeros(eltype(H), n_receivers)
    for i in 1:n_receivers
        y = receiver_locations[i, 1]
        node = findfirst(v -> v > y, y_nodes)
        Δy1 = y - y_nodes[node - 1]
        Δy2 = y_nodes[node] - y
        Erec[i] = E[node - 1] * Δy2 + E[node] * Δy1
        Hrec[i] = H[node - 1] * Δy2 + H[node] * Δy1
    end
    Erec, Hrec
end

# TE: Ex and Hy at the receivers of one surface row
function compute_fields_at_receivers_te(ω::Float64, receiver_locations::Matrix{Float64}, y_nodes::Vector{Float64},
                                        first_cell_thickness::Float64, sigma_row::AbstractVector{<:Real},
                                        electric_pair::AbstractMatrix{<:Complex})
    E, H = _mt2d_surface_fields_te(ω, y_nodes, first_cell_thickness, sigma_row, electric_pair)
    _mt2d_interpolate_receivers(receiver_locations, y_nodes, E, H)
end

# TM: Ey and Hx on the nodes of a surface row, from its two node rows
function _mt2d_surface_fields_tm(
    ω::Float64,
    y_nodes::Vector{Float64},
    first_cell_thickness::Float64,
    sigma_row::AbstractVector{<:Real},
    magnetic_pair::AbstractMatrix{<:Complex},
)
    y_lengths = diff(y_nodes)
    n_y = length(y_lengths)

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

    # tread-centre Ey of every cell from its own ρ alone: next to a topographic step a surface node
    # averages in the neighbouring air cell, whose ρ makes the nodal Ey mesh dependent
    Hmid0 = (magnetic_pair[1:end-1, 1] .+ magnetic_pair[2:end, 1]) ./ 2
    Hmid1 = (magnetic_pair[1:end-1, 2] .+ magnetic_pair[2:end, 2]) ./ 2
    Eycell = (Hmid1 .- Hmid0) ./ (first_cell_thickness .* sigma_row) .-
             1im * ω * μ₀_2D .* (0.75 .* Hmid0 .+ 0.25 .* Hmid1) .* (0.5 * first_cell_thickness)
    Eysurface, CT.(H), CT.(Eycell)
end

# TM: Ey and Hx at the receivers of one surface row
function compute_fields_at_receivers_tm(ω::Float64, receiver_locations::Matrix{Float64}, y_nodes::Vector{Float64},
                                        first_cell_thickness::Float64, sigma_row::AbstractVector{<:Real},
                                        magnetic_pair::AbstractMatrix{<:Complex})
    E, H = _mt2d_surface_fields_tm(ω, y_nodes, first_cell_thickness, sigma_row, magnetic_pair)
    _mt2d_interpolate_receivers(receiver_locations, y_nodes, E, H)
end

# fold phases into [0, 90] degrees
_phase_fold_to_0_90(phases::AbstractArray) = map(ϕ -> begin
    folded = ϕ < 0 ? ϕ + 180 : ϕ
    folded > 90 ? 180 - folded : folded
end, phases)

# node row nearest the receiver depth; model files round cell sizes, so no exact match
function _mt2d_surface_row(z_lengths, z_origin, surface_depth)
    nodes = [0.0; cumsum(z_lengths)] .- z_origin
    offset, index = findmin(abs.(nodes .- surface_depth))
    offset <= 1e-6 * last(nodes) || error("no mesh node row at the receiver depth $(surface_depth) m")
    index
end

#---------- solve ----------

# one frequency and mode: impedance at the receivers and the state the Fréchet code reuses
function _solve_mt2d(pol::Symbol, frequency::Float64, mesh::TensorMesh2D, coefficients::CoeffMat,
                     receiver_locations::Matrix{Float64}, plan)
    boundary_values = pol == :TE ? get_boundary_mt2d_te : get_boundary_mt2d_tm
    y_lengths = mesh.y_lengths
    z_lengths = mesh.z_lengths
    conductivity = mesh.conductivity
    y_nodes = [0.0; cumsum(y_lengths)] .- mesh.origin[1]
    n_y = length(y_lengths)
    n_z = length(z_lengths)
    ω = 2π * frequency

    Aii = coefficients.real_ii + 1im * ω * coefficients.imag_ii
    Aio = coefficients.real_io + 1im * ω * coefficients.imag_io
    boundary = boundary_values(frequency, y_lengths, z_lengths, conductivity)
    rhs = -Aio * boundary
    factor = lu(Aii)
    interior = factor \ rhs

    field = zeros(ComplexF64, n_z + 1, n_y + 1)
    field[1, :] = boundary[1:n_y+1]
    field[2:end, 1] = boundary[n_y+2:n_y+n_z+1]
    field[2:end, end] = boundary[n_y+n_z+2:n_y+2*n_z+1]
    field[end, 2:end-1] = boundary[n_y+2*n_z+2:end]
    field[2:end-1, 2:end-1] = copy(transpose(reshape(interior, n_y - 1, n_z - 1)))

    surface = map(plan.groups) do g
        pair = copy(transpose(field[g.z_index:g.z_index+1, :]))
        σrow = conductivity[(g.z_index - 1) * n_y + 1:g.z_index * n_y]
        _mt2d_surface_fields(pol, ω, y_nodes, z_lengths[g.z_index], σrow, pair)
    end

    u = vec(copy(transpose(field)))
    _mt2d_station_impedance(pol, plan, surface, receiver_locations, y_nodes), (; factor, Aio, u, frequency)
end

_mt2d_surface_fields(pol, ω, y_nodes, dz, σrow, pair) =
    pol == :TE ? _mt2d_surface_fields_te(ω, y_nodes, dz, σrow, pair) :
                 _mt2d_surface_fields_tm(ω, y_nodes, dz, σrow, pair)

# where each receiver reads its fields: the surface node row of its column, and, next to a
# topographic step, a dipole window of the columns within half the dipole length each side
# (its own column alone for a dipole shorter than a cell) whose tread-centre TM Ey, from each
# cell's own ρ, is averaged by width, as an electric dipole integrates E; on a staircase the
# point value swings between the convex and concave step corners, and a corner node's Ey takes
# in the air beside it (checked on the trapezoidal hill in test/TestTopography2D.jl)
function _mt2d_station_plan(mesh::MT2DMesh, t::TensorMesh2D, locations::Matrix{Float64})
    ny = length(mesh.y_cell_sizes)
    topo = mt2d_topo_air(mesh)
    surface_row(k) = mesh.n_air_cells + topo[k] + 1
    columns = mt2d_receiver_columns(mesh)
    rows = [_mt2d_surface_row(t.z_lengths, t.origin[2], locations[i, 2]) for i in axes(locations, 1)]
    all(i -> rows[i] == surface_row(columns[i]), eachindex(rows)) || error("receiver depths do not match the ground")
    windows = map(columns) do c
        half = round(Int, mesh.dipole_length / 2 / mesh.y_cell_sizes[c])
        cols = max(1, c - half):min(ny, c + half)
        flat = all(k -> topo[k] == topo[c], union(cols, max(1, c - 1):min(ny, c + 1)))
        flat ? nothing : collect(cols)
    end
    needed = sort(unique(vcat(rows, [surface_row(k) for w in windows if w !== nothing for k in w])))
    index = Dict(r => i for (i, r) in enumerate(needed))
    groups = [(z_index = r,) for r in needed]
    stations = [(group = index[rows[i]],
                 window = windows[i] === nothing ? nothing :
                          [(group = index[surface_row(k)], column = k, width = mesh.y_cell_sizes[k]) for k in windows[i]])
                for i in eachindex(rows)]
    (; groups, stations)
end

# impedance at every receiver from the nodal surface fields of each group
function _mt2d_station_impedance(pol, plan, surface, locations, y_nodes)
    T = promote_type(eltype(surface[1][1]), eltype(surface[1][2]))
    Z = zeros(T, length(plan.stations))
    for (i, p) in enumerate(plan.stations)
        Es, Hs = surface[p.group]
        y = locations[i, 1]
        node = findfirst(v -> v > y, y_nodes)
        Δy1 = y - y_nodes[node - 1]
        Δy2 = y_nodes[node] - y
        E = Es[node - 1] * Δy2 + Es[node] * Δy1
        H = Hs[node - 1] * Δy2 + Hs[node] * Δy1
        if pol == :TM && p.window !== nothing
            Ew = sum(w.width * surface[w.group][3][w.column] for w in p.window)
            E = (Δy1 + Δy2) * Ew / sum(w.width for w in p.window)
        end
        Z[i] = E / H
    end
    Z
end

# TE (Ex) and TM (Hx) systems on the tensor mesh, and the receiver locations
function _assemble_mt2d_system(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real})
    n_z = length(mesh.z_cell_sizes)
    n_y = length(mesh.y_cell_sizes)
    size(resistivity) == (n_z, n_y) || error("resistivity must be size ($(n_z), $(n_y))")

    σ = 1.0 ./ Matrix{Float64}(resistivity)
    σ[mt2d_air_mask(mesh)] .= 1 / mesh.air_resistivity

    tensor_mesh = TensorMesh2D(
        mesh.y_cell_sizes,
        mesh.z_cell_sizes;
        origin = [0.0, 0.0],
        conductivity = vec(copy(transpose(σ))),
    )
    setup_tensor_mesh_2d!(tensor_mesh)

    receiver_y = mesh.receiver_positions .- first(mesh.y_nodes)
    receiver_z = mesh.z_nodes[mesh.n_air_cells+1] - first(mesh.z_nodes) .+ mt2d_receiver_depths(mesh)
    receiver_locations = hcat(receiver_y, receiver_z)

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

# response plus, with cache_fields, the per-frequency states for Fréchet derivatives
function _mt2d_forward_cache(mesh, resistivity; mode = :TETM, cache_fields::Bool = true)
    mesh.dimension == 1 && return _mt1d_forward_cache(mesh, resistivity; mode, cache_fields)
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
    plan = _mt2d_station_plan(mesh, t, locations)
    air = mt2d_air_mask(mesh)
    nf, nr = length(mesh.frequencies), length(mesh.receiver_positions)
    arrays = (rho_xy = zeros(nf, nr), phase_xy = zeros(nf, nr), z_xy = zeros(ComplexF64, nf, nr),
              rho_yx = zeros(nf, nr), phase_yx = zeros(nf, nr), z_yx = zeros(ComplexF64, nf, nr))
    states = []
    if nr > 0
        for (i, f) in enumerate(mesh.frequencies), pol in (:TE, :TM)
            mode in (pol, :TETM) || continue
            Z, state = _solve_mt2d(pol, f, t, pol == :TE ? te : tm, locations, plan)
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
    response, (; mesh, t, locations, inside, outside, states, response, air, plan)
end

"""
    run_mt2d_forward(mesh, resistivity; mode=:TETM) -> MT2DResponse

TE and/or TM response of `resistivity` (ohm m, `(nz, ny)`, air rows included) at every
frequency and receiver of `mesh`; the fields of a mode that is not solved stay zero.
"""
function run_mt2d_forward(
    mesh::MT2DMesh,
    resistivity::AbstractMatrix{<:Real};
    mode::Symbol = :TETM,
)
    response, _ = _mt2d_forward_cache(mesh, resistivity; mode, cache_fields = false)
    response
end

Forward2D(mesh::MT2DMesh, resistivity::AbstractMatrix{<:Real}; kwargs...) = run_mt2d_forward(mesh, resistivity; kwargs...)

#---------- data ----------

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
    data_from_response2d(response; z_xy_error, z_yx_error, z_xx_error, z_yy_error,
                         impedance_error_fraction=0.05, title, site_names, x_positions,
                         z_positions, latitudes, longitudes, origin) -> DataFile2D

Data from a forward response. Missing errors default to `impedance_error_fraction`
of |Z|, missing sites to `Site001`… at x = z = 0.
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
    latitudes::AbstractVector{<:Real} = Float64[],
    longitudes::AbstractVector{<:Real} = Float64[],
    origin::AbstractVector{<:Real} = [0.0, 0.0],
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
        latitudes = Float64.(latitudes),
        longitudes = Float64.(longitudes),
        origin = Float64.(origin),
    )
    size(data.z_xy_error) == size(data.z_xy) || error("z_xy_error does not match z_xy size")
    size(data.z_yx_error) == size(data.z_yx) || error("z_yx_error does not match z_yx size")
    data
end

# a template (all impedances zero) holds error fractions of the predicted |Z|
function _resolve_forwardsolve2d_errors(template::DataFile2D, response::MT2DResponse)
    if all(iszero, template.z_xy) && all(iszero, template.z_yx)
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
    write_data2d(path, data_or_response; impedance_error_fraction=0.05) -> path

Write 2D impedances as a ModEM Full_Impedance file, ZXY (TE) and ZYX (TM) only,
exp(+iωt), [mV/km]/[nT], through `write_data_modem`, so 1D, 2D and 3D codes read
the same file. Sites keep their lat/lon, local x, y (profile) and z (depth).
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
    nf, ns = size(data.z_xy)
    Z = fill(ComplexF64(NaN, NaN), nf, 4, ns)
    Zerr = fill(ComplexF64(NaN, NaN), nf, 4, ns)
    Z[:, 2, :] .= data.z_xy
    Z[:, 3, :] .= data.z_yx
    Zerr[:, 2, :] .= data.z_xy_error
    Zerr[:, 3, :] .= data.z_yx_error
    coord(v) = isempty(v) ? zeros(ns) : Float64.(v)
    d = make_nan_data()
    d.T, d.f = Float64.(data.periods), Float64.(data.frequencies)
    d.site, d.ns, d.nf = String.(data.site_names), ns, nf
    d.responses, d.nr = ["ZXY", "ZYX"], 2
    d.loc = hcat(coord(data.latitudes), coord(data.longitudes), Float64.(data.z_positions))
    d.x, d.y, d.z = Float64.(data.x_positions), Float64.(data.receivers), Float64.(data.z_positions)
    d.Z, d.Zerr = Z, Zerr
    d.tip = fill(ComplexF64(NaN, NaN), nf, 2, ns)
    d.tiperr = copy(d.tip)
    d.zrot = zeros(nf, ns)
    d.trot = d.zrot
    d.origin = [coord(data.origin)[1:2]; 0.0]
    mkpath(dirname(abspath(path)))
    redirect_stdout(devnull) do
        write_data_modem(path, d; sign = 1, units = "[mV/km]/[nT]", include_tipper = false,
                         description = "2D MT data written by MTGeophysics.jl, ZXY = TE, ZYX = TM")
    end
    String(path)
end

"""
    load_data2d(path) -> DataFile2D

Read ZXY (TE) and ZYX (TM) from a ModEM impedance file through `load_data_modem`,
which converts units and time convention to Ohm and exp(+iωt). Frequencies come out
ascending. Local y is the position along the profile.
"""
function load_data2d(path::AbstractString)
    isfile(path) || error("data file not found: $path")
    d = redirect_stdout(() -> load_data_modem(path), devnull)
    any(isfinite, d.Z[:, 2:3, :]) || error("no ZXY or ZYX impedances in $path")
    order = sortperm(d.f)
    frequencies = d.f[order]
    pick(k) = d.Z[order, k, :]
    err(k) = abs.(d.Zerr[order, k, :])
    z_xy, z_yx = pick(2), pick(3)
    rho_phase(z) = [_impedance_to_rho_phase(z[i, j], frequencies[i]) for i in axes(z, 1), j in axes(z, 2)]
    xy, yx = rho_phase(z_xy), rho_phase(z_yx)
    DataFile2D(
        title = basename(path),
        periods = 1 ./ frequencies,
        frequencies = frequencies,
        site_names = String.(d.site),
        receivers = Float64.(d.y),
        x_positions = Float64.(d.x),
        z_positions = Float64.(d.z),
        z_xy = z_xy, z_xy_error = err(2),
        z_yx = z_yx, z_yx_error = err(3),
        z_xx = pick(1), z_xx_error = err(1),
        z_yy = pick(4), z_yy_error = err(4),
        rho_xy = getfield.(xy, :rho), phase_xy = getfield.(xy, :phase),
        rho_yx = getfield.(yx, :rho), phase_yx = getfield.(yx, :phase),
        path = String(path),
        latitudes = Float64.(d.loc[:, 1]),
        longitudes = Float64.(d.loc[:, 2]),
        origin = Float64.(d.origin[1:2]),
    )
end

"""
    chi2_rms2d(observed, predicted; components=["ZXY", "ZYX"]) -> FitSummary2D

Error-weighted impedance misfit, Re and Im counted as separate data, over the finite
observations with positive errors. Also takes two data file paths.
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
        key, errkey = component == "ZXY" ? (:z_xy, :z_xy_error) :
                      component == "ZYX" ? (:z_yx, :z_yx_error) : (nothing, nothing)
        key === nothing && continue
        zobs, zpred, err = getproperty(observed, key), getproperty(predicted, key), getproperty(observed, errkey)
        for index in eachindex(zobs)
            zo, zp, σ = zobs[index], zpred[index], err[index]
            if isfinite(zo) && isfinite(zp) && isfinite(σ) && σ > 0
                χ² += ((real(zp) - real(zo)) / σ)^2 + ((imag(zp) - imag(zo)) / σ)^2
                count += 2
            end
        end
    end

    FitSummary2D(
        chi2 = χ²,
        rms = count > 0 ? sqrt(χ² / count) : NaN,
        count = count,
    )
end

chi2_rms2d(observed_path::AbstractString, predicted_path::AbstractString; components::AbstractVector{<:AbstractString} = ["ZXY", "ZYX"]) =
    chi2_rms2d(load_data2d(observed_path), load_data2d(predicted_path); components = components)

# zero impedances with error fractions at the mesh receivers, sites Site001…
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

write_mt2d_data_template(path::AbstractString, mesh::MT2DMesh; impedance_error_fraction::Real = 0.05) =
    write_data2d(path, build_mt2d_data_template(mesh; impedance_error_fraction))

# gaussian noise on Re and Im, each with half the impedance error variance
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

#---------- file workflow ----------

"""
    ForwardSolve2D(model_path, data_path, fwd_path; output_path, add_noise=false, rng_seed=20260308) -> path

ModEM-style forward run: a ModEM-layout model, a data file giving the sites, periods
and errors, and `fwd.ctrl` giving the mode and the air. Writes the predicted data,
by default `data.pred` next to the data file, and with `Write Frechet derivative : yes`
also G as `<output>.frechet`. A data file whose impedances are all zero is a template:
its errors are then fractions of the predicted |Z|.
"""
function ForwardSolve2D(
    model_path::AbstractString,
    data_path::AbstractString,
    fwd_path::AbstractString;
    output_path::AbstractString = joinpath(dirname(abspath(data_path)), "data.pred"),
    add_noise::Bool = false,
    rng_seed::Integer = 20260308,
)
    template = load_data2d(data_path)
    fwd = ReadFwdCtrl2D(fwd_path)
    mesh, ρ = Mesh2DFromInputs(ReadModel2D(model_path), template, fwd)
    response = run_mt2d_forward(mesh, ρ; mode = fwd.mode)
    errors = _resolve_forwardsolve2d_errors(template, response)
    predicted = data_from_response2d(response; errors...,
        site_names = template.site_names, x_positions = template.x_positions, z_positions = template.z_positions,
        latitudes = template.latitudes, longitudes = template.longitudes, origin = template.origin)
    fwd.mode == :TE && (predicted.z_yx .= NaN)
    fwd.mode == :TM && (predicted.z_xy .= NaN)
    written = write_data2d(output_path, add_noise ? _apply_mt2d_noise(predicted; rng_seed) : predicted)
    fwd.write_frechet && WriteFrechet2D(splitext(written)[1] * ".frechet", mesh, ρ, predicted; mode = fwd.mode)
    written
end

"""
    WriteFrechet2D(path, mesh, ρ, data; mode=:TETM) -> path

Write G = ∂d/∂m for the impedances in `data`, in log10 resistivity of the earth cells.
One row per real datum (Re and Im of each ZXY/ZYX data line), impedance in
[mV/km]/[nT], columns in model order with z fastest from the model top.
"""
function WriteFrechet2D(path::AbstractString, mesh::MT2DMesh, ρ::AbstractMatrix{<:Real}, data::DataFile2D;
                        mode::Symbol = :TETM)
    # every earth-model cell is a column, topographic air ones are zero
    G = FrechetDerivative2D(mesh, ρ; mode, parameterization = :log10_resistivity,
                            active_cells = [i for i in CartesianIndices(ρ) if i[1] > mesh.n_air_cells])
    scale = 1 / (μ₀_2D * 1000)
    nf, ns = size(data.z_xy)
    nz, ny = length(mesh.z_cell_sizes) - mesh.n_air_cells, length(mesh.y_cell_sizes)
    comps = [(k, c) for (k, c) in ((:z_xy, "ZXY"), (:z_yx, "ZYX")) if mode in (:TETM, k == :z_xy ? :TE : :TM)]
    rows = [(k, c, i, j) for j in 1:ns for i in 1:nf for (k, c) in comps if isfinite(getproperty(data, k)[i, j])]
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# MTGeophysics 2D Fréchet derivative G = ∂d/∂m")
        println(io, "# d: impedance in [mV/km]/[nT], two rows (Re, Im) per data line, exp(+iωt)")
        println(io, "# m: log10 resistivity of the earth cells, column (iz, iy) = iz + (iy - 1) * nz, iz from the model top")
        println(io, "# > rows columns nz ny, then: period site component part g_1 ... g_columns")
        @printf(io, "> %d %d %d %d\n", 2 * length(rows), length(G.cells), nz, ny)
        for (k, c, i, j) in rows
            g = getproperty(G, k)[i + (j - 1) * nf, :] .* scale
            for (part, v) in (("Re", real.(g)), ("Im", imag.(g)))
                print(io, @sprintf("%.8e %s %s %s", data.periods[i], data.site_names[j], c, part))
                foreach(x -> print(io, @sprintf(" %.6e", x)), v)
                println(io)
            end
        end
    end
    String(path)
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
    scale[mt2d_air_mask(mesh)] .= 0
    scale
end

# sampler inputs, per surface row of the station plan: its cell row σ, then Re and Im of its
# two node rows; adjacent rows share nodes, so node indices may repeat
function _mt2d_surface_input(c, s)
    ny = c.t.grid_size[1]
    groups = c.plan.groups
    rows = mapreduce(g -> collect((g.z_index-1)*ny+1:g.z_index*ny), vcat, groups)
    nodes = mapreduce(g -> collect((g.z_index-1)*(ny+1)+1:(g.z_index+1)*(ny+1)), vcat, groups)
    vcat(c.t.conductivity[rows], real.(s.u[nodes]), imag.(s.u[nodes])), rows, nodes
end

function _mt2d_sample(c, s, x)
    groups = c.plan.groups
    ny, ng = c.t.grid_size[1], length(groups)
    n = 2*(ny+1)
    ore, oim = ng*ny, ng*ny + ng*n
    y_nodes = [0.0; cumsum(c.t.y_lengths)]
    surface = map(enumerate(groups)) do (k, g)
        σ = x[(k-1)*ny+1:k*ny]
        pair = reshape(complex.(x[ore+(k-1)*n+1:ore+k*n], x[oim+(k-1)*n+1:oim+k*n]), ny+1, 2)
        _mt2d_surface_fields(s.pol, 2π*s.frequency, y_nodes, c.t.z_lengths[g.z_index], σ, pair)
    end
    _mt2d_station_impedance(s.pol, c.plan, surface, c.locations, y_nodes)
end

function _mt2d_frechet(c, δρ)
    σ = c.t.conductivity
    δσ = -σ.^2 .* vec(permutedims(δρ))
    δσ[vec(permutedims(c.air))] .= 0
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
    x, rows, nodes = _mt2d_surface_input(c, s)
    x̂ = ForwardDiff.gradient(x -> real(dot(ẑ, _mt2d_sample(c, s, x))), x)
    nσ, n = length(rows), length(nodes)
    û = zeros(ComplexF64, length(s.u))
    for (k, i) in enumerate(nodes)
        û[i] += complex(x̂[nσ+k], x̂[nσ+n+k])
    end
    for (k, i) in enumerate(rows)
        σ̂[i] += x̂[k]
    end
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
    ρ̂[c.air] .= 0
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
  indices; requested air columns (air layers and topography) are zero
- for large models use `ApplyFrechet2D` and `ApplyFrechetTranspose2D` instead
"""
function FrechetDerivative2D(mesh::MT2DMesh, ρ::AbstractMatrix{<:Real};
                           mode::Symbol = :TETM, parameterization::Symbol = :resistivity,
                           active_cells = nothing)
    response, c = _mt2d_forward_cache(mesh, ρ; mode)
    scale = _mt2d_parameter_scale(mesh, ρ, parameterization)
    allcells = CartesianIndices(ρ)
    cells = if active_cells === nothing
        [i for i in allcells if !c.air[i]]
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

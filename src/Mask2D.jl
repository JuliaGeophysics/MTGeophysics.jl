# 2D covariance masks
# Author: @pankajkmishra
# The cov.ctrl mask of an earth model, as Mask3D does in 3D: 0 for topographic air and fixed cells, 9 for
# water, 1 for free cells; and the ground surface a model file describes

"""
    Mask2D(model::ModelFile2D; water=nothing, fixed_below_m=Inf, fixed=[]) -> Matrix{Int}

Covariance mask `(nz, ny)` of a ModEM-layout earth model:
- 0 for topographic air (cells above 1e15 ohm m), for cells whose top lies below
  `fixed_below_m` (depth below the model top) and for cells whose centre lies in a
  `fixed` box `(y_range, z_range)` (m, the model's y frame and depth below the top);
- 9 for `water`, a model-shaped Bool mask such as `Topography2D(...).mask .== 9`;
- 1 elsewhere (free).
"""
function Mask2D(model::ModelFile2D; water::Union{Nothing, AbstractMatrix{Bool}} = nothing, fixed_below_m::Real = Inf,
                fixed = NamedTuple[])
    nz, ny = size(model.resistivity)
    mask = ones(Int, nz, ny)
    ztop = vcat(0.0, cumsum(model.z_cell_sizes)[1:end-1])
    zc = ztop .+ model.z_cell_sizes ./ 2
    yc = model.origin[2] .+ cumsum(model.y_cell_sizes) .- model.y_cell_sizes ./ 2
    mask[ztop .>= fixed_below_m, :] .= 0
    for b in fixed
        for iy in 1:ny, iz in 1:nz
            b.y_range[1] <= yc[iy] <= b.y_range[2] && b.z_range[1] <= zc[iz] <= b.z_range[2] && (mask[iz, iy] = 0)
        end
    end
    if water !== nothing
        size(water) == (nz, ny) || throw(DimensionMismatch("water mask must match the model"))
        mask[water] .= MT2D_MASK_WATER
    end
    mask[model.resistivity .> MT2D_AIR_THRESHOLD] .= 0
    mask
end

"""
    mt2d_ground(model::ModelFile2D) -> (y, depth)

Column centres (the model's y frame) and the ground depth below the model top of each
column: the bottom of its topographic air.
"""
function mt2d_ground(model::ModelFile2D)
    topo = _mt2d_model_topo_air(model.resistivity)
    yc = model.origin[2] .+ cumsum(model.y_cell_sizes) .- model.y_cell_sizes ./ 2
    yc, vcat(0.0, cumsum(model.z_cell_sizes))[topo .+ 1]
end

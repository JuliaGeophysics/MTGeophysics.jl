# Project a ModEM model/data pair from one mesh onto another (coarse -> fine or
# fine -> coarse) and emit everything the new mesh needs to run. Names come from
# the inputs, so a run is self-describing on disk:
#
#   NLCGonVFSA/<yyyymmdd_HHMMSS>/
#     modelC_to_modelF.rho   resistivity resampled onto the target grid
#     dataC-W.dat            the same data with every site elevation snapped to
#                            the target grid's surface (-W = rewritten)
#     modelC_to_modelF.cov   covariance with mask 0 on air and 9 on ocean, so
#                            ModEM actually freezes topography and bathymetry
#
# The two meshes need not share extent or origin: everything is resolved in
# absolute coordinates and target cells that fall outside the source are filled
# from the nearest source cell.
#
# Standalone: include this file, it only uses the public MTGeophysics API.
#
#   include("examples/ProjectMesh3D.jl")
#   r = project_model_data("MT3DINV4/best_model.rho", "MT3DINV4/dataC.dat",
#                          "MT3DINV4/modelF.rho"; method=:linear)
#   report(r)

using MTGeophysics
using Printf
using Statistics
using Dates

#----- cell lookup ---------------------------------------------------------#

# index of the cell containing coordinate q, clamped to the grid: target cells
# outside the source extent take the edge value rather than a fill constant
function _cell_index(edges::AbstractVector{<:Real}, q::Real)
    n = length(edges) - 1
    q <= edges[1] && return 1
    q >= edges[end] && return n
    return clamp(searchsortedlast(edges, q), 1, n)
end

_edges(d::AbstractVector{<:Real}, o::Real) = vcat(0.0, cumsum(d)) .+ o

# trilinear weights over the 8 surrounding source cell centers, dropping any
# neighbour that is not usable (air, or masked) and renormalising what is left
function _interp_linear(A, cx, cy, cz, x, y, z, usable)
    i = clamp(searchsortedlast(cx, x), 1, length(cx) - 1)
    j = clamp(searchsortedlast(cy, y), 1, length(cy) - 1)
    k = clamp(searchsortedlast(cz, z), 1, length(cz) - 1)
    tx = clamp((x - cx[i]) / (cx[i+1] - cx[i]), 0.0, 1.0)
    ty = clamp((y - cy[j]) / (cy[j+1] - cy[j]), 0.0, 1.0)
    tz = clamp((z - cz[k]) / (cz[k+1] - cz[k]), 0.0, 1.0)
    num = 0.0; den = 0.0
    for (di, wi) in ((0, 1 - tx), (1, tx)),
        (dj, wj) in ((0, 1 - ty), (1, ty)),
        (dk, wk) in ((0, 1 - tz), (1, tz))
        ii, jj, kk = i + di, j + dj, k + dk
        usable[ii, jj, kk] || continue
        w = wi * wj * wk
        num += w * A[ii, jj, kk]; den += w
    end
    return den > 0 ? num / den : NaN
end

#----- masks on a loaded model ---------------------------------------------#

# air is NaN once loaded; ocean is either an explicit threshold or the exact
# seawater value the mesh was built with
function _masks(m::WS3DModel, water_log10)
    air = isnan.(m.A)
    water = falses(size(m.A))
    if !isnan(water_log10)
        @inbounds for q in eachindex(m.A)
            water[q] = !air[q] && m.A[q] <= water_log10
        end
    end
    return air, water
end

#----- covariance ----------------------------------------------------------#

const _COV_HEADER = """
+-----------------------------------------------------------------------------+
| This file defines model covariance for a recursive autoregression scheme.   |
| The model space may be divided into distinct areas using integer masks.     |
| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |
| air, ocean and the rest of the model is turned off automatically. You can   |
| also define exceptions to override smoothing between any two model areas.   |
| To turn off smoothing set it to zero. This header is 16 lines long.         |
| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
| 2. Smoothing in the X direction (NzEarth real values)                       |
| 3. Smoothing in the Y direction (NzEarth real values)                       |
| 4. Vertical smoothing (1 real value)                                        |
| 5. Number of times the smoothing should be applied (1 integer >= 0)         |
| 6. Number of exceptions (1 integer >= 0)                                    |
| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 2 & 3) |
| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
+-----------------------------------------------------------------------------+"""

"""
    write_covariance(path, air, water; smooth_xy, smooth_z, napply, exceptions)

Write a ModEM covariance file whose mask is 0 on air, 9 on ocean and 1 elsewhere.
Without this file ModEM will happily invert the topographic air cells: the model
file alone does not freeze them.
"""
function write_covariance(path::AbstractString,
                          air::AbstractArray{Bool,3}, water::AbstractArray{Bool,3};
                          smooth_xy::Real = 0.3, smooth_z::Real = 0.3,
                          napply::Integer = 1,
                          exceptions::Vector{NTuple{3,Int}} = NTuple{3,Int}[])
    nx, ny, nz = size(air)
    open(path, "w") do io
        println(io, _COV_HEADER)
        println(io)
        @printf(io, " %d %d %d\n", nx, ny, nz)
        println(io)
        for _ in 1:2
            println(io, " " * join((@sprintf("%.4g", smooth_xy) for _ in 1:nz), "  "))
            println(io)
        end
        @printf(io, " %.4g\n\n", smooth_z)
        @printf(io, " %d\n\n", napply)
        @printf(io, " %d\n", length(exceptions))
        for e in exceptions
            @printf(io, " %d %d %d.\n", e[1], e[2], e[3])
        end
        println(io)
        for k in 1:nz
            @printf(io, " %d %d\n", k, k)
            # rows run north to south: ModEM's read_iscalar takes the first row of
            # a block as x = Nx (sg_scalar.f90, `do j = Nx,1,-1`). Writing them
            # ascending mirrors the mask about x and frees half the topography.
            for i in nx:-1:1
                println(io, " " * join((air[i, j, k] ? "0" :
                                        water[i, j, k] ? "9" : "1" for j in 1:ny), " "))
            end
            println(io)
        end
    end
    return path
end

#----- the projection ------------------------------------------------------#

"""
    project_model_data(src_model, src_data, dst_mesh; kwargs...)

Resample `src_model` onto the grid of `dst_mesh`, snap the site elevations in
`src_data` to the target surface, and write a matching covariance file.

Works in either direction (coarse to fine or fine to coarse) and does not
require the two meshes to share extent, origin or rotation of the grid lines.
Target cells outside the source extent take the nearest source value.

The target mesh keeps **its own** topography and bathymetry: cells that are air
or ocean in `dst_mesh` are left at their own values and never filled from the
source, so the projected model honours the new mesh's surface rather than
inheriting a coarser one.

By default the three files land in a fresh timestamped directory beside the
target mesh and are named after their inputs, so nothing is ever overwritten and
the provenance of an NLCG start is readable from the filenames alone:

    <dst dir>/NLCGonVFSA/20260726_150322/modelC_to_modelF.rho
                                        /dataC-W.dat
                                        /modelC_to_modelF.cov

Keywords
- `outdir`     override the timestamped directory
- `stem`       override the `<src>_to_<dst>` name shared by the model and cov
- `method`     `:nearest` (default, preserves contrasts) or `:linear`
- `water_log10` log10 threshold below which a cell counts as ocean, `NaN` for none
- `snap_sites` set false to leave the data file's elevations untouched
- `smooth_xy`, `smooth_z`, `napply` covariance smoothing parameters

Returns a NamedTuple of the three output paths plus diagnostics.
"""
function project_model_data(src_model::AbstractString, src_data::AbstractString,
                            dst_mesh::AbstractString;
                            outdir::AbstractString = joinpath(dirname(dst_mesh),
                                "NLCGonVFSA", Dates.format(now(), "yyyymmdd_HHMMSS")),
                            stem::AbstractString = string(
                                first(splitext(basename(src_model))), "_to_",
                                first(splitext(basename(dst_mesh)))),
                            method::Symbol = :nearest,
                            water_log10::Real = NaN,
                            snap_sites::Bool = true,
                            smooth_xy::Real = 0.3, smooth_z::Real = 0.3,
                            napply::Integer = 1)

    method in (:nearest, :linear) || error("method must be :nearest or :linear")

    src = load_ws3d_model(src_model)
    dst = load_ws3d_model(dst_mesh)
    _, _, _, _, _, _, _, rotation = read_ws3d_model(dst_mesh, true)

    src_air, src_water = _masks(src, water_log10)
    dst_air, dst_water = _masks(dst, water_log10)
    usable = .!src_air                       # never pull a value out of the air

    sx = _edges(src.dx, src.origin[1])
    sy = _edges(src.dy, src.origin[2])
    sz = _edges(src.dz, src.origin[3])

    A = similar(dst.A)
    n_outside = 0
    n_kept = 0
    @inbounds for k in 1:dst.nz, j in 1:dst.ny, i in 1:dst.nx
        # the target mesh's own surface wins
        if dst_air[i, j, k] || dst_water[i, j, k]
            A[i, j, k] = dst.A[i, j, k]; n_kept += 1; continue
        end
        x, y, z = dst.cx[i], dst.cy[j], dst.cz[k]
        outside = x < sx[1] || x > sx[end] || y < sy[1] || y > sy[end] ||
                  z < sz[1] || z > sz[end]
        outside && (n_outside += 1)
        v = if method === :linear && !outside
            _interp_linear(src.A, src.cx, src.cy, src.cz, x, y, z, usable)
        else
            NaN
        end
        if !isfinite(v)
            # nearest cell, walking up the column if we land in source air
            ii = _cell_index(sx, x); jj = _cell_index(sy, y); kk = _cell_index(sz, z)
            while kk >= 1 && src_air[ii, jj, kk]; kk -= 1; end
            kk < 1 && (kk = findfirst(!, view(src_air, ii, jj, :)))
            v = kk === nothing ? median(filter(isfinite, src.A)) : src.A[ii, jj, kk]
        end
        A[i, j, k] = v
    end

    mkpath(outdir)
    # -W marks the data file as rewritten: same observations, new site elevations
    dstem, dext = splitext(basename(src_data))
    model_out = joinpath(outdir, stem * ".rho")
    data_out  = joinpath(outdir, dstem * "-W" * (isempty(dext) ? ".dat" : dext))
    cov_out   = joinpath(outdir, stem * ".cov")

    write_ws3d_model(model_out, dst.dx, dst.dy, dst.dz, A, dst.origin;
                     rotation=rotation, type_str="LOGE")
    write_covariance(cov_out, dst_air, dst_water;
                     smooth_xy=smooth_xy, smooth_z=smooth_z, napply=napply)

    #---------- site elevations onto the target surface ----------
    nsite, n_off_grid, dz_moved = _rewrite_elevations(src_data, data_out, dst, snap_sites)

    return (model = model_out, data = data_out, cov = cov_out,
            src_grid = (src.nx, src.ny, src.nz),
            dst_grid = (dst.nx, dst.ny, dst.nz),
            cells_outside_source = n_outside,
            surface_cells_kept = n_kept,
            n_air = count(dst_air), n_water = count(dst_water),
            sites = nsite, sites_off_grid = n_off_grid,
            elev_shift = isempty(dz_moved) ? (0.0, 0.0, 0.0) :
                (minimum(dz_moved), median(dz_moved), maximum(dz_moved)))
end

# Rewrites only the Z(m) column, streaming the file line by line. Going through
# load/write_data_modem would be lossy here: that writer demands one unit across
# every block and so rejects the usual pairing of [mV/km]/[nT] impedance with
# dimensionless tipper. Editing in place also keeps periods, error floors, signs
# and any extra blocks exactly as they were.
function _rewrite_elevations(src_data::AbstractString, out_path::AbstractString,
                             dst::WS3DModel, snap::Bool)
    ex = _edges(dst.dx, dst.origin[1])
    ey = _edges(dst.dy, dst.origin[2])
    shifts = Float64[]
    seen = Set{String}(); off = Set{String}()
    open(out_path, "w") do io
        for line in eachline(src_data)
            s = lstrip(line)
            p = split(s)
            if !snap || isempty(s) || startswith(s, "#") || startswith(s, ">") || length(p) < 8
                println(io, line); continue
            end
            x = tryparse(Float64, p[5]); y = tryparse(Float64, p[6])
            zold = tryparse(Float64, p[7])
            if x === nothing || y === nothing || zold === nothing
                println(io, line); continue
            end
            code = String(p[2])
            (x < ex[1] || x > ex[end] || y < ey[1] || y > ey[end]) && push!(off, code)
            i = _cell_index(ex, x); j = _cell_index(ey, y)
            znew = dst.Z[i, j]              # bottom edge of the last air/ocean cell
            if !(code in seen)
                push!(seen, code); push!(shifts, znew - zold)
            end
            p[7] = @sprintf("%.3f", znew)
            println(io, join(p, " "))
        end
    end
    return length(seen), length(off), shifts
end

"""
    report(r)

Print the NamedTuple returned by [`project_model_data`](@ref).
"""
function report(r::NamedTuple)
    @printf("  source grid %s  ->  target grid %s\n", r.src_grid, r.dst_grid)
    @printf("  target cells outside the source extent : %d (filled from nearest)\n",
            r.cells_outside_source)
    @printf("  target surface cells kept as-is        : %d (%d air, %d ocean)\n",
            r.surface_cells_kept, r.n_air, r.n_water)
    @printf("  sites %d (%d outside the target grid)\n", r.sites, r.sites_off_grid)
    @printf("  site elevation shift min/median/max    : %.1f / %.1f / %.1f m\n",
            r.elev_shift...)
    println("  wrote ", r.model)
    println("        ", r.data)
    println("        ", r.cov)
    return nothing
end

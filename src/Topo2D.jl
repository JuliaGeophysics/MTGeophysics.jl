# 2D MT topography
# Author: @pankajkmishra
# topo.dat (WGS84 lat, lon, elevation) along a profile, and ModEM-layout models with the
# topographic air and water cut in, following MakeMesh3D in 3D

using Printf

"""
    Topo2D

Ground elevation points along a profile: WGS84 `latitudes`, `longitudes` and
`elevations` in metres above sea level.
"""
Base.@kwdef struct Topo2D
    latitudes::Vector{Float64}
    longitudes::Vector{Float64}
    elevations::Vector{Float64}
end

"""
    ReadTopo2D(path) -> Topo2D

Read `topo.dat`: one `lat lon elevation` point per line, WGS84 degrees and metres
above sea level; `#` lines are skipped.
"""
function ReadTopo2D(path::AbstractString)
    isfile(path) || error("topography file not found: $path")
    lat, lon, elev = Float64[], Float64[], Float64[]
    for line in eachline(path)
        s = strip(first(split(line, '#')))
        isempty(s) && continue
        p = split(s)
        length(p) >= 3 || error("$path: expected 'lat lon elevation', got '$s'")
        push!(lat, parse(Float64, p[1])); push!(lon, parse(Float64, p[2])); push!(elev, parse(Float64, p[3]))
    end
    length(lat) >= 2 || error("$path: at least two topography points are needed")
    Topo2D(latitudes = lat, longitudes = lon, elevations = elev)
end

"""
    WriteTopo2D(path, topo::Topo2D) -> path

Write `topo.dat`: a `#` header, then `lat lon elevation` per point.
"""
function WriteTopo2D(path::AbstractString, t::Topo2D)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# 2D topography written by MTGeophysics.jl, WGS84")
        println(io, "# Lat(deg) Lon(deg) Elevation(m a.s.l.)")
        for i in eachindex(t.latitudes)
            @printf(io, "%11.6f %11.6f %10.2f\n", t.latitudes[i], t.longitudes[i], t.elevations[i])
        end
    end
    String(path)
end

#---------- projection onto the profile ----------

# local metres (east, north) around the survey origin
function _mt2d_local_xy(data, lat, lon)
    trans = Proj.Transformation("EPSG:4326", _local_tm_proj_string(data.origin[1], data.origin[2]); always_xy = true)
    p = [trans(lon[i], lat[i]) for i in eachindex(lat)]
    [q[1] for q in p], [q[2] for q in p]
end

"""
    mt2d_profile_topography(topo, data) -> (y, elevation)

Profile position y (the data's local y) of each topography point, sorted along y. Each
point goes to the nearest point of the station polyline, so curved lines stay curved;
the end segments extend beyond the first and last station.
"""
function mt2d_profile_topography(topo::Topo2D, data)
    ns = length(data.receivers)
    ns >= 2 || throw(ArgumentError("projecting topography needs at least two stations"))
    length(data.latitudes) == ns && !all(iszero, data.latitudes) ||
        throw(ArgumentError("the data file has no station lat/lon to place the topography"))
    order = sortperm(data.receivers)
    sy = data.receivers[order]
    se, sn = _mt2d_local_xy(data, data.latitudes[order], data.longitudes[order])
    te, tn = _mt2d_local_xy(data, topo.latitudes, topo.longitudes)
    y = map(eachindex(te)) do k
        best, ybest = Inf, NaN
        for i in 1:ns-1
            de, dn = se[i+1] - se[i], sn[i+1] - sn[i]
            len2 = de^2 + dn^2
            len2 > 0 || continue
            t = ((te[k] - se[i]) * de + (tn[k] - sn[i]) * dn) / len2
            t = clamp(t, i == 1 ? -Inf : 0.0, i == ns - 1 ? Inf : 1.0)
            d = hypot(te[k] - se[i] - t * de, tn[k] - sn[i] - t * dn)
            d < best && ((best, ybest) = (d, sy[i] + t * (sy[i+1] - sy[i])))
        end
        ybest
    end
    p = sortperm(y)
    y[p], topo.elevations[p]
end

# linear interpolation, constant beyond the ends
function _mt2d_interp(x::AbstractVector, v::AbstractVector, xi::Real)
    xi <= x[1] && return v[1]
    xi >= x[end] && return v[end]
    j = searchsortedlast(x, xi)
    x[j+1] == x[j] && return v[j]
    v[j] + (v[j+1] - v[j]) * (xi - x[j]) / (x[j+1] - x[j])
end

#---------- model with topography ----------

"""
    Topography2D(model::ModelFile2D, data::DataFile2D, topo::Topo2D;
                 water=[], water_resistivity=100.0) -> (; model, mask, data, datum, ground)

Cut topography into a ModEM-layout earth model, as `MakeMesh3D` does in 3D.
- The model top becomes the datum: the highest station elevation, read from `topo` at
  the stations. Ground higher than the datum is flattened to the model top.
- Cells whose centre lies above the ground are air, written as 1e17 ohm m with mask 0.
- In each station column the ground moves, up or down, to the cell boundary nearest
  the station (the shallowest station's when several share a column), so every station
  sits on its column's ground within half a cell.
- `water`: lakes or sea as `(y_range = (y1, y2), level = m a.s.l.)`. Cells between the
  level and the ground are water, `water_resistivity`, mask 9. No station may stand
  over water.

Returns the new model, the covariance mask (0 air, 9 water, 1 free), the data with
Z = datum − station elevation, the datum (m a.s.l.) and the ground depth of each column.
"""
function Topography2D(model::ModelFile2D, data, topo::Topo2D;
                      water = NamedTuple[], water_resistivity::Real = 100.0)
    model.n_air_cells == 0 || throw(ArgumentError("model must hold earth cells only, read it with ReadModel2D"))
    ty, te = mt2d_profile_topography(topo, data)
    elevation(y) = _mt2d_interp(ty, te, y)
    station_elevation = elevation.(data.receivers)
    datum = maximum(station_elevation)

    nz, ny = size(model.resistivity)
    y_nodes = model.origin[2] .+ vcat(0.0, cumsum(model.y_cell_sizes))
    yc = (y_nodes[1:end-1] .+ y_nodes[2:end]) ./ 2
    zn = vcat(0.0, cumsum(model.z_cell_sizes))
    zc = (zn[1:end-1] .+ zn[2:end]) ./ 2

    # earth cells start below the ground; a station column takes the boundary nearest its station, up or
    # down, so a station in a dip narrower than a column is not left under the ground (the shallowest wins
    # when stations share a column)
    ground = [count(<(datum - elevation(y)), zc) for y in yc]
    columns = [searchsortedlast(y_nodes, y) for y in data.receivers]
    snapped = Dict{Int, Int}()
    for (i, iy) in enumerate(columns)
        1 <= iy <= ny || throw(ArgumentError("station $(data.site_names[i]) lies outside the model"))
        k = argmin(abs.(zn .- (datum - station_elevation[i]))) - 1
        snapped[iy] = min(get(snapped, iy, k), k)
    end
    foreach(((iy, k),) -> ground[iy] = k, snapped)
    ground = min.(ground, nz - 1)

    # water between its level and the ground; air above
    top = copy(ground)
    for w in water
        level = datum - Float64(w.level)
        for iy in eachindex(yc)
            w.y_range[1] <= yc[iy] <= w.y_range[2] || continue
            iy in columns && throw(ArgumentError("water at y = $(yc[iy]) m lies under a station; place stations on land"))
            top[iy] = min(top[iy], count(<(level), zc))
        end
    end

    ρ = copy(model.resistivity)
    mask = ones(Int, nz, ny)
    for iy in 1:ny
        ρ[1:top[iy], iy] .= MT2D_AIR_TAG
        mask[1:top[iy], iy] .= 0
        ρ[top[iy]+1:ground[iy], iy] .= water_resistivity
        mask[top[iy]+1:ground[iy], iy] .= MT2D_MASK_WATER
    end

    out = deepcopy(data)
    out.z_positions .= datum .- station_elevation
    newmodel = ModelFile2D(title = model.title, x_cell_sizes = model.x_cell_sizes, y_cell_sizes = model.y_cell_sizes,
                           z_cell_sizes = model.z_cell_sizes, resistivity = ρ, n_air_cells = 0, origin = model.origin,
                           rotation = model.rotation, format = model.format, path = "")
    (; model = newmodel, mask, data = out, datum, ground = zn[ground .+ 1])
end

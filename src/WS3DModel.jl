# WS3D format 3D resistivity model I/O (log10 internal representation).
# Author: @pankajkmishra
# Unlike ModEMModel (linear Ω·m), WS3DModel stores resistivity in log10(Ω·m).
# This is the native format used by the 3D VFSA inversion engine.

mutable struct WS3DModel
    A::Array{Float64,3}        # log10 resistivity
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    nx::Int
    ny::Int
    nz::Int
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    origin::Vector{Float64}
    cx::Vector{Float64}
    cy::Vector{Float64}
    cz::Vector{Float64}
    X::Matrix{Float64}
    Y::Matrix{Float64}
    Xc::Matrix{Float64}
    Yc::Matrix{Float64}
    Z::Matrix{Float64}
    npad::NTuple{2,Int}
    name::String
    niter::String
end

"""
    read_ws3d_model(fname; block=true)

Low-level parser for WS3D model files returning a tuple of raw arrays.
Handles both LOGE and LINEAR formats, converting to log10 internally.
"""
function read_ws3d_model(fname::AbstractString, block::Bool=true)
    lines = readlines(fname)
    i = 1
    while i ≤ length(lines) && (isempty(lines[i]) || startswith(strip(lines[i]), "#"))
        i += 1
    end
    hdr = strip(lines[i]); i += 1
    toks = split(hdr)
    ints = Int[]
    for t in toks
        v = tryparse(Int, t)
        v === nothing || push!(ints, v)
    end
    nx, ny, nz, nzAir = ints[1], ints[2], ints[3], ints[4]
    type = occursin("LOGE", hdr) ? "LOGE" : "LINEAR"

    function take_floats!(n::Int)
        vals = Float64[]
        while length(vals) < n && i ≤ length(lines)
            ts = split(strip(lines[i]))
            for s in ts
                v = tryparse(Float64, s)
                if v !== nothing
                    push!(vals, v)
                    length(vals) == n && break
                end
            end
            i += 1
        end
        return vals
    end

    dx = take_floats!(nx)
    dy = take_floats!(ny)
    dz = take_floats!(nz)

    A = zeros(nx, ny, nz)
    if block
        for k in 1:nz
            vals = take_floats!(nx*ny)
            tmp = reshape(vals, nx, ny)
            A[:,:,k] = reverse(tmp, dims=1)
        end
    else
        kmax = 0
        while kmax < nz && i ≤ length(lines)
            s = strip(lines[i]); i += 1
            isempty(s) && continue
            ks = split(s)
            kk = Int[]
            for t in ks
                v = tryparse(Int, t)
                v === nothing || push!(kk, v)
            end
            isempty(kk) && continue
            length(kk) == 1 && (kk = [kk[1], kk[1]])
            vals = take_floats!(nx*ny)
            tmp = reshape(vals, nx, ny)
            for k in kk[1]:kk[2]
                A[:,:,k] = tmp
            end
            kmax = max(kmax, kk[2])
        end
    end

    origin = [0.0, 0.0, 0.0]
    rotation = 0.0
    while i ≤ length(lines)
        s = strip(lines[i]); i += 1
        isempty(s) && continue
        ts = split(s)
        nums = Float64[]
        for t in ts
            v = tryparse(Float64, t)
            v === nothing || push!(nums, v)
        end
        if length(nums) == 3
            origin = nums
        elseif length(nums) == 1
            rotation = nums[1]
        end
    end

    if type == "LOGE"
        A .= A ./ log(10.0)   # natural log → log10
    else
        @. A = log10(A)        # linear → log10
    end

    return dx, dy, dz, A, nzAir, type, origin, rotation
end

"""
    load_ws3d_model(name::AbstractString) -> WS3DModel

Load a WS3D model file into a `WS3DModel` struct with computed mesh coordinates,
cell centers, topography surface, and padding estimates.
"""
function load_ws3d_model(name::AbstractString)
    dx, dy, dz, A, nzAir, type, origin, rotation = read_ws3d_model(name, true)
    m = WS3DModel(
        A, dx, dy, dz, size(A,1), size(A,2), size(A,3),
        Float64[], Float64[], Float64[], collect(origin),
        Float64[], Float64[], Float64[],
        zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0),
        (0,0), String(name), ""
    )

    m.y = vcat(0.0, cumsum(m.dy)) .+ m.origin[2]
    m.x = vcat(0.0, cumsum(m.dx)) .+ m.origin[1]
    m.z = vcat(0.0, cumsum(m.dz)) .+ m.origin[3]

    m.A[m.A .> 15.0] .= NaN   # anything > 10^15 Ω·m is air

    m.cx = (m.x[1:end-1] .+ m.x[2:end]) ./ 2
    m.cy = (m.y[1:end-1] .+ m.y[2:end]) ./ 2
    m.cz = (m.z[1:end-1] .+ m.z[2:end]) ./ 2
    m.X, m.Y = meshgrid(m.y, m.x)
    m.Xc, m.Yc = meshgrid(m.cy, m.cx)
    m.nx = length(m.cx)
    m.ny = length(m.cy)
    m.nz = length(m.cz)

    Z = zeros(m.nx, m.ny)
    @inbounds for i in 1:m.nx, j in 1:m.ny
        col_log10 = view(m.A, i, j, :)
        col_lin = @. 10.0^(col_log10)
        indair = findlast(isnan, col_log10)
        indocean = findlast(x -> abs(x - 0.3) < 1e-5, col_lin)
        ind = maximum((indair === nothing ? 0 : indair,
                       indocean === nothing ? 0 : indocean))
        ind == m.nz && (ind = m.nz - 1)
        Z[i,j] = m.z[ind+1]
    end
    m.Z = Z

    cx0 = m.dx[clamp(round(Int, m.nx/2), 1, length(m.dx))]
    cy0 = m.dy[clamp(round(Int, m.ny/2), 1, length(m.dy))]
    nx_core = count(==(cx0), m.dx)
    ny_core = count(==(cy0), m.dy)
    m.npad = (Int((m.nx - nx_core) ÷ 2), Int((m.ny - ny_core) ÷ 2))

    return m
end

"""
    write_ws3d_model(fname, dx, dy, dz, A_log10, origin; rotation=0.0, type_str="LOGE")

Write a WS3D model file from log10 resistivity array.
"""
function write_ws3d_model(fname::AbstractString,
                          dx::Vector{Float64}, dy::Vector{Float64}, dz::Vector{Float64},
                          A_log10::Array{Float64,3},
                          origin::AbstractVector{<:Real}=[0.0,0.0,0.0];
                          rotation::Real=0.0,
                          type_str::String="LOGE")
    nx, ny, nz = size(A_log10)
    nzAir = 0
    ln10 = log(10.0)
    ln1e17 = 17.0 * ln10

    vals = if type_str == "LOGE"
        map(x -> isfinite(x) ? (ln10 * x) : ln1e17, A_log10)
    else
        map(x -> isfinite(x) ? 10.0^x : 1e17, A_log10)
    end

    open(fname, "w") do io
        println(io, "# Written by MTGeophysics.jl write_ws3d_model")
        @printf(io, "%d %d %d %d %s\n", nx, ny, nz, nzAir, type_str)
        for v in dx; @printf(io, "%G ", v); end; println(io)
        for v in dy; @printf(io, "%G ", v); end; println(io)
        for v in dz; @printf(io, "%G ", v); end; println(io)
        for k in 1:nz
            println(io)
            for j in 1:ny
                for i in nx:-1:1
                    @printf(io, "%15.5E", vals[i,j,k])
                end
                println(io)
            end
        end
        @printf(io, "%g %g %g\n", origin[1], origin[2], origin[3])
        @printf(io, "%g\n", rotation)
    end
end



mutable struct ModEMModel
    A::Array{Float64,3}
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

meshgrid(ax::AbstractVector, ay::AbstractVector) = (repeat(reshape(ax, 1, :), length(ay), 1), repeat(reshape(ay, :, 1), 1, length(ax)))

function read_mackie3d_model(fname::AbstractString, block::Bool=true)
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
        if v !== nothing
            push!(ints, v)
        end
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
                    if length(vals) == n
                        break
                    end
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
            if isempty(s)
                continue
            end
            ks = split(s)
            kk = Int[]
            for t in ks
                v = tryparse(Int, t)
                if v !== nothing
                    push!(kk, v)
                end
            end
            if isempty(kk)
                continue
            end
            if length(kk) == 1
                kk = [kk[1], kk[1]]
            end
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
        if isempty(s)
            continue
        end
        ts = split(s)
        nums = Float64[]
        for t in ts
            v = tryparse(Float64, t)
            if v !== nothing
                push!(nums, v)
            end
        end
        if length(nums) == 3
            origin = nums
        elseif length(nums) == 1
            rotation = nums[1]
        end
    end
    return dx, dy, dz, A, nzAir, type, origin, rotation
end

function load_model_modem(name::AbstractString)
    # This part is inspired by what was originally written in MATLAB repo (10.5281/zenodo.6784201)
    dx, dy, dz, A, nzAir, type, origin, rotation = read_mackie3d_model(name, true)
    m = ModEMModel(A, dx, dy, dz, size(A,1), size(A,2), size(A,3), Float64[], Float64[], Float64[], collect(origin), Float64[], Float64[], Float64[], zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0), (0,0), String(name), "")
    if type == "LOGE"
        m.A = exp.(m.A)
    end
    m.y = vcat(0.0, cumsum(m.dy)) .+ m.origin[2]
    m.x = vcat(0.0, cumsum(m.dx)) .+ m.origin[1]
    m.z = vcat(0.0, cumsum(m.dz)) .+ m.origin[3]
    m.A[m.A .> 1e15] .= NaN
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
        col = view(m.A, i, j, :)
        indair = findlast(isnan, col)
        indocean = findlast(x -> abs(x - 0.3) < 1e-5, col)
        ind = maximum((indair === nothing ? 0 : indair, indocean === nothing ? 0 : indocean))
        if ind == 0
            ind = 0
        end
        if ind == m.nz
            ind = m.nz - 1
        end
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



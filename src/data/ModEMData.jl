# LinearAlgebra and Statistics imported from main module

mutable struct ModEMData
    T::Vector{Float64}
    f::Vector{Float64}
    zrot::Matrix{Float64}
    trot::Matrix{Float64}
    site::Vector{String}
    loc::Matrix{Float64}
    ns::Int
    nf::Int
    nr::Int
    responses::Vector{String}
    Z::Array{ComplexF64,3}
    Zerr::Array{ComplexF64,3}
    rho::Array{Float64,3}
    rhoerr::Array{Float64,3}
    pha::Array{Float64,3}
    phaerr::Array{Float64,3}
    tip::Array{ComplexF64,3}
    tiperr::Array{ComplexF64,3}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    origin::Vector{Float64}
    niter::String
    name::String
end

function make_nan_data()
    return ModEMData(
        Float64[], Float64[], Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0),
        String[], Matrix{Float64}(undef, 0, 0), 0, 0, 0, String[],
        Array{ComplexF64}(undef, 0, 0, 0), Array{ComplexF64}(undef, 0, 0, 0),
        Array{Float64}(undef, 0, 0, 0), Array{Float64}(undef, 0, 0, 0),
        Array{Float64}(undef, 0, 0, 0), Array{Float64}(undef, 0, 0, 0),
        Array{ComplexF64}(undef, 0, 0, 0), Array{ComplexF64}(undef, 0, 0, 0),
        Float64[], Float64[], Float64[], Float64[], "", ""
    )
end

# -------------------- small helpers --------------------
strip_gt(s::AbstractString) = begin
    t = strip(s)
    startswith(t, ">") ? strip(t[2:end]) : t
end

safeparsefloat(s) = something(tryparse(Float64, s), 0.0)
safeparseint(s)   = something(tryparse(Int, s), 0)

is_dataline(line::AbstractString) = begin
    t = strip(line)
    !isempty(t) &&
    !startswith(t, "#") &&
    !startswith(t, ">") &&
    occursin(r"^[\s\+\-\.0-9]", t)
end

# -------------------- main loader ----------------------
function load_data_modem(path::AbstractString)
    println("Loading ModEM Data File: $path")

    μ0 = 4π * 1e-7

    d = make_nan_data()
    d.niter = ""
    d.name  = path

    # collectors
    T          = Float64[]
    site       = String[]
    loc_lat    = Float64[]
    loc_lon    = Float64[]
    loc_elev   = Float64[]
    xvals      = Float64[]
    yvals      = Float64[]
    zvals      = Float64[]
    responses  = String[]
    D          = ComplexF64[]
    Derr       = ComplexF64[]

    header_units = String[]
    isign        = Int[]
    rotation     = Float64[]
    origins      = Vector{NTuple{2,Float64}}()

    hb = 0

    open(path, "r") do io
        while !eof(io)
            line = readline(io)
            t = strip(line)
            if isempty(t)
                continue
            end
            if startswith(t, "#")
                continue
            elseif startswith(t, ">")
                hb += 1
                hdr = String[t]
                for _ in 1:5
                    if !eof(io)
                        push!(hdr, readline(io))
                    end
                end
                dataType = strip_gt(hdr[1])
                sstr     = strip_gt(hdr[2])
                push!(isign, occursin("-", lowercase(sstr)) ? -1 : 1)
                ustr     = strip_gt(hdr[3])
                push!(header_units, ustr)
                rtok = split(strip_gt(hdr[4]))
                push!(rotation, length(rtok)>=1 ? safeparsefloat(rtok[1]) : 0.0)
                otok = split(strip_gt(hdr[5]))
                olat = length(otok)>=1 ? safeparsefloat(otok[1]) : 0.0
                olon = length(otok)>=2 ? safeparsefloat(otok[2]) : 0.0
                push!(origins, (olat, olon))
                # nf/ns line exists at hdr[6], but we don't need it here
            elseif is_dataline(t)
                parts = split(t)
                if length(parts) >= 11
                    push!(T,     parse(Float64, parts[1]))
                    push!(site,  parts[2])
                    push!(loc_lat,  parse(Float64, parts[3]))
                    push!(loc_lon,  parse(Float64, parts[4]))
                    push!(loc_elev, parse(Float64, parts[7]))
                    push!(xvals, parse(Float64, parts[5]))
                    push!(yvals, parse(Float64, parts[6]))
                    push!(zvals, parse(Float64, parts[7]))
                    push!(responses, parts[8])
                    re = parse(Float64, parts[9]); im = parse(Float64, parts[10])
                    σ  = parse(Float64, parts[11])
                    push!(D,    ComplexF64(re, im))
                    push!(Derr, ComplexF64(σ, σ))
                end
            end
        end
    end

    # unique sites (stable)
    uniq_sites = unique(site)
    inds = [findfirst(==(s), site) for s in uniq_sites]
    d.site = uniq_sites
    d.ns   = length(uniq_sites)

    d.loc = hcat(loc_lat, loc_lon, loc_elev)[inds, :]
    d.x   = xvals[inds]
    d.y   = yvals[inds]
    d.z   = zvals[inds]

    d.T  = sort(unique(T))
    d.f  = 1.0 ./ d.T
    d.nf = length(d.f)

    d.responses = unique(responses)
    d.nr        = length(d.responses)

    # allocate
    d.Z      = fill(ComplexF64(NaN, NaN), d.nf, 4, d.ns)
    d.Zerr   = fill(ComplexF64(NaN, NaN), d.nf, 4, d.ns)
    d.rho    = fill(NaN, d.nf, 4, d.ns)
    d.rhoerr = fill(NaN, d.nf, 4, d.ns)
    d.pha    = fill(NaN, d.nf, 4, d.ns)
    d.phaerr = fill(NaN, d.nf, 4, d.ns)
    d.tip    = fill(ComplexF64(NaN, NaN), d.nf, 2, d.ns)
    d.tiperr = fill(ComplexF64(NaN, NaN), d.nf, 2, d.ns)

    compidx = Dict("ZXX"=>1, "ZXY"=>2, "ZYX"=>3, "ZYY"=>4, "TX"=>5, "TY"=>6)

    for i in eachindex(T)
        ifreq = findfirst(==(T[i]), d.T)
        isite = findfirst(==(site[i]), d.site)
        icomp = get(compidx, responses[i], 0)
        if !(isnothing(ifreq) || isnothing(isite) || icomp==0)
            if icomp >= 5
                d.tip[ifreq, icomp-4, isite]    = D[i]
                d.tiperr[ifreq, icomp-4, isite] = Derr[i]
            else
                d.Z[ifreq, icomp, isite]    = D[i]
                d.Zerr[ifreq, icomp, isite] = Derr[i]
            end
        end
    end

    # Units: convert mV/km/nT to SI if present (case/space-insensitive)
    for u in header_units
        lu = lowercase(replace(strip(u), " " => ""))
        if lu == "[mv/km]/[nt]"
            d.Z    .*= (μ0 * 1000)
            d.Zerr .*= (μ0 * 1000)
            break
        elseif lu == "[v/m]/[t]" || lu == "[]"
            break
        end
    end

    # Consistency checks
    if !isempty(isign) && !all(==(isign[1]), isign)
        error("Sign convention is inconsistent between data blocks.")
    end
    if !isempty(rotation) && !all(==(rotation[1]), rotation)
        error("Rotation is inconsistent between data blocks.")
    end
    if !isempty(origins)
        olats = [o[1] for o in origins]
        olons = [o[2] for o in origins]
        if !(all(==(olats[1]), olats) && all(==(olons[1]), olons))
            error("Mesh origin is inconsistent between data blocks.")
        end
    end

    # store rotations/origin
    if !isempty(rotation)
        if rotation[1] != 0
            @warn "Caution: non-zero rotation in ModEM data file may cause issues." rotation=rotation[1]
        end
        d.zrot = fill(rotation[1], d.nf, d.ns)
        d.trot = d.zrot
    else
        d.zrot = zeros(d.nf, d.ns)
        d.trot = d.zrot
    end

    if !isempty(origins)
        d.origin = [origins[1][1], origins[1][2], 0.0]
    else
        d.origin = [0.0, 0.0, 0.0]
    end

    # sign convention
    if !isempty(isign)
        if isign[1] == -1
            d.Z   .= conj.(d.Z)
            d.tip .= conj.(d.tip)
        elseif isign[1] != 1
            error("Sign convention must be e^{+iωt} or e^{-iωt}")
        end
    end

    # Huge errors → NaN (check only finite magnitudes)
    magsZ = abs.(d.Zerr)
    finiteZ = isfinite.(magsZ)
    if any(finiteZ) && all(magsZ[finiteZ] .> 1e10)
        d.Zerr .= ComplexF64(NaN, NaN)
    end
    magsT = abs.(d.tiperr)
    finiteT = isfinite.(magsT)
    if any(finiteT) && all(magsT[finiteT] .> 1e10)
        d.tiperr .= ComplexF64(NaN, NaN)
    end

    d.rho, d.pha, d.rhoerr, d.phaerr = calc_rho_pha(d.Z, d.Zerr, d.T)
    return d
end

function calc_rho_pha(Z::Array{ComplexF64,3}, Zerr::Array{ComplexF64,3}, T::Vector{Float64})
    nf, ncomp, ns = size(Z)
    μ0 = 4π * 1e-7

    rho    = fill(NaN, nf, ncomp, ns)
    pha    = fill(NaN, nf, ncomp, ns)
    rhoerr = fill(NaN, nf, ncomp, ns)
    phaerr = fill(NaN, nf, ncomp, ns)

    for i in 1:nf
        ω = 2π / T[i]
        for j in 1:ncomp
            for k in 1:ns
                zval = Z[i,j,k]
                if !isnan(real(zval)) && !isnan(imag(zval))
                    ρa  = abs(zval)^2 / (ω * μ0)
                    phi = rad2deg(angle(zval))
                    rho[i,j,k] = ρa
                    pha[i,j,k] = phi

                    zerr = Zerr[i,j,k]
                    if !isnan(real(zerr)) && !isnan(imag(zerr))
                        Z_abs    = abs(zval)
                        Zerr_abs = abs(zerr)
                        if Z_abs > 0
                            rhoerr[i,j,k] = 2 * ρa * (Zerr_abs / Z_abs)
                            phaerr[i,j,k] = rad2deg(Zerr_abs / Z_abs)
                        end
                    end
                end
            end
        end
    end
    return rho, pha, rhoerr, phaerr
end

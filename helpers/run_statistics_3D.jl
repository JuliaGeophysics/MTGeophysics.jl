#!/usr/bin/env julia
# ensemble statistics over VFSA chains
#
# walks a project directory, collects best_model.rho from every chain, checks the
# grids agree, and writes per-cell mean/median/std/p16/p84 plus a text summary.
# stats are computed in LOGE (the file's native unit), so the "mean" is a
# geometric mean of resistivity -- the right average for a log-parameterized
# inversion.
#
#   julia ensemble_stats.jl                  # PROJECT_DIR below
#   julia ensemble_stats.jl /path/to/project
#   julia ensemble_stats.jl /path/to/project --min-iter=500 --rms-tol=1.15
#
# chains still running are fine: a half-written best_model.rho is skipped with a
# warning rather than aborting the run.

using Printf

#---------- config ----------

const PROJECT_DIR  = "/projappl/project_2011796/cascadia"
const CHAIN_GLOBS  = ["run_*", "00runs/run_*"]   # dirs searched for best_model.rho
const MODEL_NAME   = "best_model.rho"
const LOG_NAME     = "0vfsa3DMT.log"

# convergence gates -- chains at wildly different iteration counts are not
# comparable, and a stalled chain drags the ensemble mean toward the prior
const MIN_ITER     = 0      # drop chains that ran fewer iterations than this
const RMS_TOL      = Inf    # drop chains with best_rms > RMS_TOL * ensemble_best

const GRID_ATOL    = 1e-8
const GRID_RTOL    = 1e-8

#---------- ws3d io ----------

struct WS3D
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    A::Array{Float64,3}          # LOGE
    origin::Vector{Float64}
    rotation::Float64
    type_str::String
end

function read_ws3d(path::AbstractString)
    lines = readlines(path)
    i = 1
    while i <= length(lines) && (isempty(strip(lines[i])) || startswith(strip(lines[i]), "#"))
        i += 1
    end
    hdr = split(strip(lines[i]))
    nx, ny, nz = parse.(Int, hdr[1:3])
    type_str = length(hdr) >= 5 ? String(hdr[5]) : "LOGE"
    i += 1

    # everything numeric from here; the trailer (origin, rotation) is whatever
    # follows the nx*ny*nz block
    nums = Float64[]
    sizehint!(nums, nx + ny + nz + nx*ny*nz + 4)
    for ln in lines[i:end], t in split(ln)
        s = strip(t)
        isempty(s) && continue
        v = (uppercase(s) == "NAN" || s == "*") ? NaN : tryparse(Float64, s)
        v === nothing || push!(nums, v)
    end

    need = nx + ny + nz + nx*ny*nz
    length(nums) >= need || error("$(basename(path)): expected >= $need numbers, got $(length(nums))")

    dx = nums[1:nx]
    dy = nums[nx+1:nx+ny]
    dz = nums[nx+ny+1:nx+ny+nz]
    off = nx + ny + nz
    A = reshape(nums[off+1:off+nx*ny*nz], nx, ny, nz)

    tail = nums[need+1:end]
    origin = length(tail) >= 3 ? tail[1:3] : zeros(3)
    rotation = length(tail) >= 4 ? tail[4] : 0.0
    WS3D(dx, dy, dz, A, origin, rotation, type_str)
end

function write_ws3d(path::AbstractString, dx, dy, dz, A, origin, rotation, type_str, comment)
    nx, ny, nz = size(A)
    open(path, "w") do io
        println(io, "# ", comment)
        @printf(io, "%d %d %d 0 %s\n", nx, ny, nz, type_str)
        for v in (dx, dy, dz)
            println(io, join((@sprintf("%.6G", x) for x in v), " "))
        end
        for k in 1:nz
            println(io)
            for j in 1:ny
                println(io, join((@sprintf("%15.5E", A[i,j,k]) for i in 1:nx), ""))
            end
        end
        println(io, join((@sprintf("%.6G", o) for o in origin), " "))
        println(io, @sprintf("%.6G", rotation))
    end
end

#---------- chain discovery ----------

struct Chain
    dir::String
    model::String
    seed::String
    iter::Int
    best_rms::Float64
    last_rms::Float64
end

# scan a vfsa log for seed and the final iteration row
function parse_log(dir::AbstractString)
    path = joinpath(dir, LOG_NAME)
    seed = "?"; iter = 0; best = NaN; last = NaN
    isfile(path) || return (seed, iter, best, last)
    for ln in eachline(path)
        s = strip(ln)
        if startswith(s, "#")
            occursin("seed", s) && (parts = split(s); length(parts) >= 4 && (seed = parts[4]))
            continue
        end
        f = split(s)
        length(f) >= 4 || continue
        it = tryparse(Int, f[1]); it === nothing && continue
        r  = tryparse(Float64, f[3]); b = tryparse(Float64, f[4])
        (r === nothing || b === nothing) && continue
        iter = it; last = r; best = b
    end
    (seed, iter, best, last)
end

function expand_glob(root::AbstractString, pat::AbstractString)
    parts = split(pat, '/')
    dirs = [root]
    for p in parts
        nxt = String[]
        for d in dirs
            isdir(d) || continue
            if occursin('*', p)
                rx = Regex("^" * replace(replace(p, "." => "\\."), "*" => ".*") * "\$")
                for e in readdir(d)
                    occursin(rx, e) && isdir(joinpath(d, e)) && push!(nxt, joinpath(d, e))
                end
            else
                isdir(joinpath(d, p)) && push!(nxt, joinpath(d, p))
            end
        end
        dirs = nxt
    end
    dirs
end

function discover_chains(root::AbstractString)
    seen = Set{String}()
    out = Chain[]
    for pat in CHAIN_GLOBS, d in expand_glob(root, pat)
        d in seen && continue
        push!(seen, d)
        m = joinpath(d, MODEL_NAME)
        isfile(m) || continue
        seed, iter, best, last = parse_log(d)
        push!(out, Chain(d, m, seed, iter, best, last))
    end
    sort!(out, by = c -> c.dir)
    out
end

#---------- grid checks ----------

maxabsdiff(a, b) = length(a) == length(b) ? maximum(abs.(a .- b)) : Inf

function grids_match(ref::WS3D, m::WS3D)
    size(ref.A) == size(m.A) || return "size $(size(m.A)) != $(size(ref.A))"
    for (nm, a, b) in (("dx", ref.dx, m.dx), ("dy", ref.dy, m.dy), ("dz", ref.dz, m.dz))
        d = maxabsdiff(a, b)
        d <= GRID_ATOL + GRID_RTOL * maximum(abs.(a)) || return "$nm differs by $d"
    end
    maximum(abs.(ref.origin .- m.origin)) <= 1e-6 || return "origin differs"
    ""
end

#---------- statistics ----------

# linear-interpolated quantile on an already-sorted view
function quantile_sorted(v::AbstractVector{Float64}, n::Int, q::Float64)
    n == 1 && return v[1]
    h = (n - 1) * q + 1
    lo = floor(Int, h); hi = min(lo + 1, n)
    v[lo] + (h - lo) * (v[hi] - v[lo])
end

function ensemble_stats(cubes::Vector{Array{Float64,3}})
    nx, ny, nz = size(cubes[1])
    nc = length(cubes)
    mean_  = Array{Float64}(undef, nx, ny, nz)
    med_   = Array{Float64}(undef, nx, ny, nz)
    std_   = Array{Float64}(undef, nx, ny, nz)
    p16_   = Array{Float64}(undef, nx, ny, nz)
    p84_   = Array{Float64}(undef, nx, ny, nz)
    vals   = Vector{Float64}(undef, nc)

    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        cnt = 0
        for c in 1:nc
            v = cubes[c][i,j,k]
            isfinite(v) && (cnt += 1; vals[cnt] = v)
        end
        if cnt == 0
            mean_[i,j,k] = med_[i,j,k] = std_[i,j,k] = p16_[i,j,k] = p84_[i,j,k] = NaN
            continue
        end
        μ = sum(view(vals, 1:cnt)) / cnt
        mean_[i,j,k] = μ
        std_[i,j,k]  = cnt == 1 ? 0.0 :
                       sqrt(sum((view(vals, 1:cnt) .- μ).^2) / (cnt - 1))
        sort!(view(vals, 1:cnt))
        med_[i,j,k] = quantile_sorted(vals, cnt, 0.5)
        p16_[i,j,k] = quantile_sorted(vals, cnt, 0.16)
        p84_[i,j,k] = quantile_sorted(vals, cnt, 0.84)
    end
    (mean = mean_, median = med_, std = std_, p16 = p16_, p84 = p84_)
end

#---------- driver ----------

function main(args)
    root = PROJECT_DIR
    min_iter = MIN_ITER
    rms_tol  = RMS_TOL
    for a in args
        if startswith(a, "--min-iter=")
            min_iter = parse(Int, split(a, '=')[2])
        elseif startswith(a, "--rms-tol=")
            rms_tol = parse(Float64, split(a, '=')[2])
        elseif !startswith(a, "--")
            root = abspath(a)
        end
    end
    isdir(root) || error("project dir not found: $root")

    chains = discover_chains(root)
    globs = join(CHAIN_GLOBS, ", ")
    isempty(chains) && error("no $MODEL_NAME found under $root in $globs")

    println("project: ", root)
    println("found ", length(chains), " chain(s) with a $MODEL_NAME\n")

    # gate on convergence before loading anything heavy
    finite = filter(c -> isfinite(c.best_rms), chains)
    ens_best = isempty(finite) ? NaN : minimum(c.best_rms for c in finite)
    keep = Chain[]; drop = Tuple{Chain,String}[]
    for c in chains
        if c.iter < min_iter
            push!(drop, (c, "iter $(c.iter) < $min_iter"))
        elseif isfinite(rms_tol) && isfinite(c.best_rms) && isfinite(ens_best) &&
               c.best_rms > rms_tol * ens_best
            push!(drop, (c, @sprintf("rms %.3f > %.2f x %.3f", c.best_rms, rms_tol, ens_best)))
        else
            push!(keep, c)
        end
    end

    @printf("%-28s %10s %8s %10s %10s\n", "chain", "seed", "iter", "best_rms", "last_rms")
    println("-"^70)
    for c in chains
        mark = any(d -> d[1].dir == c.dir, drop) ? "  DROPPED" : ""
        @printf("%-28s %10s %8d %10.4f %10.4f%s\n",
                basename(c.dir), c.seed, c.iter, c.best_rms, c.last_rms, mark)
    end
    for (c, why) in drop
        println("  dropped ", basename(c.dir), ": ", why)
    end
    println()

    length(keep) >= 2 || error("need >= 2 chains after filtering, have $(length(keep))")

    # load, skipping anything a live chain is mid-write on
    models = WS3D[]; used = Chain[]
    for c in keep
        try
            push!(models, read_ws3d(c.model)); push!(used, c)
        catch e
            @warn "skipping $(basename(c.dir)): $(sprint(showerror, e))"
        end
    end
    length(models) >= 2 || error("need >= 2 readable models, have $(length(models))")

    ref = models[1]
    for (m, c) in zip(models[2:end], used[2:end])
        why = grids_match(ref, m)
        isempty(why) || error("grid mismatch in $(basename(c.dir)): $why")
    end

    st = ensemble_stats([m.A for m in models])
    nx, ny, nz = size(ref.A)

    stamp = readchomp(`date +%Y%m%d_%H%M%S`)
    outdir = joinpath(root, "0stats_$stamp")
    mkpath(outdir)

    tag = "ensemble of $(length(models)) chains, $(basename(root))"
    for (name, A) in (("model.mean", st.mean), ("model.median", st.median),
                      ("model.std", st.std), ("model.p16", st.p16), ("model.p84", st.p84))
        write_ws3d(joinpath(outdir, name), ref.dx, ref.dy, ref.dz, A,
                   ref.origin, ref.rotation, ref.type_str, "$name -- $tag")
        println("wrote ", joinpath(outdir, name))
    end

    # per-layer agreement: where do the chains actually agree with depth?
    # std is LOGE; /log(10) puts it in decades, which is how it reads on a plot
    depth = cumsum(ref.dz)
    summary = IOBuffer()
    println(summary, "# ensemble stats -- ", tag)
    println(summary, "# generated ", stamp)
    println(summary, "#")
    println(summary, "# chains used:")
    for (m, c) in zip(models, used)
        @printf(summary, "#   %-26s seed %-10s iter %6d  best_rms %.4f\n",
                basename(c.dir), c.seed, c.iter, c.best_rms)
    end
    for (c, why) in drop
        @printf(summary, "#   %-26s DROPPED (%s)\n", basename(c.dir), why)
    end
    println(summary, "#")
    println(summary, "# per-layer inter-chain spread, whole grid, decades of log10(rho)")
    @printf(summary, "# %4s %12s %10s %10s %10s\n", "iz", "depth_m", "mean_std", "max_std", "p84-p16")
    for k in 1:nz
        s  = view(st.std, :, :, k)
        sp = view(st.p84, :, :, k) .- view(st.p16, :, :, k)
        fs = filter(isfinite, vec(s)); fp = filter(isfinite, vec(sp))
        isempty(fs) && continue
        @printf(summary, "  %4d %12.1f %10.4f %10.4f %10.4f\n",
                k, depth[k], sum(fs)/length(fs)/log(10),
                maximum(fs)/log(10), sum(fp)/length(fp)/log(10))
    end
    txt = String(take!(summary))
    write(joinpath(outdir, "0summary.txt"), txt)
    println("wrote ", joinpath(outdir, "0summary.txt"))
    println()
    print(txt)
end

# guard so this file can also be `include`d to reuse read_ws3d / write_ws3d
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

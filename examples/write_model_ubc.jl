using MTGeophysics
using Printf

function write_model_ubc(modfile::AbstractString, m::WS3DModel)
    nx, ny, nz = m.nx, m.ny, m.nz
    x_south = minimum(m.x)
    y_west  = minimum(m.y)
    z_top   = minimum(m.z)

    meshfile = modfile * ".mesh"
    confile  = modfile * ".con"

    open(meshfile, "w") do io
        @printf(io, "%d %d %d \n", ny, nx, nz)
        @printf(io, "%d %d %d \n", round(Int, y_west), round(Int, x_south), -round(Int, z_top))
        for i in 1:ny
            @printf(io, "%d ", round(Int, m.dy[i]))
        end
        @printf(io, "\n")
        for i in 1:nx
            @printf(io, "%d ", round(Int, m.dx[i]))
        end
        @printf(io, "\n")
        for i in 1:nz
            @printf(io, "%d ", round(Int, m.dz[i]))
        end
        @printf(io, "\n")
    end

    A = 10.0 .^ m.A
    A[isnan.(A)] .= 1.0e17
    C = 1.0 ./ A

    open(confile, "w") do io
        for i in 1:nx, j in 1:ny, k in 1:nz
            @printf(io, "%E\n", C[i, j, k])
        end
    end

    return meshfile, confile
end

function convert_models_ubc(outdir::AbstractString, models::AbstractVector{<:AbstractString})
    isdir(outdir) || mkpath(outdir)
    written = String[]
    for path in models
        isfile(path) || error("model file not found: $path")
        m = load_ws3d_model(path)
        stem = first(splitext(basename(path)))
        meshfile, confile = write_model_ubc(joinpath(outdir, stem), m)
        nair = count(isnan, m.A)
        @printf("%-18s %3d x %3d x %3d   %7d cells (%6d air)\n",
                stem, m.nx, m.ny, m.nz, m.nx * m.ny * m.nz, nair)
        println("    ", meshfile)
        println("    ", confile)
        push!(written, meshfile, confile)
    end
    return written
end

const DEFAULT_ROOT   = "/projappl/project_2011796/MT3DINV4"
const DEFAULT_OUTDIR = joinpath(DEFAULT_ROOT, "ubc_models")
const DEFAULT_MODELS = [joinpath(DEFAULT_ROOT, "NLCG-3-1", "mF.rho"),
                        joinpath(DEFAULT_ROOT, "NLCG-3-1", "I_NLCG_232.rho")]

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_OUTDIR
    models = length(ARGS) >= 2 ? ARGS[2:end] : DEFAULT_MODELS
    convert_models_ubc(outdir, models)
end

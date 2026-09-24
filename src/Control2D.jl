# 2D control files
# Author: @pankajkmishra
# fwd.ctrl, inv.ctrl (GN, NLCG), the VFSA control, the model covariance and mask.ctrl of the 2D workflows
# Control files are "Key : value" lines; unknown, duplicate and missing required keys are errors

using Printf

#---------- key-value parsing ----------

_ctrl_key(s) = lowercase(join(split(strip(s)), " "))

function _ctrl_parse_bool(s)
    v = lowercase(strip(s))
    v in ("yes", "true", "1") && return true
    v in ("no", "false", "0") && return false
    throw(ArgumentError("expected yes or no, got '$s'"))
end

function _ctrl_parse_pair(s)
    t = split(strip(s))
    length(t) == 2 || throw(ArgumentError("expected two numbers, got '$s'"))
    (parse(Float64, t[1]), parse(Float64, t[2]))
end

# spec: display key => (field, parser, required)
function _read_ctrl(path::AbstractString, spec, kind::AbstractString)
    isfile(path) || error("$kind file not found: $path")
    lookup = Dict(_ctrl_key(k) => (k, v...) for (k, v) in spec)
    values = Dict{Symbol, Any}()
    for (n, raw) in enumerate(eachline(path))
        line = strip(first(split(raw, '#')))
        isempty(line) && continue
        occursin(':', line) || error("$path:$n: expected 'Key : value', got '$line'")
        key, value = strip.(split(line, ':'; limit = 2))
        entry = get(lookup, _ctrl_key(key), nothing)
        entry === nothing && error("$path:$n: unknown $kind key '$key'")
        _, field, parser, _ = entry
        haskey(values, field) && error("$path:$n: duplicate $kind key '$key'")
        values[field] = try
            parser(value)
        catch err
            error("$path:$n: bad value for '$key': $(sprint(showerror, err))")
        end
    end
    for (k, (field, _, required)) in spec
        required && !haskey(values, field) && error("$path: missing required $kind key '$k'")
    end
    values
end

function _write_ctrl(io, rows)
    width = maximum(length ∘ first, rows)
    foreach(((k, v),) -> println(io, rpad(k, width), " : ", v), rows)
end

_ctrl_mode(s) = (m = Symbol(uppercase(strip(s))); m in (:TE, :TM, :TETM) ? m :
                 throw(ArgumentError("mode must be TE, TM or TETM")))
_ctrl_algorithm(s) = (a = Symbol(lowercase(strip(s))); a in (:gn, :nlcg) ? a :
                      throw(ArgumentError(a == :vfsa ? "VFSA runs through VFSA2D with its own control and mask.ctrl" :
                                          "algorithm must be GN or NLCG")))
_ctrl_yesno(b::Bool) = b ? "yes" : "no"

#---------- fwd.ctrl ----------

"""
    FwdCtrl2D

Forward control, read from `fwd.ctrl`. There are no built-in defaults for the air.
- `mode`: `:TE`, `:TM` or `:TETM`
- `air_layers`, `air_thickness`: number and total thickness (m) of the air layers
- `air_growth`: thickness ratio between successive air layers going up, 1 = uniform
- `air_resistivity`: ohm metres
- `write_frechet`: also write G next to the predicted data
- `dipole_length`: electric dipole length (m); next to a topographic step TM Ey is
  averaged over it, flat stations are unaffected
"""
Base.@kwdef struct FwdCtrl2D
    mode::Symbol
    air_layers::Int
    air_thickness::Float64
    air_growth::Float64
    air_resistivity::Float64
    write_frechet::Bool = false
    dipole_length::Float64 = 100.0
end

const _FWD_CTRL_SPEC = (
    "Mode"                     => (:mode, _ctrl_mode, true),
    "Air layers"               => (:air_layers, s -> parse(Int, s), true),
    "Air thickness (m)"        => (:air_thickness, s -> parse(Float64, s), true),
    "Air growth factor"        => (:air_growth, s -> parse(Float64, s), true),
    "Air resistivity (ohm m)"  => (:air_resistivity, s -> parse(Float64, s), true),
    "Write Frechet derivative" => (:write_frechet, _ctrl_parse_bool, false),
    "Dipole length (m)"        => (:dipole_length, s -> parse(Float64, s), false),
)

function _validate_ctrl(c::FwdCtrl2D)
    c.air_layers > 0 || throw(ArgumentError("air layers must be positive"))
    c.air_thickness > 0 || throw(ArgumentError("air thickness must be positive"))
    c.air_growth >= 1 || throw(ArgumentError("air growth factor must be at least 1"))
    c.air_resistivity > 0 || throw(ArgumentError("air resistivity must be positive"))
    c.dipole_length >= 0 || throw(ArgumentError("dipole length must be nonnegative"))
    c
end

"""
    ReadFwdCtrl2D(path) -> FwdCtrl2D

Read a 2D forward control file. Every air key is required.
"""
ReadFwdCtrl2D(path::AbstractString) = _validate_ctrl(FwdCtrl2D(; _read_ctrl(path, _FWD_CTRL_SPEC, "fwd.ctrl")...))

"""
    WriteFwdCtrl2D(path, ctrl::FwdCtrl2D) -> path

Write a 2D forward control file.
"""
function WriteFwdCtrl2D(path::AbstractString, c::FwdCtrl2D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        _write_ctrl(io, [
            ("Mode", string(c.mode)),
            ("Air layers", string(c.air_layers)),
            ("Air thickness (m)", @sprintf("%.6g", c.air_thickness)),
            ("Air growth factor", @sprintf("%.6g", c.air_growth)),
            ("Air resistivity (ohm m)", @sprintf("%.6g", c.air_resistivity)),
            ("Write Frechet derivative", _ctrl_yesno(c.write_frechet)),
            ("Dipole length (m)", @sprintf("%.6g", c.dipole_length)),
        ])
    end
    String(path)
end

#---------- inv.ctrl ----------

"""
    InvCtrl2D

Deterministic inversion control, read from `inv.ctrl`. The lambda, rms and iteration
keys follow ModEM's NLCG control file; the rest are MTGeophysics extras with defaults.
VFSA has its own control, `VFSACtrl2D`.
- `algorithm`: `:gn` or `:nlcg`
- `lambda`: regularization weight β, fixed through the run
- `target_rms`, `max_iter`: stopping controls
- `mode`, `log_bounds`, `max_step`, `max_linesearch`: shared driver controls
- `smallness`, `smooth_y`, `smooth_z`: weights of the gradient regularizer
- `gn_damping`, `nlcg_restart`, `nlcg_precondition`: per-algorithm settings
"""
Base.@kwdef struct InvCtrl2D
    algorithm::Symbol
    lambda::Float64
    target_rms::Float64
    max_iter::Int
    mode::Symbol = :TETM
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    max_step::Float64 = 0.5
    max_linesearch::Int = 12
    smallness::Float64 = 0.01
    smooth_y::Float64 = 1.0
    smooth_z::Float64 = 1.0
    gn_damping::Float64 = 1e-2
    nlcg_restart::Int = 30
    nlcg_precondition::Bool = true
end

_ctrl_float(s) = parse(Float64, s)
_ctrl_int(s) = parse(Int, s)

const _INV_CTRL_SPEC = (
    "Algorithm"                          => (:algorithm, _ctrl_algorithm, true),
    "Initial damping factor lambda"      => (:lambda, _ctrl_float, true),
    "Exit search when rms is less than"  => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"       => (:max_iter, _ctrl_int, true),
    "Mode"                               => (:mode, _ctrl_mode, false),
    "Log10 resistivity bounds"           => (:log_bounds, _ctrl_parse_pair, false),
    "Max log10 step"                     => (:max_step, _ctrl_float, false),
    "Max line search steps"              => (:max_linesearch, _ctrl_int, false),
    "Smallness weight"                   => (:smallness, _ctrl_float, false),
    "Smoothing weight y"                 => (:smooth_y, _ctrl_float, false),
    "Smoothing weight z"                 => (:smooth_z, _ctrl_float, false),
    "GN damping"                         => (:gn_damping, _ctrl_float, false),
    "NLCG restart"                       => (:nlcg_restart, _ctrl_int, false),
    "NLCG precondition"                  => (:nlcg_precondition, _ctrl_parse_bool, false),
)

function _validate_ctrl(c::InvCtrl2D)
    c.algorithm in (:gn, :nlcg) || throw(ArgumentError("algorithm must be GN or NLCG"))
    c.lambda >= 0 || throw(ArgumentError("lambda must be nonnegative"))
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c.log_bounds[1] < c.log_bounds[2] || throw(ArgumentError("log10 bounds must be increasing"))
    c
end

"""
    ReadInvCtrl2D(path) -> InvCtrl2D

Read a GN or NLCG inversion control file.
"""
ReadInvCtrl2D(path::AbstractString) = _validate_ctrl(InvCtrl2D(; _read_ctrl(path, _INV_CTRL_SPEC, "inv.ctrl")...))

"""
    WriteInvCtrl2D(path, ctrl::InvCtrl2D) -> path

Write a GN or NLCG inversion control file with the shared keys and those of its
algorithm; the other algorithm's keys are left out and read back as defaults.
"""
function WriteInvCtrl2D(path::AbstractString, c::InvCtrl2D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g(x) = @sprintf("%.6g", x)
    rows = [
        ("Algorithm", uppercase(string(c.algorithm))),
        ("Initial damping factor lambda", g(c.lambda)),
        ("Exit search when rms is less than", g(c.target_rms)),
        ("Maximum number of iterations", string(c.max_iter)),
        ("Mode", string(c.mode)),
        ("Log10 resistivity bounds", "$(g(c.log_bounds[1])) $(g(c.log_bounds[2]))"),
        ("Max log10 step", g(c.max_step)),
        ("Max line search steps", string(c.max_linesearch)),
        ("Smallness weight", g(c.smallness)),
        ("Smoothing weight y", g(c.smooth_y)),
        ("Smoothing weight z", g(c.smooth_z)),
    ]
    c.algorithm == :gn ? push!(rows, ("GN damping", g(c.gn_damping))) :
        append!(rows, [("NLCG restart", string(c.nlcg_restart)),
                       ("NLCG precondition", _ctrl_yesno(c.nlcg_precondition))])
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

#---------- vfsa control ----------

"""
    VFSACtrl2D

VFSA control, read from `InvCtrl.VFSA`. No regularization and no prior: the start model
is the centre of the search and `mask.ctrl` says which cells move.
- `target_rms`, `max_iter`: a chain stops at the target rms or after `max_iter` iterations
- `mode`, `log_bounds`: data mode and the log10 ρ box
- `chains`, `control_points`, `trials`: independent chains, RBF control points per chain,
  proposals per iteration
- `step_scale`, `temperature`, `cooling`: proposal width as a share of the box, starting
  temperature and its ratio at `max_iter`
- `rbf_y`, `rbf_z`: RBF widths in cells
- `seed`, `snapshots`: random seed, best-model snapshot interval (0 = off)
"""
Base.@kwdef struct VFSACtrl2D
    target_rms::Float64
    max_iter::Int
    mode::Symbol = :TETM
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    chains::Int = 2
    control_points::Int = 400
    trials::Int = 1
    step_scale::Float64 = 0.11
    temperature::Float64 = 1.0
    cooling::Float64 = 1e-3
    rbf_y::Float64 = 2.0
    rbf_z::Float64 = 2.5
    seed::Int = 20260308
    snapshots::Int = 0
end

const _VFSA_CTRL_SPEC = (
    "Exit search when rms is less than" => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"      => (:max_iter, _ctrl_int, true),
    "Mode"                              => (:mode, _ctrl_mode, false),
    "Log10 resistivity bounds"          => (:log_bounds, _ctrl_parse_pair, false),
    "Number of chains"                  => (:chains, _ctrl_int, false),
    "Control points"                    => (:control_points, _ctrl_int, false),
    "Trials per iteration"              => (:trials, _ctrl_int, false),
    "Step scale"                        => (:step_scale, _ctrl_float, false),
    "Starting temperature"              => (:temperature, _ctrl_float, false),
    "Cooling ratio"                     => (:cooling, _ctrl_float, false),
    "RBF width y (cells)"               => (:rbf_y, _ctrl_float, false),
    "RBF width z (cells)"               => (:rbf_z, _ctrl_float, false),
    "Random seed"                       => (:seed, _ctrl_int, false),
    "Snapshot interval"                 => (:snapshots, _ctrl_int, false),
)

function _validate_ctrl(c::VFSACtrl2D)
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c.log_bounds[1] < c.log_bounds[2] || throw(ArgumentError("log10 bounds must be increasing"))
    c.chains >= 1 && c.control_points >= 1 && c.trials >= 1 ||
        throw(ArgumentError("chains, control points and trials must be at least 1"))
    c.temperature > 0 && 0 < c.cooling <= 1 || throw(ArgumentError("need temperature > 0 and 0 < cooling ratio ≤ 1"))
    c
end

"""
    ReadVFSACtrl2D(path) -> VFSACtrl2D

Read a VFSA control file.
"""
ReadVFSACtrl2D(path::AbstractString) = _validate_ctrl(VFSACtrl2D(; _read_ctrl(path, _VFSA_CTRL_SPEC, "VFSA control")...))

"""
    WriteVFSACtrl2D(path, ctrl::VFSACtrl2D) -> path

Write a VFSA control file.
"""
function WriteVFSACtrl2D(path::AbstractString, c::VFSACtrl2D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g(x) = @sprintf("%.6g", x)
    rows = [
        ("Exit search when rms is less than", g(c.target_rms)),
        ("Maximum number of iterations", string(c.max_iter)),
        ("Mode", string(c.mode)),
        ("Log10 resistivity bounds", "$(g(c.log_bounds[1])) $(g(c.log_bounds[2]))"),
        ("Number of chains", string(c.chains)),
        ("Control points", string(c.control_points)),
        ("Trials per iteration", string(c.trials)),
        ("Step scale", g(c.step_scale)),
        ("Starting temperature", g(c.temperature)),
        ("Cooling ratio", g(c.cooling)),
        ("RBF width y (cells)", g(c.rbf_y)),
        ("RBF width z (cells)", g(c.rbf_z)),
        ("Random seed", string(c.seed)),
        ("Snapshot interval", string(c.snapshots)),
    ]
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

#---------- model covariance ----------

"""
    Cov2D

Model covariance, read from the covariance file. The smoothing values keep the
ModEM layout; the inversion reads only the mask.
- `sy`: horizontal smoothing, one value per earth layer, top to bottom
- `sz`: vertical smoothing
- `n_smooth`: number of smoothing passes
- `exceptions`: `(class_a, class_b, smoothing)` rules between mask classes
- `mask`: `(nz, ny)` earth cells, top row first; 0 = air or fixed, 1 = free,
  9 = water (fixed, no regularization), other values = exception classes, inverted like 1
"""
Base.@kwdef struct Cov2D
    sy::Vector{Float64}
    sz::Float64
    n_smooth::Int
    exceptions::Vector{Tuple{Int, Int, Float64}} = Tuple{Int, Int, Float64}[]
    mask::Matrix{Int}
end

function _validate_ctrl(c::Cov2D)
    nz, ny = size(c.mask)
    length(c.sy) == nz || throw(ArgumentError("cov: $(length(c.sy)) y smoothing values for $nz layers"))
    all(s -> 0 <= s <= 1, c.sy) && 0 <= c.sz <= 1 || throw(ArgumentError("cov: smoothing must lie in [0, 1]"))
    c.n_smooth >= 0 || throw(ArgumentError("cov: number of smoothing passes must be nonnegative"))
    all(>=(0), c.mask) || throw(ArgumentError("cov: mask values must be nonnegative"))
    c
end

"""
    Cov2D(nz, ny; smoothing=0.3, n_smooth=1)

Uniform covariance for an `nz × ny` earth grid, every cell free.
"""
Cov2D(nz::Integer, ny::Integer; smoothing::Real = 0.3, n_smooth::Integer = 1) =
    Cov2D(sy = fill(Float64(smoothing), nz), sz = Float64(smoothing), n_smooth = n_smooth, mask = ones(Int, nz, ny))

"""
    ReadCov2D(path) -> Cov2D

Read a 2D model covariance file. The numbers are read in order: `ny nz`, `nz` y
smoothings, the z smoothing, the number of passes, the number of exception rules and
the rules, then `nz` mask rows of `ny` integers.
"""
function ReadCov2D(path::AbstractString)
    isfile(path) || error("covariance file not found: $path")
    tokens = String[]
    for raw in eachline(path)
        append!(tokens, split(first(split(raw, '#'))))
    end
    k = 0
    next() = (k += 1; k <= length(tokens) || error("$path: covariance file ends early"); tokens[k])
    ny, nz = parse(Int, next()), parse(Int, next())
    sy = [parse(Float64, next()) for _ in 1:nz]
    sz = parse(Float64, next())
    n_smooth = parse(Int, next())
    exceptions = [(parse(Int, next()), parse(Int, next()), parse(Float64, next())) for _ in 1:parse(Int, next())]
    mask = Matrix{Int}(undef, nz, ny)
    for iz in 1:nz, iy in 1:ny
        mask[iz, iy] = parse(Int, next())
    end
    k == length(tokens) || error("$path: $(length(tokens) - k) extra values after the mask")
    _validate_ctrl(Cov2D(; sy, sz, n_smooth, exceptions, mask))
end

"""
    WriteCov2D(path, cov::Cov2D) -> path

Write a 2D model covariance file.
"""
function WriteCov2D(path::AbstractString, c::Cov2D)
    _validate_ctrl(c)
    nz, ny = size(c.mask)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        @printf(io, "%d %d\n", ny, nz)
        _write_vector_lines(io, c.sy)
        @printf(io, "%.6g\n", c.sz)
        println(io, c.n_smooth)
        println(io, length(c.exceptions))
        for (a, b, s) in c.exceptions
            @printf(io, "%d %d %.6g\n", a, b, s)
        end
        for iz in 1:nz
            println(io, join(c.mask[iz, :], ' '))
        end
    end
    String(path)
end

#---------- mask.ctrl ----------

"""
    ReadMask2D(path) -> Matrix{Int}

Read a `mask.ctrl`: `ny nz`, then `nz` rows of `ny` integers, top row first, with the
covariance mask's codes (0 = air or fixed, 9 = water, others free). VFSA reads it in
place of the covariance file.
"""
function ReadMask2D(path::AbstractString)
    isfile(path) || error("mask file not found: $path")
    tokens = String[]
    for raw in eachline(path)
        append!(tokens, split(first(split(raw, '#'))))
    end
    length(tokens) >= 2 || error("$path: expected ny nz and the mask")
    ny, nz = parse(Int, tokens[1]), parse(Int, tokens[2])
    length(tokens) == 2 + ny * nz || error("$path: expected $(ny * nz) mask values for $ny × $nz cells, found $(length(tokens) - 2)")
    mask = permutedims(reshape(parse.(Int, tokens[3:end]), ny, nz))
    all(>=(0), mask) || error("$path: mask values must be nonnegative")
    mask
end

"""
    WriteMask2D(path, mask) -> path

Write a `(nz, ny)` integer mask as `mask.ctrl`.
"""
function WriteMask2D(path::AbstractString, mask::AbstractMatrix{<:Integer})
    all(>=(0), mask) || throw(ArgumentError("mask values must be nonnegative"))
    nz, ny = size(mask)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        @printf(io, "%d %d\n", ny, nz)
        foreach(iz -> println(io, join(mask[iz, :], ' ')), 1:nz)
    end
    String(path)
end

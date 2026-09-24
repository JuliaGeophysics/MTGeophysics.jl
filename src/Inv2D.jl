# 2D MT deterministic inversion framework
# Author: @pankajkmishra
# Algorithm-agnostic problem setup, objective, regularization, line search, and driver loop
# Each algorithm lives in its own Inv2D_<Name>.jl and plugs in through the interface below

using LinearAlgebra
using Printf
using SparseArrays
using Statistics

#---------- algorithm interface ----------
#
# an algorithm is a config struct `A <: AbstractInversion2D` plus these methods:
#
#   inv2d_init(alg, problem, state)                    -> work   (mutable scratch)
#   inv2d_prepare!(alg, work, problem, state)          -> gradient, once per iteration
#   inv2d_direction(alg, work, problem, state, grad)   -> search direction in log10 ρ
#   inv2d_reject!(alg, work)                           -> true = retry with a new direction
#   inv2d_accept!(alg, work, problem, old, new)        -> nothing
#   inv2d_info(alg, work)                              -> NamedTuple added to history
#   inv2d_tag(alg)                                     -> short label for the log, e.g. "gn"
#   inv2d_validate(alg)                                -> throws on bad settings
#
# the driver owns stopping tests, step capping and the Armijo line search,
# so a new algorithm only decides which direction to take

"""
    AbstractInversion2D

Supertype of 2D deterministic inversion algorithms. See `Inv2D.jl` for the interface.
"""
abstract type AbstractInversion2D end

inv2d_reject!(::AbstractInversion2D, work) = false
inv2d_accept!(::AbstractInversion2D, work, problem, old, new) = nothing
inv2d_info(::AbstractInversion2D, work) = (;)
inv2d_validate(::AbstractInversion2D) = nothing

#---------- options and result ----------

"""
    Inv2DOptions(; kwargs...)

Algorithm-independent inversion controls.
- `mode`: `:TE`, `:TM`, or `:TETM` impedances to fit
- `max_iter`: maximum accepted iterations; 0 evaluates the start model only
- `beta`: regularization weight about the reference model
- `smallness`, `smooth_y`, `smooth_z`: weights of the regularization terms
- `max_step`: cap on the largest log10 change per iteration
- `max_linesearch`: backtracking halvings per direction
- `target_rms`: stop at this rms; 0 = off
- `gradient_tolerance`, `step_tolerance`, `objective_tolerance`: remaining stopping tests
- `verbose`: print one line per accepted iteration
"""
Base.@kwdef struct Inv2DOptions
    mode::Symbol = :TETM
    max_iter::Int = 15
    beta::Float64 = 1.0
    smallness::Float64 = 0.01
    smooth_y::Float64 = 1.0
    smooth_z::Float64 = 1.0
    max_step::Float64 = 0.5
    max_linesearch::Int = 12
    target_rms::Float64 = 1.0
    gradient_tolerance::Float64 = 1e-6
    step_tolerance::Float64 = 1e-5
    objective_tolerance::Float64 = 1e-6
    verbose::Bool = true
end

"""
    Inv2DResult

Result of `Invert2D`.
- `resistivity`: recovered model in ohm metres, air included
- `response`, `fit`: predicted response and its chi2/rms
- `history`: accepted iterations, iteration zero first
- `active_cells`: inverted cells in model `(z, y)` order
- `converged`: false for `:max_iter` and `:line_search_failed`
- `reason`: termination reason
- `algorithm`: the algorithm config that produced it
"""
struct Inv2DResult{A<:AbstractInversion2D}
    resistivity::Matrix{Float64}
    response::MT2DResponse
    fit::FitSummary2D
    history::Vector{NamedTuple}
    active_cells::Vector{CartesianIndex{2}}
    converged::Bool
    reason::Symbol
    algorithm::A
end

#---------- problem ----------

"""
    Inv2DProblem

Fixed part of an inversion: mesh, weighted data rows, active cells, regularization
operator, reference model, and options. Built once by `Invert2D`.
"""
struct Inv2DProblem
    mesh::MT2DMesh
    rows::Vector{NamedTuple}
    cells::Vector{CartesianIndex{2}}
    R::SparseMatrixCSC{Float64, Int}          # all model cells
    R_active::SparseMatrixCSC{Float64, Int}   # active columns only
    mref::Matrix{Float64}
    options::Inv2DOptions
end

function _inv2d_validate(opt::Inv2DOptions)
    opt.max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    opt.mode in (:TE, :TM, :TETM) || throw(ArgumentError("invalid mode"))
    all(x -> isfinite(x) && x >= 0, (opt.beta, opt.smallness, opt.smooth_y, opt.smooth_z,
        opt.target_rms, opt.gradient_tolerance, opt.step_tolerance, opt.objective_tolerance)) ||
        throw(ArgumentError("regularization and stopping controls must be finite and nonnegative"))
    isfinite(opt.max_step) && opt.max_step > 0 || throw(ArgumentError("max_step must be finite and positive"))
    opt.max_linesearch > 0 || throw(ArgumentError("max_linesearch must be positive"))
    nothing
end

# active cells: default all earth cells; mask, cartesian, or linear indices otherwise
function _inv2d_cells(mesh, rho, active)
    indices = CartesianIndices(rho)
    air = mt2d_air_mask(mesh)
    cells = if active === nothing
        [i for i in indices if !air[i]]
    elseif active isa AbstractArray{Bool}
        size(active) == size(rho) || throw(DimensionMismatch("active mask must match model"))
        findall(active)
    else
        [i isa CartesianIndex{2} ? i : indices[i] for i in active]
    end
    isempty(cells) && throw(ArgumentError("at least one active earth cell is required"))
    all(i -> checkbounds(Bool, rho, i), cells) || throw(BoundsError(rho, cells))
    all(i -> !air[i], cells) || throw(ArgumentError("air cells cannot be inverted"))
    length(unique(cells)) == length(cells) || throw(ArgumentError("active cells must be unique"))
    cells
end

# one row per valid complex datum; the real and imaginary parts each use sigma,
# as chi2_rms2d does. invalid observations are dropped here once, never later
function _inv2d_data(mesh, data, mode)
    data.frequencies ≈ mesh.frequencies || throw(ArgumentError("data frequencies must match the mesh, including order"))
    data.receivers ≈ mesh.receiver_positions || throw(ArgumentError("data receivers must match the mesh, including order"))
    rows = NamedTuple[]
    dims = (length(mesh.frequencies), length(mesh.receiver_positions))
    for (pol, key, errkey) in ((:TE, :z_xy, :z_xy_error), (:TM, :z_yx, :z_yx_error))
        mode in (pol, :TETM) || continue
        z, err = getproperty(data, key), getproperty(data, errkey)
        size(z) == size(err) == dims || throw(DimensionMismatch("data/error arrays must match survey"))
        for i in eachindex(z)
            if isfinite(z[i]) && isfinite(err[i]) && err[i] > 0
                push!(rows, (; key, index = i, value = z[i], sigma = err[i]))
            end
        end
    end
    isempty(rows) && throw(ArgumentError("no valid impedance observations with positive errors"))
    rows
end

# finite-volume regularization of log10(m/mref): gradient terms approximate the
# area integral, smallness uses cell area over the median earth-cell area; excluded
# cells (air, water) take no term and no smoothing pair
function _inv2d_regularizer(mesh, rho, opt; excluded = mt2d_air_mask(mesh))
    nz, ny = size(rho)
    linear = LinearIndices(rho)
    I, J, V = Int[], Int[], Float64[]
    row = 0
    area0 = median(mesh.z_cell_sizes[mesh.n_air_cells+1:end]) * median(mesh.y_cell_sizes)
    for iy in 1:ny, iz in mesh.n_air_cells+1:nz
        excluded[iz, iy] && continue
        cell = linear[iz, iy]
        if opt.smallness > 0
            row += 1
            push!(I, row); push!(J, cell)
            push!(V, sqrt(opt.smallness * mesh.z_cell_sizes[iz] * mesh.y_cell_sizes[iy] / area0))
        end
        if iy < ny && opt.smooth_y > 0 && !excluded[iz, iy+1]
            row += 1
            w = sqrt(opt.smooth_y * mesh.z_cell_sizes[iz] / ((mesh.y_cell_sizes[iy] + mesh.y_cell_sizes[iy+1]) / 2))
            append!(I, (row, row)); append!(J, (cell, linear[iz, iy+1])); append!(V, (-w, w))
        end
        if iz < nz && opt.smooth_z > 0 && !excluded[iz+1, iy]
            row += 1
            w = sqrt(opt.smooth_z * mesh.y_cell_sizes[iy] / ((mesh.z_cell_sizes[iz] + mesh.z_cell_sizes[iz+1]) / 2))
            append!(I, (row, row)); append!(J, (cell, linear[iz+1, iy])); append!(V, (-w, w))
        end
    end
    sparse(I, J, V, row, length(rho))
end

#---------- state: one evaluated model ----------

function _inv2d_residual(response, rows)
    r = zeros(2 * length(rows))
    for (j, row) in enumerate(rows)
        delta = (getproperty(response, row.key)[row.index] - row.value) / row.sigma
        r[2j-1], r[2j] = real(delta), imag(delta)
    end
    r
end

# evaluate log10 model m on the active cells, other cells taken from rho0
function _inv2d_evaluate(problem::Inv2DProblem, rho0::AbstractMatrix, m::AbstractVector)
    opt = problem.options
    rho = copy(rho0)
    rho[problem.cells] = 10.0 .^ m
    fullm = log10.(rho)
    response, cache = _mt2d_forward_cache(problem.mesh, rho; mode = opt.mode)
    r = _inv2d_residual(response, problem.rows)
    reg = problem.R * vec(fullm - problem.mref)
    objective = (sum(abs2, r) + opt.beta * sum(abs2, reg)) / 2
    (; m = Vector{Float64}(m), rho, response, cache, r, reg, objective)
end

#---------- derivatives shared by all algorithms ----------

"""
    inv2d_frechet(problem, state; method=:auto)

Data-weighted Fréchet derivative C_D^{-1/2} G in log10 resistivity, `2 n_data × n_active`, one row pair
(real, imaginary) per datum.
- `:adjoint`: one transpose solve Gᵗ δd̂ per real datum, in that datum's frequency
- `:forward`: one tangent linear solve G δm per active cell over all frequencies and modes
- `:auto`: adjoint when there are fewer real data than active cells
"""
function inv2d_frechet(problem::Inv2DProblem, state; method::Symbol = :auto)
    method in (:auto, :adjoint, :forward) || throw(ArgumentError("method must be :auto, :adjoint, or :forward"))
    rows, cells = problem.rows, problem.cells
    if method == :adjoint || (method == :auto && 2 * length(rows) < length(cells))
        # real part: weight 1/σ, imaginary part: weight i/σ, since real(conj(i/σ) z) = imag(z)/σ
        adjoint_rows = NamedTuple[]
        for row in rows
            push!(adjoint_rows, (; row.key, row.index, weight = complex(1 / row.sigma)))
            push!(adjoint_rows, (; row.key, row.index, weight = 1im / row.sigma))
        end
        Gρ = _mt2d_frechet_rows(state.cache, adjoint_rows)
        linear = LinearIndices(state.rho)
        scale = [log(10.0) * state.rho[c] for c in cells]
        return Gρ[:, [linear[c] for c in cells]] .* scale'
    end
    G = zeros(2 * length(rows), length(cells))
    direction = zeros(size(state.rho))
    for (j, cell) in enumerate(cells)
        direction[cell] = log(10.0) * state.rho[cell]
        dz = _mt2d_frechet(state.cache, direction)
        for (i, row) in enumerate(rows)
            value = getproperty(dz, row.key)[row.index] / row.sigma
            G[2i-1, j], G[2i, j] = real(value), imag(value)
        end
        direction[cell] = 0
    end
    G
end

"""
    inv2d_gradient(problem, state[, G])

Gradient of the objective in active log10 resistivity. With `G` it is `Gᵗr + βRᵗreg`;
without it the data term comes from one adjoint solve per frequency and mode, which is
what gradient-only algorithms should use.
"""
function inv2d_gradient(problem::Inv2DProblem, state, G::AbstractMatrix)
    G' * state.r + problem.options.beta * (problem.R_active' * state.reg)
end

function inv2d_gradient(problem::Inv2DProblem, state)
    # d(chi2/2) = Σ real(conj(w) dz), w = (r_re + i r_im)/σ, the Gᵗ pairing
    zbar = (z_xy = zeros(ComplexF64, size(state.response.z_xy)),
            z_yx = zeros(ComplexF64, size(state.response.z_yx)))
    for (j, row) in enumerate(problem.rows)
        getproperty(zbar, row.key)[row.index] += complex(state.r[2j-1], state.r[2j]) / row.sigma
    end
    grho = _mt2d_frechet_transpose(state.cache, zbar)
    gdata = [grho[c] * log(10.0) * state.rho[c] for c in problem.cells]
    gdata + problem.options.beta * (problem.R_active' * state.reg)
end

#---------- line search ----------

# armijo backtracking along a capped direction, unbounded (bounds are a VFSA setting); nothing when
# no step is accepted
function _inv2d_linesearch(problem::Inv2DProblem, state, direction, gradient)
    opt = problem.options
    d = direction .* min(1.0, opt.max_step / max(norm(direction, Inf), eps()))
    for ls in 0:opt.max_linesearch-1
        alpha = 0.5^ls
        candidate = state.m + alpha .* d
        delta = candidate - state.m
        slope = dot(gradient, delta)
        slope < 0 || continue
        trial = _inv2d_evaluate(problem, state.rho, candidate)
        if isfinite(trial.objective) && trial.objective <= state.objective + 1e-4 * slope
            return trial, norm(delta, Inf), alpha
        end
    end
    nothing
end

#---------- driver ----------

"""
    Invert2D(mesh, initial_resistivity, observed::DataFile2D;
             algorithm=GaussNewton2DConfig(), options=Inv2DOptions(),
             active_cells=nothing, reference_resistivity=initial_resistivity)

Invert TE/TM complex impedances for earth-cell log10 resistivity, minimizing
`0.5*sum(abs2, Wd*(F(m)-d)) + 0.5*beta*sum(abs2, R*(m-mref))`.

- `algorithm`: `GaussNewton2DConfig`, `NLCG2DConfig` or any `AbstractInversion2D`
- `options`: `Inv2DOptions`, shared by every algorithm
- `active_cells`: model-shaped Boolean mask or vector of cartesian/linear indices;
  default all earth cells. Other cells stay fixed, air is always fixed
- `reference_resistivity`: model the regularization pulls toward; default the start
- `water_cells`: model-shaped Boolean mask of fixed water, left out of the default
  active cells; like air, it takes no regularization term and no smoothing pair

Returns `Inv2DResult`; no files are written by this method, the six-file method below
writes a run directory.
"""
function Invert2D(mesh::MT2DMesh, initial_resistivity::AbstractMatrix{<:Real}, observed::DataFile2D;
                  algorithm::AbstractInversion2D = GaussNewton2DConfig(),
                  options::Inv2DOptions = Inv2DOptions(),
                  active_cells = nothing, reference_resistivity = initial_resistivity,
                  water_cells::Union{Nothing, AbstractMatrix{Bool}} = nothing)
    opt = options
    _inv2d_validate(opt)
    inv2d_validate(algorithm)

    #---------- setup ----------
    rho0 = Matrix{Float64}(initial_resistivity)
    size(reference_resistivity) == size(rho0) || throw(DimensionMismatch("reference must match initial model"))
    air = mt2d_air_mask(mesh)
    all(x -> isfinite(x) && x > 0, reference_resistivity[.!air]) ||
        throw(ArgumentError("reference earth resistivity must be positive and finite"))
    cells = _inv2d_cells(mesh, rho0, active_cells)
    rows = _inv2d_data(mesh, observed, opt.mode)
    all(x -> isfinite(x) && x > 0, rho0[.!air]) ||
        throw(ArgumentError("initial earth resistivity must be positive and finite"))
    excluded = copy(air)
    if water_cells !== nothing
        size(water_cells) == size(rho0) || throw(DimensionMismatch("water mask must match initial model"))
        active_cells === nothing && filter!(i -> !water_cells[i], cells)
        any(i -> water_cells[i], cells) && throw(ArgumentError("water cells cannot be inverted"))
        excluded .|= water_cells
    end

    # air takes part in neither the parameters nor the regularization, water only in the forward
    rho0[air] .= mesh.air_resistivity
    reference = Matrix{Float64}(reference_resistivity)
    reference[air] .= mesh.air_resistivity
    R = _inv2d_regularizer(mesh, rho0, opt; excluded)
    R_active = R[:, LinearIndices(rho0)[cells]]
    problem = Inv2DProblem(mesh, rows, cells, R, R_active, log10.(reference), opt)

    m0 = log10.(rho0[cells])
    state = _inv2d_evaluate(problem, rho0, m0)
    all(isfinite, state.r) || error("initial forward response is not finite")
    work = inv2d_init(algorithm, problem, state)

    #---------- history ----------
    history = NamedTuple[]
    function record(iteration, step, alpha, info)
        chi2 = sum(abs2, state.r)
        entry = (; iteration, objective = state.objective, chi2, rms = sqrt(chi2 / length(state.r)),
                   regularization = sum(abs2, state.reg) / 2, step, alpha, info...)
        push!(history, entry)
        opt.verbose && @printf("%s %3d  RMS %.5g  objective %.6g  step %.3g\n",
                               uppercase(inv2d_tag(algorithm)), iteration, entry.rms, entry.objective, step)
        opt.verbose && flush(stdout)      # batch logs are block buffered otherwise
    end
    record(0, 0.0, 0.0, inv2d_info(algorithm, work))

    #---------- iterations ----------
    reason = :max_iter
    for iteration in 1:opt.max_iter
        if opt.target_rms > 0 && history[end].rms <= opt.target_rms
            reason = :target_rms
            break
        end
        gradient = inv2d_prepare!(algorithm, work, problem, state)
        if norm(gradient, Inf) <= opt.gradient_tolerance
            reason = :gradient_tolerance
            break
        end

        # the algorithm may retry with a new direction (e.g. more damping) after a failed search
        found = nothing
        while true
            direction = inv2d_direction(algorithm, work, problem, state, gradient)
            found = _inv2d_linesearch(problem, state, direction, gradient)
            found === nothing || break
            inv2d_reject!(algorithm, work) || break
        end
        if found === nothing
            reason = :line_search_failed
            break
        end

        old = state
        state, step, alpha = found
        record(iteration, step, alpha, inv2d_info(algorithm, work))
        inv2d_accept!(algorithm, work, problem, old, state)

        if step <= opt.step_tolerance
            reason = :step_tolerance
            break
        elseif old.objective - state.objective <= opt.objective_tolerance * max(old.objective, eps())
            reason = :objective_tolerance
            break
        end
    end
    if opt.target_rms > 0 && history[end].rms <= opt.target_rms
        reason = :target_rms
    end

    chi2 = sum(abs2, state.r)
    fit = FitSummary2D(chi2 = chi2, rms = sqrt(chi2 / length(state.r)), count = length(state.r))
    Inv2DResult(state.rho, state.response, fit, history, cells,
                reason ∉ (:max_iter, :line_search_failed), reason, algorithm)
end

#---------- file workflow ----------

# shared options and the algorithm config from inv.ctrl
function _inv2d_from_ctrl(c::InvCtrl2D; mode::Symbol = c.mode)
    options = Inv2DOptions(mode = mode, max_iter = c.max_iter, beta = c.lambda, smallness = c.smallness,
                           smooth_y = c.smooth_y, smooth_z = c.smooth_z,
                           max_step = c.max_step, max_linesearch = c.max_linesearch, target_rms = c.target_rms)
    algorithm = c.algorithm == :gn ? GaussNewton2DConfig(damping = c.gn_damping) :
                NLCG2DConfig(restart = c.nlcg_restart, precondition = c.nlcg_precondition)
    options, algorithm
end

# run_YYYYmmdd_HHMMSS next to the data, suffixed _2, _3 when started in the same second
function _inv2d_run_dir(data_path::AbstractString)
    base = joinpath(dirname(abspath(data_path)), "run_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    dir, k = base, 1
    while ispath(dir)
        k += 1
        dir = "$(base)_$k"
    end
    mkpath(dir)
    dir
end

_same_grid(a::ModelFile2D, b::ModelFile2D) =
    size(a.resistivity) == size(b.resistivity) && a.y_cell_sizes ≈ b.y_cell_sizes && a.z_cell_sizes ≈ b.z_cell_sizes

# predicted data on the observed survey, with the observed errors and coordinates
function _inv2d_predicted(response::MT2DResponse, observed::DataFile2D, mode::Symbol)
    predicted = data_from_response2d(response; z_xy_error = observed.z_xy_error, z_yx_error = observed.z_yx_error,
        site_names = observed.site_names, x_positions = observed.x_positions, z_positions = observed.z_positions,
        latitudes = observed.latitudes, longitudes = observed.longitudes, origin = observed.origin)
    mode == :TE && (predicted.z_yx .= NaN)
    mode == :TM && (predicted.z_xy .= NaN)
    predicted
end

# mesh, start model, active and water cells of the file workflows from the model and its mask
# (cov.ctrl or mask.ctrl): 0 = air or fixed, 9 = water, others free
function _inv2d_file_setup(start::ModelFile2D, observed::DataFile2D, fwd::FwdCtrl2D, mask::AbstractMatrix{<:Integer}, kind)
    size(mask) == size(start.resistivity) ||
        error("the $kind is $(size(mask)) cells, the model is $(size(start.resistivity))")
    mesh, ρ0 = Mesh2DFromInputs(start, observed, fwd)
    na = mesh.n_air_cells
    air = mt2d_air_mask(mesh)
    bad = findall(air[na+1:end, :] .& (mask .!= 0))
    isempty(bad) || error("$(length(bad)) topographic air cell(s) have $kind ≠ 0, first at (row, column) $(Tuple(bad[1]))")
    water = falses(size(ρ0))
    water[na+1:end, :] .= mask .== MT2D_MASK_WATER
    wet = intersect(mt2d_receiver_columns(mesh), findall(vec(any(water; dims = 1))))
    isempty(wet) || error("stations stand over water (mask $(MT2D_MASK_WATER)) in model column(s) $wet; place stations on land")
    active = falses(size(ρ0))
    active[na+1:end, :] .= (mask .!= 0) .& (mask .!= MT2D_MASK_WATER)
    active .&= .!air
    any(active) || error("the $kind fixes every cell")
    (; mesh, ρ0, active, water, offsets = mt2d_station_offsets(mesh, observed))
end

function _inv2d_write_history(path, history)
    open(path, "w") do io
        println(io, join(string.(keys(first(history))), ','))
        foreach(h -> println(io, join(values(h), ',')), history)
    end
    path
end

function _inv2d_open_run(run_dir, data_path, inputs)
    dir = run_dir === nothing ? _inv2d_run_dir(data_path) : mkpath(String(run_dir))
    mkpath(joinpath(dir, "inputs"))
    foreach(p -> cp(p, joinpath(dir, "inputs", basename(p)); force = true), inputs)
    println("Run directory: ", dir)
    dir
end

# cells and station snapping, the tail of every Summary.txt
function _inv2d_summary_cells(io, s)
    println(io, "Active cells: ", count(s.active))
    println(io, "Topographic air cells: ", sum(mt2d_topo_air(s.mesh)))
    println(io, "Water cells: ", count(s.water))
    any(o -> o.offset != 0, s.offsets) || return
    println(io, "Station snapping (Z = data depth below the model top, ground = mesh surface of the column, m):")
    @printf(io, "  %-12s %10s %10s %10s %10s\n", "site", "Z", "ground", "offset", "tolerance")
    for o in s.offsets
        @printf(io, "  %-12s %10.2f %10.2f %10.2f %10.2f%s\n", o.site, o.z, o.surface, o.offset, o.tolerance,
                abs(o.offset) > o.tolerance ? "  beyond tolerance" : "")
    end
end

_inv2d_components(mode) = mode == :TETM ? ["ZXY", "ZYX"] : mode == :TE ? ["ZXY"] : ["ZYX"]

"""
    Invert2D(start_path, data_path, fwd_path, inv_path, cov_path, prior_path; run_dir=nothing)

ModEM-style deterministic inversion (GN or NLCG) from six files: the start model, the
observed data, `fwd.ctrl`, `inv.ctrl`, the model covariance and the prior model. The
algorithm comes from `inv.ctrl`, the air from `fwd.ctrl`, and the inverted cells from the
covariance mask (0 = air or fixed, 9 = water, fixed and unregularized, others free).
Topographic air (model cells above 1e15 ohm m) must carry mask 0, and no station may
stand over water. The prior is the reference model of the regularization. VFSA has its
own five-file entry, `VFSA2D`.

Everything is written to `run_dir`, by default `run_YYYYmmdd_HHMMSS/` next to the data:
`model.rho` (ModEM layout, restartable), `data.pred`, `History.csv`, `Summary.txt`, and
copies of the inputs in `inputs/`. Returns a named tuple with the run directory, mesh,
models, data, history, rms and termination reason.
"""
function Invert2D(start_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString,
                  inv_path::AbstractString, cov_path::AbstractString, prior_path::AbstractString;
                  run_dir::Union{Nothing, AbstractString} = nothing)
    observed = load_data2d(data_path)
    fwd, ctrl, cov = ReadFwdCtrl2D(fwd_path), ReadInvCtrl2D(inv_path), ReadCov2D(cov_path)
    start, prior = ReadModel2D(start_path), ReadModel2D(prior_path)
    _same_grid(start, prior) || error("the prior model must be on the start model's grid")
    s = _inv2d_file_setup(start, observed, fwd, cov.mask, "covariance mask")
    priormesh, ρref = Mesh2DFromInputs(prior, observed, fwd; warn = false)
    mt2d_topo_air(priormesh) == mt2d_topo_air(s.mesh) || error("the prior model's topography differs from the start model's")
    dir = _inv2d_open_run(run_dir, data_path, (start_path, data_path, fwd_path, inv_path, cov_path, prior_path))

    options, algorithm = _inv2d_from_ctrl(ctrl)
    result = Invert2D(s.mesh, s.ρ0, observed; algorithm, options, active_cells = s.active, reference_resistivity = ρref,
                      water_cells = any(s.water) ? s.water : nothing)
    predicted = _inv2d_predicted(result.response, observed, ctrl.mode)
    fit = chi2_rms2d(observed, predicted; components = _inv2d_components(ctrl.mode))

    WriteModel2D(joinpath(dir, "model.rho"), s.mesh, result.resistivity)
    write_data2d(joinpath(dir, "data.pred"), predicted)
    _inv2d_write_history(joinpath(dir, "History.csv"), result.history)
    open(joinpath(dir, "Summary.txt"), "w") do io
        println(io, "Algorithm: ", uppercase(string(ctrl.algorithm)))
        println(io, "Termination: ", result.reason)
        println(io, "Converged: ", result.converged)
        @printf(io, "RMS: %.6f\n", fit.rms)
        println(io, "Real data count: ", fit.count)
        println(io, "Accepted iterations: ", length(result.history) - 1)
        _inv2d_summary_cells(io, s)
    end
    (; run_dir = dir, algorithm = ctrl.algorithm, ctrl, fwd, mesh = s.mesh, observed, predicted, start = s.ρ0,
       prior = ρref, final = result.resistivity, active = s.active, water = s.water, history = result.history,
       vfsa = nothing, rms = fit.rms, reason = result.reason, converged = result.converged)
end

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
# the driver owns stopping tests, step capping, bounds, and the Armijo line search,
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
- `log_bounds`: box bounds on active log10 resistivity
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
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
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
    lo, hi = opt.log_bounds
    isfinite(lo) && isfinite(hi) && -300 < lo < hi < 300 || throw(ArgumentError("invalid log_bounds"))
    nothing
end

# active cells: default all earth cells; mask, cartesian, or linear indices otherwise
function _inv2d_cells(mesh, rho, active)
    indices = CartesianIndices(rho)
    cells = if active === nothing
        [i for i in indices if i[1] > mesh.n_air_cells]
    elseif active isa AbstractArray{Bool}
        size(active) == size(rho) || throw(DimensionMismatch("active mask must match model"))
        findall(active)
    else
        [i isa CartesianIndex{2} ? i : indices[i] for i in active]
    end
    isempty(cells) && throw(ArgumentError("at least one active earth cell is required"))
    all(i -> checkbounds(Bool, rho, i), cells) || throw(BoundsError(rho, cells))
    all(i -> i[1] > mesh.n_air_cells, cells) || throw(ArgumentError("air cells cannot be inverted"))
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
# area integral, smallness uses cell area over the median earth-cell area
function _inv2d_regularizer(mesh, rho, opt)
    nz, ny = size(rho)
    linear = LinearIndices(rho)
    I, J, V = Int[], Int[], Float64[]
    row = 0
    area0 = median(mesh.z_cell_sizes[mesh.n_air_cells+1:end]) * median(mesh.y_cell_sizes)
    for iy in 1:ny, iz in mesh.n_air_cells+1:nz
        cell = linear[iz, iy]
        if opt.smallness > 0
            row += 1
            push!(I, row); push!(J, cell)
            push!(V, sqrt(opt.smallness * mesh.z_cell_sizes[iz] * mesh.y_cell_sizes[iy] / area0))
        end
        if iy < ny && opt.smooth_y > 0
            row += 1
            w = sqrt(opt.smooth_y * mesh.z_cell_sizes[iz] / ((mesh.y_cell_sizes[iy] + mesh.y_cell_sizes[iy+1]) / 2))
            append!(I, (row, row)); append!(J, (cell, linear[iz, iy+1])); append!(V, (-w, w))
        end
        if iz < nz && opt.smooth_z > 0
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

# projected armijo backtracking along a capped direction; nothing when no step is accepted
function _inv2d_linesearch(problem::Inv2DProblem, state, direction, gradient)
    opt = problem.options
    lo, hi = opt.log_bounds
    d = direction .* min(1.0, opt.max_step / max(norm(direction, Inf), eps()))
    for ls in 0:opt.max_linesearch-1
        alpha = 0.5^ls
        candidate = clamp.(state.m + alpha .* d, lo, hi)
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

Returns `Inv2DResult`; no files are written by this method, the six-file method below
writes a run directory.
"""
function Invert2D(mesh::MT2DMesh, initial_resistivity::AbstractMatrix{<:Real}, observed::DataFile2D;
                  algorithm::AbstractInversion2D = GaussNewton2DConfig(),
                  options::Inv2DOptions = Inv2DOptions(),
                  active_cells = nothing, reference_resistivity = initial_resistivity)
    opt = options
    _inv2d_validate(opt)
    inv2d_validate(algorithm)

    #---------- setup ----------
    rho0 = Matrix{Float64}(initial_resistivity)
    size(reference_resistivity) == size(rho0) || throw(DimensionMismatch("reference must match initial model"))
    all(x -> isfinite(x) && x > 0, reference_resistivity[mesh.n_air_cells+1:end, :]) ||
        throw(ArgumentError("reference earth resistivity must be positive and finite"))
    cells = _inv2d_cells(mesh, rho0, active_cells)
    rows = _inv2d_data(mesh, observed, opt.mode)
    all(x -> isfinite(x) && x > 0, rho0[mesh.n_air_cells+1:end, :]) ||
        throw(ArgumentError("initial earth resistivity must be positive and finite"))

    # air takes part in neither the parameters nor the regularization
    rho0[1:mesh.n_air_cells, :] .= mesh.air_resistivity
    reference = Matrix{Float64}(reference_resistivity)
    reference[1:mesh.n_air_cells, :] .= mesh.air_resistivity
    R = _inv2d_regularizer(mesh, rho0, opt)
    R_active = R[:, LinearIndices(rho0)[cells]]
    problem = Inv2DProblem(mesh, rows, cells, R, R_active, log10.(reference), opt)

    lo, hi = opt.log_bounds
    m0 = log10.(rho0[cells])
    all(x -> lo <= x <= hi, m0) || throw(ArgumentError("active starting resistivities are outside log_bounds"))
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
        if norm(state.m - clamp.(state.m - gradient, lo, hi), Inf) <= opt.gradient_tolerance
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
function _inv2d_from_ctrl(c::InvCtrl2D)
    options = Inv2DOptions(mode = c.mode, max_iter = c.max_iter, beta = c.lambda, smallness = c.smallness,
                           smooth_y = c.smooth_y, smooth_z = c.smooth_z, log_bounds = c.log_bounds,
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

# vfsa keeps its own driver and file layout until it moves behind this interface
function _invert2d_vfsa(dir, mesh, ρ0, observed, ctrl, fwd, cov)
    any(==(0), cov.mask) && @warn "VFSA does not use the covariance mask yet, every earth cell is perturbed"
    fwd.air_resistivity == 1e9 || @warn "VFSA uses 1e9 ohm m air, not the $(fwd.air_resistivity) in fwd.ctrl"
    vdir = joinpath(dir, "vfsa")
    mkpath(vdir)
    config = VFSA2DMTConfig(n_chains = ctrl.vfsa_chains, n_ctrl = ctrl.vfsa_control_points, max_iter = ctrl.max_iter,
                            n_trials = ctrl.vfsa_trials, log_bounds = ctrl.log_bounds, step_scale = ctrl.vfsa_step_scale,
                            seed = ctrl.vfsa_seed, target_rms = ctrl.target_rms, keep_models = true)
    vfsa = VFSA2DMT(write_model2d(joinpath(vdir, "start_with_air.rho"), mesh, ρ0),
                    write_data2d(joinpath(vdir, "observed.dat"), observed); run_dir = vdir, config)
    final = load_model2d(vfsa.best_chain.best_model_path).resistivity
    final, _inv2d_predicted(vfsa.best_chain.best_response, observed, ctrl.mode)
end

"""
    Invert2D(start_path, data_path, fwd_path, inv_path, cov_path, prior_path; run_dir=nothing)

ModEM-style inversion from six files: the start model, the observed data, `fwd.ctrl`,
`inv.ctrl`, the model covariance and the prior model. The algorithm (GN, NLCG or VFSA)
comes from `inv.ctrl`, the air from `fwd.ctrl`, and the inverted cells from the
covariance mask (0 = fixed). The prior is the reference model of the regularization.

Everything is written to `run_dir`, by default `run_YYYYmmdd_HHMMSS/` next to the data:
`model.rho` (ModEM layout, restartable), `data.pred`, `History.csv` (GN and NLCG),
`Summary.txt`, and copies of the inputs in `inputs/`. `examples/model_2D_to_SEGY.jl`
turns `model.rho` into a SEG-Y section. Returns a named tuple with the run directory,
mesh, models, data, history, rms and termination reason.
"""
function Invert2D(start_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString,
                  inv_path::AbstractString, cov_path::AbstractString, prior_path::AbstractString;
                  run_dir::Union{Nothing, AbstractString} = nothing)
    observed = load_data2d(data_path)
    fwd, ctrl, cov = ReadFwdCtrl2D(fwd_path), ReadInvCtrl2D(inv_path), ReadCov2D(cov_path)
    start, prior = ReadModel2D(start_path), ReadModel2D(prior_path)
    _same_grid(start, prior) || error("the prior model must be on the start model's grid")
    size(cov.mask) == size(start.resistivity) ||
        error("the covariance mask is $(size(cov.mask)) cells, the model is $(size(start.resistivity))")
    mesh, ρ0 = Mesh2DFromInputs(start, observed, fwd)
    _, ρref = Mesh2DFromInputs(prior, observed, fwd)
    active = falses(size(ρ0))
    active[mesh.n_air_cells+1:end, :] .= cov.mask .!= 0
    any(active) || error("the covariance mask fixes every cell")

    dir = run_dir === nothing ? _inv2d_run_dir(data_path) : String(run_dir)
    inputs = joinpath(dir, "inputs")
    mkpath(inputs)
    for p in (start_path, data_path, fwd_path, inv_path, cov_path, prior_path)
        cp(p, joinpath(inputs, basename(p)); force = true)
    end
    println("Run directory: ", dir)

    history, reason, converged = nothing, :max_iter, false
    if ctrl.algorithm == :vfsa
        final, predicted = _invert2d_vfsa(dir, mesh, ρ0, observed, ctrl, fwd, cov)
    else
        options, algorithm = _inv2d_from_ctrl(ctrl)
        result = Invert2D(mesh, ρ0, observed; algorithm, options, active_cells = active, reference_resistivity = ρref)
        final, history, reason, converged = result.resistivity, result.history, result.reason, result.converged
        predicted = _inv2d_predicted(result.response, observed, ctrl.mode)
    end
    components = ctrl.mode == :TETM ? ["ZXY", "ZYX"] : ctrl.mode == :TE ? ["ZXY"] : ["ZYX"]
    fit = chi2_rms2d(observed, predicted; components)
    if ctrl.algorithm == :vfsa
        converged = fit.rms <= ctrl.target_rms
        reason = converged ? :target_rms : :max_iter
    end

    WriteModel2D(joinpath(dir, "model.rho"), mesh, final)
    write_data2d(joinpath(dir, "data.pred"), predicted)
    if history !== nothing
        open(joinpath(dir, "History.csv"), "w") do io
            println(io, join(string.(keys(first(history))), ','))
            foreach(h -> println(io, join(values(h), ',')), history)
        end
    end
    open(joinpath(dir, "Summary.txt"), "w") do io
        println(io, "Algorithm: ", uppercase(string(ctrl.algorithm)))
        println(io, "Termination: ", reason)
        println(io, "Converged: ", converged)
        @printf(io, "RMS: %.6f\n", fit.rms)
        println(io, "Real data count: ", fit.count)
        history === nothing || println(io, "Accepted iterations: ", length(history) - 1)
        println(io, "Active cells: ", count(active))
    end
    (; run_dir = dir, algorithm = ctrl.algorithm, ctrl, fwd, mesh, observed, predicted, start = ρ0, prior = ρref,
       final, active, history, rms = fit.rms, reason, converged)
end

# 2D MT nonlinear conjugate-gradient inversion
# Author: @pankajkmishra
# Preconditioned Polak–Ribière+ NLCG direction for the Inv2D driver
# Gradient only (Gᵗ by adjoint solves), G never formed, so it scales to large active models

using LinearAlgebra
using SparseArrays

#---------- configuration ----------

"""
    NLCG2DConfig(; restart=30, precondition=true, epsilon=1e-2)

Nonlinear conjugate-gradient algorithm for `Invert2D`.
- `restart`: reset to preconditioned steepest descent every `restart` iterations
- `precondition`: apply `P = (βR'R + εI)⁻¹` to the gradient, the model-covariance
  preconditioner of MT NLCG codes; smooths updates the same way the regularization does
- `epsilon`: `ε` relative to the mean diagonal of `βR'R`

Each iteration costs one adjoint gradient (one solve per frequency and mode) plus the
line-search forward solves, so NLCG needs more iterations than Gauss–Newton but each
one is far cheaper and memory stays linear in the number of cells.
"""
Base.@kwdef struct NLCG2DConfig <: AbstractInversion2D
    restart::Int = 30
    precondition::Bool = true
    epsilon::Float64 = 1e-2
end

const NLCG2DResult = Inv2DResult{NLCG2DConfig}

inv2d_tag(::NLCG2DConfig) = "nlcg"

function inv2d_validate(alg::NLCG2DConfig)
    alg.restart > 0 || throw(ArgumentError("restart must be positive"))
    isfinite(alg.epsilon) && alg.epsilon > 0 || throw(ArgumentError("epsilon must be finite and positive"))
    nothing
end

#---------- interface ----------

mutable struct NLCG2DWork
    P::Any                      # preconditioner factor, nothing = identity
    g::Vector{Float64}          # this iteration: gradient,
    h::Vector{Float64}          # preconditioned gradient,
    d::Vector{Float64}          # search direction
    g_prev::Vector{Float64}     # last accepted iteration
    h_prev::Vector{Float64}
    d_prev::Vector{Float64}
    scale::Float64              # step length that the last accepted step implied
    k::Int                      # iterations since the last restart
    steepest::Bool              # current d is -h
end

function inv2d_init(alg::NLCG2DConfig, problem, state)
    P = nothing
    beta = problem.options.beta
    if alg.precondition && beta > 0
        A = beta * (problem.R_active' * problem.R_active)
        ε = alg.epsilon * max(mean(diag(A)), eps())
        P = cholesky(Symmetric(A + ε * I))
    end
    e = Float64[]
    NLCG2DWork(P, e, e, e, e, e, e, 0.0, 0, true)
end

inv2d_info(::NLCG2DConfig, work::NLCG2DWork) = (; restart = work.steepest ? 1 : 0)

inv2d_prepare!(alg::NLCG2DConfig, work::NLCG2DWork, problem, state) = inv2d_gradient(problem, state)

function inv2d_direction(alg::NLCG2DConfig, work::NLCG2DWork, problem, state, gradient)
    g = gradient
    h = work.P === nothing ? copy(g) : work.P \ g

    # polak-ribière+; restart on schedule or when the result is not a descent direction
    d = -h
    work.steepest = true
    if work.k > 0 && work.k % alg.restart != 0
        β = max(0.0, dot(g, h - work.h_prev) / dot(work.g_prev, work.h_prev))
        candidate = -h + β * work.d_prev
        if β > 0 && dot(candidate, g) < 0
            d = candidate
            work.steepest = false
        end
    end
    work.g, work.h, work.d = g, h, d

    # initial step length: the first iteration takes max_step, later ones carry over the
    # last step's slope, α = α_prev (g_prev·d_prev)/(g·d), doubled so a step that was cut
    # short can grow back; the driver only caps and backtracks, it never expands
    s = if work.scale > 0 && !isempty(work.d_prev)
        2 * work.scale * dot(work.g_prev, work.d_prev) / dot(g, d)
    else
        problem.options.max_step / max(norm(d, Inf), eps())
    end
    s = isfinite(s) && s > 0 ? s : problem.options.max_step / max(norm(d, Inf), eps())
    s .* d
end

# a failed conjugate direction gets one retry as steepest descent, then the run ends
function inv2d_reject!(alg::NLCG2DConfig, work::NLCG2DWork)
    work.steepest && return false
    work.k = 0
    work.scale = 0.0
    true
end

function inv2d_accept!(alg::NLCG2DConfig, work::NLCG2DWork, problem, old, new)
    # step actually taken along d, after capping, backtracking, and the bounds
    delta = new.m - old.m
    work.scale = max(dot(delta, work.d) / dot(work.d, work.d), 0.0)
    work.g_prev, work.h_prev, work.d_prev = work.g, work.h, work.d
    work.k = work.steepest ? 1 : work.k + 1
    nothing
end

#---------- entry points ----------

"""
    NLCG2D(mesh, initial_resistivity, observed; config=NLCG2DConfig(), options=Inv2DOptions(), ...)
    NLCG2D(model_path, data_path; output_dir=nothing, config=..., options=..., ...)

NLCG shorthand for `Invert2D(...; algorithm=config)`. The file method writes
`model_nlcg.rho`, `data_nlcg.dat`, `history_nlcg.csv`, and `summary_nlcg.txt`.
"""
NLCG2D(mesh::MT2DMesh, initial_resistivity::AbstractMatrix{<:Real}, observed::DataFile2D;
       config::NLCG2DConfig = NLCG2DConfig(), kwargs...) =
    Invert2D(mesh, initial_resistivity, observed; algorithm = config, kwargs...)

NLCG2D(model_path::AbstractString, data_path::AbstractString;
       config::NLCG2DConfig = NLCG2DConfig(), kwargs...) =
    Invert2D(model_path, data_path; algorithm = config, kwargs...)

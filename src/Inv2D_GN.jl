# 2D MT Gauss–Newton inversion
# Author: @pankajkmishra
# Damped (Levenberg–Marquardt) Gauss–Newton direction for the Inv2D driver
# Explicit Fréchet derivative G from implicit differentiation, dense normal equations

using LinearAlgebra

#---------- configuration ----------

"""
    GaussNewton2DConfig(; damping=1e-2, max_damping_trials=6)

Gauss–Newton algorithm for `Invert2D`.
- `damping`: initial Levenberg–Marquardt damping, relative to the column scaling of
  `GᵗG + βRᵗR`; divided by 3 after an accepted step, multiplied by 10 after a failed search
- `max_damping_trials`: failed searches per iteration before `:line_search_failed`

The explicit G costs one tangent linear solve per active cell and the step solves a dense
`n_active × n_active` system, so it suits small to medium active models.
"""
Base.@kwdef struct GaussNewton2DConfig <: AbstractInversion2D
    damping::Float64 = 1e-2
    max_damping_trials::Int = 6
end

const GaussNewton2DResult = Inv2DResult{GaussNewton2DConfig}

inv2d_tag(::GaussNewton2DConfig) = "gn"

function inv2d_validate(alg::GaussNewton2DConfig)
    isfinite(alg.damping) && alg.damping > 0 || throw(ArgumentError("damping must be finite and positive"))
    alg.max_damping_trials > 0 || throw(ArgumentError("max_damping_trials must be positive"))
    nothing
end

#---------- interface ----------

mutable struct GaussNewton2DWork
    damping::Float64
    trials::Int
    H::Matrix{Float64}          # GᵗG + βRᵗR
    scaling::Vector{Float64}    # sqrt of its diagonal
end

inv2d_init(alg::GaussNewton2DConfig, problem, state) =
    GaussNewton2DWork(alg.damping, 0, zeros(0, 0), Float64[])

inv2d_info(::GaussNewton2DConfig, work::GaussNewton2DWork) = (; damping = work.damping)

function inv2d_prepare!(alg::GaussNewton2DConfig, work::GaussNewton2DWork, problem, state)
    beta = problem.options.beta
    G = inv2d_frechet(problem, state)
    work.H = G' * G + beta * Matrix(problem.R_active' * problem.R_active)
    work.scaling = sqrt.(max.(diag(work.H), 1e-12))
    work.trials = 0
    inv2d_gradient(problem, state, G)
end

# (GᵗG + βRᵗR + λ S²) d = -g, S = column scaling; same step as the damped augmented least squares
function inv2d_direction(alg::GaussNewton2DConfig, work::GaussNewton2DWork, problem, state, gradient)
    A = copy(work.H)
    for i in axes(A, 1)
        A[i, i] += work.damping * work.scaling[i]^2
    end
    -(cholesky!(Symmetric(A)) \ gradient)
end

function inv2d_reject!(alg::GaussNewton2DConfig, work::GaussNewton2DWork)
    work.damping *= 10
    work.trials += 1
    work.trials < alg.max_damping_trials
end

function inv2d_accept!(alg::GaussNewton2DConfig, work::GaussNewton2DWork, problem, old, new)
    work.damping = max(work.damping / 3, 1e-12)
    nothing
end

#---------- entry points ----------

"""
    GaussNewton2D(mesh, initial_resistivity, observed; config=GaussNewton2DConfig(),
                  options=Inv2DOptions(), active_cells=nothing,
                  reference_resistivity=initial_resistivity)

Gauss–Newton shorthand for `Invert2D(...; algorithm=config)`.
"""
GaussNewton2D(mesh::MT2DMesh, initial_resistivity::AbstractMatrix{<:Real}, observed::DataFile2D;
              config::GaussNewton2DConfig = GaussNewton2DConfig(), kwargs...) =
    Invert2D(mesh, initial_resistivity, observed; algorithm = config, kwargs...)

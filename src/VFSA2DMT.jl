# This script runs a full 2D VFSA inversion workflow with RBF-controlled models and shared 2D file I/O.

if !isdefined(@__MODULE__, :plot_mt2d_model)
    include(joinpath(@__DIR__, "MTGeophysics2D.jl"))
end

using CairoMakie
using Dates
using Printf
using Random
using Statistics

const air_resistivity = 1e9

Base.@kwdef struct VFSA2DMTConfig
    n_chains::Int = 2
    n_ctrl::Int = 400
    max_iter::Int = 3000
    # n_trials = 1 is classic vfsa: one proposal, one metropolis test per iteration
    n_trials::Int = 1
    log_bounds::Tuple{Float64, Float64} = (0.0, 5.0)
    frac_update_controls::Float64 = 1.0
    step_scale::Float64 = 0.11
    # one schedule drives both proposal width and acceptance (same as the 3D driver),
    # cooling from temp_kappa to cool_ratio*temp_kappa at max_iter; pick temp_kappa
    # on the scale of the typical uphill dE = (rms² - rms₀²)/rms₀²
    temp_kappa::Float64 = 1.0
    cool_ratio::Float64 = 1e-3
    target_rms::Float64 = 1.0
    seed::Int = 20260308
    pad_tolerance::Float64 = 0.20
    padding_decay_length::Float64 = 8.0
    rbf_sigma_scale_y::Float64 = 2.0
    rbf_sigma_scale_z::Float64 = 2.5
    trunc_sigmas::Float64 = 3.0
    perturb_depth_m::Float64 = Inf
    keep_models::Bool = true
    snapshot_interval::Int = 10
    output_root::String = dirname(@__DIR__)
end

Base.@kwdef struct VFSA2DMTParams
    script_path::String
    start_model_path::String
    data_path::String
    model_path::String = ""
    run_name::String = "run_VFSA2DMT"
    config::VFSA2DMTConfig = VFSA2DMTConfig()
end

Base.@kwdef struct MT2DIterationRecord
    chain::Int
    iteration::Int
    accepted_trial::Int
    temperature::Float64
    proposal_chi2::Float64
    proposal_rms::Float64
    current_chi2::Float64
    current_rms::Float64
    best_chi2::Float64
    best_rms::Float64
    delta_chi2::Float64
    accepted::Bool
end

Base.@kwdef struct MT2DChainResult
    chain_id::Int
    chain_dir::String
    best_model_path::String
    best_data_path::String
    final_model_path::String
    final_data_path::String
    best_chi2::Float64
    best_rms::Float64
    best_resistivity::Matrix{Float64}
    best_response::MT2DResponse
    best_data::DataFile2D
    final_chi2::Float64
    final_rms::Float64
    final_resistivity::Matrix{Float64}
    final_response::MT2DResponse
    final_data::DataFile2D
    iterations::Vector{MT2DIterationRecord}
    control_y_m::Vector{Float64}
    control_z_m::Vector{Float64}
end

Base.@kwdef struct MT2DEnsembleSummary
    chain_count::Int
    mean_resistivity::Matrix{Float64}
    median_resistivity::Matrix{Float64}
    std_log10_resistivity::Matrix{Float64}
    p05_log10_resistivity::Matrix{Float64}
    p95_log10_resistivity::Matrix{Float64}
    mean_response::MT2DResponse
    mean_data::DataFile2D
    mean_fit::FitSummary2D
end

struct MT2DRBFMap
    ny::Int
    nz::Int
    control_y::Vector{Int}
    control_z::Vector{Int}
    ptr::Vector{Int}
    nbrs::Vector{Int}
    weights::Vector{Float64}
    y_centers::Vector{Float64}
    z_centers::Vector{Float64}
end

core_cell_indices(cell_sizes::AbstractVector{<:Real}; tol::Real = 0.20) = begin
    minimum_size = minimum(Float64.(cell_sizes))
    threshold = minimum_size * (1 + Float64(tol))
    indices = findall(size -> Float64(size) <= threshold + 1e-9, cell_sizes)
    isempty(indices) && error("no core cells were detected")
    first(indices):last(indices)
end

function _ensure_file_in_dir(path::AbstractString, target_dir::AbstractString; target_name::AbstractString = basename(path))
    mkpath(target_dir)
    target_path = joinpath(target_dir, target_name)
    abspath(path) == abspath(target_path) || cp(path, target_path; force = true)
    target_path
end

_vfsa2d_timestamp_slug(now_value::DateTime = now()) = Dates.format(now_value, "yyyymmdd_HHMMSS")

_vfsa2d_timestamp_label(now_value::DateTime = now()) = Dates.format(now_value, "yyyy-mm-dd HH:MM:SS")

function _create_mt2d_run_dir(
    root_dir::AbstractString;
    run_name::AbstractString = "VFSA2DMT",
    now_value::DateTime = now(),
)
    timestamp_slug = _vfsa2d_timestamp_slug(now_value)
    timestamp_label = _vfsa2d_timestamp_label(now_value)
    run_dir = joinpath(root_dir, "$(run_name)_$(timestamp_slug)")
    mkpath(run_dir)
    (
        run_dir = run_dir,
        workflow = String(run_name),
        timestamp_slug = timestamp_slug,
        timestamp_label = timestamp_label,
    )
end

function _maybe_copy_mt2d_true_model(observed_data_path::AbstractString, run_dir::AbstractString)
    stem, _ = splitext(basename(observed_data_path))
    candidate = joinpath(dirname(observed_data_path), "$(stem).true")
    if !isfile(candidate)
        candidate = joinpath(dirname(observed_data_path), "$(stem).rho")
    end
    isfile(candidate) || return ""
    _ensure_file_in_dir(candidate, run_dir; target_name = _vfsa2d_true_model_output_filename())
end

function _vfsa_output_dirs(run_dir::AbstractString)
    mkpath(run_dir)
    (
        logs = run_dir,
        plots = run_dir,
        models = run_dir,
        data = run_dir,
    )
end

function _vfsa2d_chain_slug(chain_id::Integer)
    @sprintf("chain_%02d", chain_id)
end

_vfsa2d_trials_log_filename() = "0vfsa2DMT.details"

_vfsa2d_iterations_log_filename() = "0vfsa2DMT.log"

_vfsa2d_observed_output_filename() = "data.obs"

_vfsa2d_true_model_output_filename() = "model.true"

_vfsa2d_mean_model_filename() = "model.mean"

_vfsa2d_median_model_filename() = "model.median"

_vfsa2d_std_model_filename() = "model.std"

_vfsa2d_convergence_plot_filename() = "plot_convergence.png"

_vfsa2d_mean_model_plot_filename() = "plot_model_mean.png"

_vfsa2d_median_model_plot_filename() = "plot_model_median.png"

_vfsa2d_observed_maps_plot_filename() = "plot_data_obs_maps.png"

_vfsa2d_best_maps_plot_filename() = "plot_data_best_maps.png"

_vfsa2d_observed_curves_plot_filename() = "plot_data_obs_curves.png"

_vfsa2d_best_curves_plot_filename() = "plot_data_best_curves.png"

function _vfsa2d_trial_model_filename(iteration::Integer, trial::Integer)
    @sprintf("itr_%05d_trial_%02d.rho", iteration, trial)
end

function _vfsa2d_best_model_filename(chain_id::Integer)
    "model.c$(chain_id)best"
end

function _vfsa2d_best_prediction_filename(chain_id::Integer)
    "data.c$(chain_id)best"
end

function _padding_background_log10(log10_resistivity::AbstractMatrix{<:Real}, mesh::MT2DMesh, core_y::UnitRange{Int})
    values = Float64[]
    for iy in eachindex(mesh.y_cell_sizes)
        iy in core_y && continue
        for iz in (mesh.n_air_cells + 1):size(log10_resistivity, 1)
            push!(values, Float64(log10_resistivity[iz, iy]))
        end
    end
    isempty(values) && return median(vec(Float64.(log10_resistivity[(mesh.n_air_cells + 1):end, core_y])))
    median(values)
end

function _smooth_padding_decay_y!(
    full_log10::AbstractMatrix{<:Real},
    mesh::MT2DMesh,
    core_y::UnitRange{Int},
    background_log10::Real,
    decay_length::Real,
)
    y_centers = mt2d_y_centers(mesh)
    left_boundary = first(core_y)
    right_boundary = last(core_y)
    decay_scale = Float64(decay_length) * median(mesh.y_cell_sizes[core_y])

    for iy in eachindex(mesh.y_cell_sizes)
        iy in core_y && continue
        boundary_index = iy < left_boundary ? left_boundary : right_boundary
        distance = abs(y_centers[iy] - y_centers[boundary_index])
        weight = exp(-distance / max(decay_scale, eps(Float64)))
        for iz in (mesh.n_air_cells + 1):length(mesh.z_cell_sizes)
            boundary_value = full_log10[iz, boundary_index]
            full_log10[iz, iy] = boundary_value * weight + Float64(background_log10) * (1 - weight)
        end
    end
    full_log10[1:mesh.n_air_cells, :] .= log10(air_resistivity)
    full_log10
end

function _perturbable_ground_indices(mesh::MT2DMesh, perturb_depth_m::Real)
    ground_top = mesh.n_air_cells + 1
    ground_bottom = length(mesh.z_cell_sizes)
    ground_top > ground_bottom && error("mesh does not contain ground cells")
    isfinite(perturb_depth_m) || return ground_top:ground_bottom

    limited_depth = max(Float64(perturb_depth_m), 0.0)
    ground_cell_bottoms = mesh.z_nodes[(ground_top + 1):(ground_bottom + 1)]
    allowed = findall(z_bottom -> z_bottom <= limited_depth + 1e-9, ground_cell_bottoms)
    deepest_offset = isempty(allowed) ? 1 : last(allowed)
    ground_top:(ground_top + deepest_offset - 1)
end

function _smooth_padding_decay_z!(
    full_log10::AbstractMatrix{<:Real},
    mesh::MT2DMesh,
    core_z::UnitRange{Int},
    background_log10::Real,
    decay_length::Real,
)
    deepest_boundary = last(core_z)
    deepest_boundary == length(mesh.z_cell_sizes) && return full_log10

    z_centers = mt2d_z_centers(mesh)
    decay_scale = Float64(decay_length) * median(mesh.z_cell_sizes[core_z])
    for iz in (deepest_boundary + 1):length(mesh.z_cell_sizes)
        distance = abs(z_centers[iz] - z_centers[deepest_boundary])
        weight = exp(-distance / max(decay_scale, eps(Float64)))
        for iy in eachindex(mesh.y_cell_sizes)
            boundary_value = full_log10[deepest_boundary, iy]
            full_log10[iz, iy] = boundary_value * weight + Float64(background_log10) * (1 - weight)
        end
    end
    full_log10[1:mesh.n_air_cells, :] .= log10(air_resistivity)
    full_log10
end

function build_mt2d_rbf_map(
    mesh::MT2DMesh,
    core_y::UnitRange{Int},
    core_z::UnitRange{Int},
    n_ctrl::Integer,
    rng::AbstractRNG;
    sigma_scale_y::Real,
    sigma_scale_z::Real,
    trunc_sigmas::Real,
)
    y_core = mt2d_y_centers(mesh)[core_y]
    z_core = mt2d_z_centers(mesh)[core_z]
    n_y = length(y_core)
    n_z = length(z_core)
    total_cells = n_y * n_z
    n_selected = min(Int(n_ctrl), total_cells)
    chosen = randperm(rng, total_cells)[1:n_selected]

    control_y = Vector{Int}(undef, n_selected)
    control_z = similar(control_y)
    for (control_index, flat_index) in enumerate(chosen)
        iz = Int(ceil(flat_index / n_y))
        iy = flat_index - (iz - 1) * n_y
        control_y[control_index] = iy
        control_z[control_index] = iz
    end

    sigma_y = Float64(sigma_scale_y) * median(mesh.y_cell_sizes[core_y])
    sigma_z = Float64(sigma_scale_z) * median(mesh.z_cell_sizes[core_z])
    truncation = Float64(trunc_sigmas)

    ptr = Vector{Int}(undef, total_cells + 1)
    neighbors = Int[]
    weights = Float64[]
    row = 1
    ptr[1] = 1

    for iz in 1:n_z, iy in 1:n_y
        local_neighbors = Int[]
        local_weights = Float64[]
        for control_index in eachindex(control_y)
            Δy = y_core[iy] - y_core[control_y[control_index]]
            Δz = z_core[iz] - z_core[control_z[control_index]]
            r2 = (Δy / sigma_y)^2 + (Δz / sigma_z)^2
            if r2 <= truncation^2
                push!(local_neighbors, control_index)
                push!(local_weights, exp(-0.5 * r2))
            end
        end

        if isempty(local_neighbors)
            nearest = argmin([
                (y_core[iy] - y_core[control_y[control_index]])^2 + (z_core[iz] - z_core[control_z[control_index]])^2 for
                control_index in eachindex(control_y)
            ])
            push!(local_neighbors, nearest)
            push!(local_weights, 1.0)
        end

        weight_sum = sum(local_weights)
        for local_index in eachindex(local_weights)
            local_weights[local_index] /= weight_sum
        end

        append!(neighbors, local_neighbors)
        append!(weights, local_weights)
        ptr[row + 1] = length(neighbors) + 1
        row += 1
    end

    MT2DRBFMap(
        n_y,
        n_z,
        control_y,
        control_z,
        ptr,
        neighbors,
        weights,
        y_core,
        z_core,
    )
end

function apply_mt2d_rbf_map!(delta_values::AbstractMatrix{<:Real}, map::MT2DRBFMap, delta_params::AbstractVector{<:Real})
    size(delta_values) == (map.nz, map.ny) || error("delta_values must be size ($(map.nz), $(map.ny))")
    row = 1
    for iz in 1:map.nz, iy in 1:map.ny
        value = 0.0
        start_index = map.ptr[row]
        stop_index = map.ptr[row + 1] - 1
        for pointer_index in start_index:stop_index
            control_index = map.nbrs[pointer_index]
            value += map.weights[pointer_index] * delta_params[control_index]
        end
        delta_values[iz, iy] = value
        row += 1
    end
    delta_values
end

@inline function _vfsa_y(u::Float64, temperature::Float64)
    sign = u >= 0.5 ? 1.0 : -1.0
    sign * temperature * ((1 + 1 / temperature) ^ abs(2u - 1.0) - 1.0)
end

function _propose_controls!(
    delta_params::Vector{Float64},
    temperature::Float64,
    lower_bound::Float64,
    upper_bound::Float64,
    base_control_values::Vector{Float64},
    n_selected::Int,
    rng::AbstractRNG;
    step_scale::Float64,
)
    n_parameters = length(delta_params)
    chosen = randperm(rng, n_parameters)[1:n_selected]
    proposal_span = (upper_bound - lower_bound) * step_scale
    for parameter_index in chosen
        perturbation = _vfsa_y(rand(rng), temperature) * proposal_span
        trial_delta = delta_params[parameter_index] + perturbation
        absolute_value = base_control_values[parameter_index] + trial_delta
        if absolute_value < lower_bound
            trial_delta = lower_bound - base_control_values[parameter_index]
        elseif absolute_value > upper_bound
            trial_delta = upper_bound - base_control_values[parameter_index]
        end
        delta_params[parameter_index] = trial_delta
    end
    chosen
end

function _print_vfsa_progress_bar(chain_id::Int, iteration::Int, total_iterations::Int, current_rms::Real, best_rms::Real)
    total_iterations > 0 || return nothing
    bar_width = 32
    fraction_done = clamp(Float64(iteration) / total_iterations, 0.0, 1.0)
    filled_width = min(bar_width, round(Int, fraction_done * bar_width))
    filled_width = max(filled_width, iteration > 0 ? 1 : 0)
    filled = repeat("=", max(filled_width - 1, 0))
    head = fraction_done < 1 ? ">" : "="
    empty = repeat(".", max(bar_width - filled_width, 0))
    bar = filled_width == 0 ? empty : string(filled, head, empty)
    message = @sprintf(
        "\rVFSA2DMT chain %d [%s] %d/%d  curr_rms=%.4f  best_rms=%.4f",
        chain_id,
        bar,
        iteration,
        total_iterations,
        Float64(current_rms),
        Float64(best_rms),
    )
    print(stdout, message)
    flush(stdout)
    if iteration == total_iterations
        println()
    end
    nothing
end

function _reconstruct_trial_model(
    start_log10::AbstractMatrix{<:Real},
    start_core_log10::AbstractMatrix{<:Real},
    delta_buffer::AbstractMatrix{<:Real},
    map::MT2DRBFMap,
    delta_params::AbstractVector{<:Real},
    mesh::MT2DMesh,
    core_y::UnitRange{Int},
    core_z::UnitRange{Int},
    background_log10::Real,
    config::VFSA2DMTConfig,
)
    apply_mt2d_rbf_map!(delta_buffer, map, delta_params)
    trial_core = clamp.(start_core_log10 .+ delta_buffer, config.log_bounds[1], config.log_bounds[2])
    full_log10 = copy(Float64.(start_log10))
    full_log10[core_z, core_y] .= trial_core
    _smooth_padding_decay_y!(full_log10, mesh, core_y, background_log10, config.padding_decay_length)
    _smooth_padding_decay_z!(full_log10, mesh, core_z, background_log10, config.padding_decay_length)
    full_log10[1:mesh.n_air_cells, :] .= log10(air_resistivity)
    resistivity = 10 .^ full_log10
    resistivity[1:mesh.n_air_cells, :] .= air_resistivity
    (log10_resistivity = full_log10, resistivity = resistivity)
end

function _predicted_data2d(response::MT2DResponse, observed_data::DataFile2D)
    data_from_response2d(
        response;
        z_xy_error = observed_data.z_xy_error,
        z_yx_error = observed_data.z_yx_error,
        z_xx_error = observed_data.z_xx_error,
        z_yy_error = observed_data.z_yy_error,
        site_names = observed_data.site_names,
        x_positions = observed_data.x_positions,
        z_positions = observed_data.z_positions,
        title = "Predicted 2D response",
    )
end

function _write_trials_header(path::AbstractString, config::VFSA2DMTConfig, timestamp::AbstractString)
    open(path, "w") do io
        println(io, "# VFSA 2D detailed trials - ", timestamp)
        println(
            io,
            "# chains=", config.n_chains,
            "  n_ctrl=", config.n_ctrl,
            "  frac_update_controls=", config.frac_update_controls,
            "  temp_kappa=", config.temp_kappa,
            "  cool_ratio=", config.cool_ratio,
            "  target_rms=", config.target_rms,
            "  n_trials=", config.n_trials,
        )
        println(io, repeat("-", 270))
        # one row per trial; only the winning trial carries a metropolis draw,
        # the rest print an em dash (same scheme as the 3D driver)
        @printf(
            io,
            "%6s %8s %8s %12s %14s %11s %14s %12s %12s %10s %10s %10s %10s %12s %11s %5s %7s %s %s\n",
            "Chain",
            "Iter",
            "Trial",
            "Temp",
            "Chi2Before",
            "RMSBefore",
            "Chi2Prop",
            "dChi2",
            "RMSProp",
            "dRMS",
            "dE_rms2",
            "Pacc",
            "Uacc",
            "Chi2Best",
            "RMSBest",
            "Acc",
            "UpdCtrl",
            "Model",
            "Data",
        )
        println(io, repeat("-", 270))
    end
end

function _append_trial_row(
    path::AbstractString;
    chain::Int,
    iteration::Int,
    trial::Int,
    T::Float64,
    chi2_before::Float64,
    rms_before::Float64,
    chi2_prop::Float64,
    dchi2::Float64,
    rms_prop::Float64,
    drms::Float64,
    dE::Float64,
    p_acc::Float64,
    u_acc::Float64,
    chi2_best::Float64,
    rms_best::Float64,
    acc::Int,
    updated_controls::Int,
    model_rel::String,
    data_rel::String,
)
    pacc_str = isnan(p_acc) ? lpad("—", 10) : @sprintf("%10.5f", p_acc)
    uacc_str = isnan(u_acc) ? lpad("—", 10) : @sprintf("%10.5f", u_acc)
    open(path, "a") do io
        @printf(
            io,
            "%6d %8d %8d %12.4g %14.4f %11.5f %14.4f %12.4f %12.5f %10.5f %10.5f %s %s %12.4f %11.5f %5d %7d %s %s\n",
            chain,
            iteration,
            trial,
            T,
            chi2_before,
            rms_before,
            chi2_prop,
            dchi2,
            rms_prop,
            drms,
            dE,
            pacc_str,
            uacc_str,
            chi2_best,
            rms_best,
            acc,
            updated_controls,
            isempty(model_rel) ? "-" : model_rel,
            isempty(data_rel) ? "-" : data_rel,
        )
    end
end

function _write_iteration_header(path::AbstractString, config::VFSA2DMTConfig, timestamp::AbstractString)
    open(path, "w") do io
        println(io, "# VFSA 2D iteration best - ", timestamp)
        println(
            io,
            "# chains=", config.n_chains,
            "  n_ctrl=", config.n_ctrl,
            "  frac_update_controls=", config.frac_update_controls,
            "  temp_kappa=", config.temp_kappa,
            "  cool_ratio=", config.cool_ratio,
            "  target_rms=", config.target_rms,
            "  n_trials=", config.n_trials,
        )
        println(io, repeat("-", 200))
        # current accepted state after the iteration; Nacc = 1 when the winning
        # trial passed its metropolis test, 0 otherwise
        @printf(
            io,
            "%6s %8s %12s %14s %12s %12s %11s %12s %11s %6s %s\n",
            "Chain",
            "Iter",
            "Temp",
            "Chi2Curr",
            "dChi2",
            "RMSCurr",
            "CurrRMS",
            "Chi2Best",
            "RMSBest",
            "Nacc",
            "Model",
        )
        println(io, repeat("-", 200))
    end
end

function _append_iteration_row(
    path::AbstractString;
    chain::Int,
    iteration::Int,
    T::Float64,
    chi2_curr::Float64,
    dchi2::Float64,
    rms_curr::Float64,
    current_rms::Float64,
    chi2_best::Float64,
    rms_best::Float64,
    nacc::Int,
    model_rel::String,
)
    open(path, "a") do io
        @printf(
            io,
            "%6d %8d %12.4g %14.4f %12.4f %12.5f %11.5f %12.4f %11.5f %6d %s\n",
            chain,
            iteration,
            T,
            chi2_curr,
            dchi2,
            rms_curr,
            current_rms,
            chi2_best,
            rms_best,
            nacc,
            isempty(model_rel) ? "-" : model_rel,
        )
    end
end

function _summarize_mt2d_resistivity_ensemble(
    resistivities::AbstractVector{<:AbstractMatrix{<:Real}},
    mesh::MT2DMesh,
    observed_data::DataFile2D,
)
    isempty(resistivities) && error("at least one model is required to build an ensemble summary")

    n_z, n_y = size(first(resistivities))
    n_chains = length(resistivities)
    log10_stack = Array{Float64}(undef, n_z, n_y, n_chains)

    for (chain_index, resistivity) in enumerate(resistivities)
        size(resistivity) == (n_z, n_y) || error("all ensemble models must share the same shape")
        log10_stack[:, :, chain_index] .= log10.(Float64.(resistivity))
    end

    mean_log10 = dropdims(mean(log10_stack; dims = 3), dims = 3)
    std_log10 = dropdims(std(log10_stack; dims = 3, corrected = false), dims = 3)
    median_log10 = similar(mean_log10)
    p05_log10 = similar(mean_log10)
    p95_log10 = similar(mean_log10)
    samples = Vector{Float64}(undef, n_chains)

    for iz in 1:n_z, iy in 1:n_y
        for chain_index in 1:n_chains
            samples[chain_index] = log10_stack[iz, iy, chain_index]
        end
        median_log10[iz, iy] = median(samples)
        p05_log10[iz, iy] = quantile(samples, 0.05)
        p95_log10[iz, iy] = quantile(samples, 0.95)
    end

    mean_resistivity = 10 .^ mean_log10
    median_resistivity = 10 .^ median_log10
    mean_response = run_mt2d_forward(mesh, mean_resistivity)
    mean_data = _predicted_data2d(mean_response, observed_data)
    mean_fit = chi2_rms2d(observed_data, mean_data)

    MT2DEnsembleSummary(
        chain_count = n_chains,
        mean_resistivity = mean_resistivity,
        median_resistivity = median_resistivity,
        std_log10_resistivity = std_log10,
        p05_log10_resistivity = p05_log10,
        p95_log10_resistivity = p95_log10,
        mean_response = mean_response,
        mean_data = mean_data,
        mean_fit = mean_fit,
    )
end

function _summarize_best_chain_models(
    chains::AbstractVector{MT2DChainResult},
    mesh::MT2DMesh,
    observed_data::DataFile2D,
)
    _summarize_mt2d_resistivity_ensemble(getfield.(chains, :best_resistivity), mesh, observed_data)
end

function _write_mt2d_uncertainty_table(
    path::AbstractString,
    mesh::MT2DMesh,
    ensemble::MT2DEnsembleSummary,
)
    mkpath(dirname(path))
    y_centers = mt2d_y_centers(mesh)
    z_centers = mt2d_z_centers(mesh)
    open(path, "w") do io
        println(
            io,
            "z_index,y_index,z_m,y_m,std_log10_rho,p05_log10_rho,p95_log10_rho,mean_rho_ohm_m,median_rho_ohm_m",
        )
        for iz in axes(ensemble.mean_resistivity, 1), iy in axes(ensemble.mean_resistivity, 2)
            @printf(
                io,
                "%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                iz,
                iy,
                z_centers[iz],
                y_centers[iy],
                ensemble.std_log10_resistivity[iz, iy],
                ensemble.p05_log10_resistivity[iz, iy],
                ensemble.p95_log10_resistivity[iz, iy],
                ensemble.mean_resistivity[iz, iy],
                ensemble.median_resistivity[iz, iy],
            )
        end
    end
    path
end

function AnalyseEnsemble2D(
    run_dir::AbstractString;
    model_paths::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    observed_data_path::Union{Nothing, AbstractString} = nothing,
    copy_script::Bool = false,
)
    output_dirs = _vfsa_output_dirs(run_dir)
    models_dir = output_dirs.models
    data_dir = output_dirs.data
    resolved_model_paths = model_paths === nothing ?
        sort(filter(path -> occursin(r"^model\.c\d+best$", basename(path)), readdir(models_dir; join = true))) :
        String.(model_paths)
    isempty(resolved_model_paths) && error("no best chain models were found in $models_dir")

    observed_path = if observed_data_path !== nothing
        String(observed_data_path)
    else
        candidates = sort(filter(path -> endswith(lowercase(path), ".obs"), readdir(data_dir; join = true)))
        isempty(candidates) && error("no observed-data file was found in $data_dir")
        first(candidates)
    end

    observed_data = load_data2d(observed_path)
    models = load_model2d.(resolved_model_paths)
    mesh = build_mesh_from_model2d(
        first(models);
        frequencies = observed_data.frequencies,
        receiver_positions = observed_data.receivers,
    )
    ensemble = _summarize_mt2d_resistivity_ensemble(getfield.(models, :resistivity), mesh, observed_data)

    mean_model_path = write_model2d(
        joinpath(models_dir, _vfsa2d_mean_model_filename()),
        mesh,
        ensemble.mean_resistivity;
        title = "Best-chain ensemble mean model",
    )
    median_model_path = write_model2d(
        joinpath(models_dir, _vfsa2d_median_model_filename()),
        mesh,
        ensemble.median_resistivity;
        title = "Best-chain ensemble median model",
    )
    uncertainty_table_path = _write_mt2d_uncertainty_table(
        joinpath(models_dir, _vfsa2d_std_model_filename()),
        mesh,
        ensemble,
    )
    summary_path = joinpath(output_dirs.logs, "RunStatistics2D.md")
    open(summary_path, "w") do io
        println(io, "# 2D chain-best ensemble summary")
        println(io)
        println(io, "- chain_count: ", ensemble.chain_count)
        println(io, "- observed_data: ", basename(observed_path))
        println(io, "- mean_model: ", basename(mean_model_path))
        println(io, "- median_model: ", basename(median_model_path))
        println(io, "- uncertainty_table: ", basename(uncertainty_table_path))
        println(io, "- ensemble_mean_chi_square: ", round(ensemble.mean_fit.chi2, digits = 4))
        println(io, "- ensemble_mean_rms: ", round(ensemble.mean_fit.rms, digits = 4))
    end

    if copy_script
        cp(joinpath(@__DIR__, "VFSA2DMT.jl"), joinpath(run_dir, "VFSA2DMT.jl"); force = true)
    end

    (
        chain_count = ensemble.chain_count,
        model_paths = resolved_model_paths,
        observed_data_path = observed_path,
        mean_model_path = mean_model_path,
        median_model_path = median_model_path,
        mean_data_path = "",
        uncertainty_table_path = uncertainty_table_path,
        summary_path = summary_path,
        ensemble = ensemble,
    )
end

function _plot_mt2d_vfsa_convergence(path::AbstractString, chains::AbstractVector{MT2DChainResult})
    CairoMakie.activate!()
    mkpath(dirname(path))

    figure = Figure(size = (1400, 950))
    ax1 = Axis(figure[1, 1], xlabel = "VFSA iteration", ylabel = "Best χ²", yscale = log10, title = "Best objective")
    ax2 = Axis(figure[1, 2], xlabel = "VFSA iteration", ylabel = "Current RMS", title = "Current fit")
    ax3 = Axis(figure[2, 1], xlabel = "VFSA iteration", ylabel = "Temperature", yscale = log10, title = "Annealing schedule")
    ax4 = Axis(figure[2, 2], xlabel = "VFSA iteration", ylabel = "Rolling acceptance", title = "Acceptance")

    offset = 0
    plotted = false
    for chain in chains
        isempty(chain.iterations) && continue
        x = offset .+ getfield.(chain.iterations, :iteration)
        best_chi2 = getfield.(chain.iterations, :best_chi2)
        current_rms = getfield.(chain.iterations, :current_rms)
        temperature = getfield.(chain.iterations, :temperature)
        accepted = Float64.(getfield.(chain.iterations, :accepted))
        rolling = [mean(accepted[max(1, index - 14):index]) for index in eachindex(accepted)]
        label = "Chain $(chain.chain_id)"

        lines!(ax1, x, best_chi2, linewidth = 3, label = label)
        lines!(ax2, x, current_rms, linewidth = 3, label = label)
        lines!(ax3, x, temperature, linewidth = 3, label = label)
        lines!(ax4, x, rolling, linewidth = 3, label = label)
        offset += length(chain.iterations)
        plotted = true
    end

    if plotted
        axislegend(ax1, position = :rt)
    else
        text!(ax1, 0.5, 0.5, text = "No VFSA iterations", space = :relative, align = (:center, :center))
    end
    save(path, figure)
    path
end

function _run_mt2d_chain(
    chain_id::Int,
    mesh::MT2DMesh,
    start_resistivity::AbstractMatrix{<:Real},
    observed_data::DataFile2D,
    output_dirs,
    config::VFSA2DMTConfig,
)
    rng = MersenneTwister(config.seed + 1000 * (chain_id - 1))
    chain_slug = _vfsa2d_chain_slug(chain_id)
    chain_dir = joinpath(output_dirs.models, chain_slug)
    mkpath(chain_dir)

    trials_log = joinpath(chain_dir, _vfsa2d_trials_log_filename())
    iteration_log = joinpath(chain_dir, _vfsa2d_iterations_log_filename())
    _write_trials_header(trials_log, config, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    _write_iteration_header(iteration_log, config, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))

    core_y = core_cell_indices(mesh.y_cell_sizes; tol = config.pad_tolerance)
    core_z = _perturbable_ground_indices(mesh, config.perturb_depth_m)

    # Keep the prior in the sparse-control subspace. A dense random cell-by-cell
    # start model defeats the control-point parameterisation used by the legacy 2D
    # workflow and the 3D driver.
    start_log10 = log10.(Float64.(start_resistivity))
    start_core_log10 = Array(start_log10[core_z, core_y])
    background_log10 = _padding_background_log10(start_log10, mesh, core_y)

    rbf_map = build_mt2d_rbf_map(
        mesh,
        core_y,
        core_z,
        config.n_ctrl,
        rng;
        sigma_scale_y = config.rbf_sigma_scale_y,
        sigma_scale_z = config.rbf_sigma_scale_z,
        trunc_sigmas = config.trunc_sigmas,
    )
    control_y_m = rbf_map.y_centers[rbf_map.control_y]
    control_z_m = rbf_map.z_centers[rbf_map.control_z]

    n_controls = length(rbf_map.control_y)
    base_control_values = Vector{Float64}(undef, n_controls)
    for control_index in 1:n_controls
        base_control_values[control_index] = start_core_log10[rbf_map.control_z[control_index], rbf_map.control_y[control_index]]
    end

    current_delta = zeros(Float64, n_controls)
    delta_buffer = zeros(Float64, rbf_map.nz, rbf_map.ny)
    n_update_controls = max(1, round(Int, config.frac_update_controls * n_controls))

    start_state_model = _reconstruct_trial_model(
        start_log10,
        start_core_log10,
        delta_buffer,
        rbf_map,
        current_delta,
        mesh,
        core_y,
        core_z,
        background_log10,
        config,
    )
    current_resistivity = start_state_model.resistivity
    current_response = run_mt2d_forward(mesh, current_resistivity)
    current_data = _predicted_data2d(current_response, observed_data)
    current_fit = chi2_rms2d(observed_data, current_data)

    best_resistivity = copy(current_resistivity)
    best_response = current_response
    best_data = current_data
    best_fit = current_fit
    best_model_path = joinpath(output_dirs.models, _vfsa2d_best_model_filename(chain_id))
    best_data_path = joinpath(output_dirs.data, _vfsa2d_best_prediction_filename(chain_id))
    write_model2d(best_model_path, mesh, best_resistivity; title = "Best $(chain_slug) model")
    write_data2d(best_data_path, best_data)

    iteration_records = MT2DIterationRecord[]
    T0 = config.temp_kappa
    ak = _resolve_ak(config.cool_ratio, config.max_iter)

    for iteration in 1:config.max_iter
        chi2_before_iter = current_fit.chi2
        rms_before_iter = current_fit.rms
        # one schedule drives both proposal width and acceptance
        T = _T_schedule(iteration; T0 = T0, ak = ak)

        #---------- generate and test the trials ----------
        # all trials branch from the iteration-start state, the lowest-rms trial
        # gets one deferred metropolis test (same scheme as the 3D driver)
        n_accepted_iter = 0
        last_accepted_trial = 0
        trial_cache = NamedTuple[]
        best_trial_rms = Inf
        best_trial_idx = 0
        best_trial_delta = similar(current_delta)
        best_trial_resistivity = current_resistivity
        best_trial_response = current_response
        best_trial_data = current_data
        best_trial_fit = current_fit
        for trial in 1:config.n_trials
            trial_delta = copy(current_delta)
            updated = _propose_controls!(
                trial_delta,
                T,
                config.log_bounds[1],
                config.log_bounds[2],
                base_control_values,
                n_update_controls,
                rng;
                step_scale = config.step_scale,
            )
            trial_model = _reconstruct_trial_model(
                start_log10,
                start_core_log10,
                delta_buffer,
                rbf_map,
                trial_delta,
                mesh,
                core_y,
                core_z,
                background_log10,
                config,
            )
            trial_response = run_mt2d_forward(mesh, trial_model.resistivity)
            trial_data = _predicted_data2d(trial_response, observed_data)
            trial_fit = chi2_rms2d(observed_data, trial_data)

            model_rel = ""
            if config.keep_models
                model_rel = joinpath(chain_slug, _vfsa2d_trial_model_filename(iteration, trial))
                write_model2d(joinpath(chain_dir, _vfsa2d_trial_model_filename(iteration, trial)), mesh, trial_model.resistivity; title = "VFSA $(chain_slug) trial model")
            end

            # relative reduced-χ² (rms²) energy against the iteration-start state;
            # a failed forward gives non-finite rms and can never win the iteration
            dE = (trial_fit.rms^2 - rms_before_iter^2) / max(rms_before_iter^2, eps())
            push!(
                trial_cache,
                (
                    trial = trial,
                    updated_controls = length(updated),
                    fit = trial_fit,
                    dE = dE,
                    p_acc = NaN,
                    u_acc = NaN,
                    accepted = false,
                    chi2_best = best_fit.chi2,
                    rms_best = best_fit.rms,
                    model_rel = model_rel,
                ),
            )
            if trial_fit.rms < best_trial_rms
                best_trial_rms = trial_fit.rms
                best_trial_idx = trial
                best_trial_delta .= trial_delta
                best_trial_resistivity = trial_model.resistivity
                best_trial_response = trial_response
                best_trial_data = trial_data
                best_trial_fit = trial_fit
            end
        end

        #---------- single metropolis test on the best trial ----------
        dE_best = (best_trial_rms^2 - rms_before_iter^2) / max(rms_before_iter^2, eps())
        u_acc_best = rand(rng)
        p_acc_best = dE_best <= 0 ? 1.0 : exp(-dE_best / max(T, 1e-12))
        accept_best = isfinite(best_trial_rms) && (u_acc_best < p_acc_best)
        if accept_best
            current_delta .= best_trial_delta
            current_resistivity = best_trial_resistivity
            current_response = best_trial_response
            current_data = best_trial_data
            current_fit = best_trial_fit
            n_accepted_iter = 1
            last_accepted_trial = best_trial_idx
            # record the global best on improvement, never feed it back
            if current_fit.chi2 < best_fit.chi2
                best_resistivity = copy(current_resistivity)
                best_response = current_response
                best_data = current_data
                best_fit = current_fit
                write_model2d(best_model_path, mesh, best_resistivity; title = "Best $(chain_slug) model")
                write_data2d(best_data_path, best_data)
            end
        end
        # backfill the winning trial's metropolis draw into its log row
        if best_trial_idx > 0
            tr = trial_cache[best_trial_idx]
            trial_cache[best_trial_idx] = merge(tr, (p_acc = p_acc_best, u_acc = u_acc_best, accepted = accept_best))
        end

        # Save best-so-far snapshot at configured intervals (for convergence GIF).
        if config.snapshot_interval > 0 && mod(iteration, config.snapshot_interval) == 0
            snap_file = joinpath(chain_dir, @sprintf("best_iter_%05d.rho", iteration))
            write_model2d(snap_file, mesh, best_resistivity; title = "$(chain_slug) best @ iter $(iteration)")
        end

        for trial in trial_cache
            _append_trial_row(
                trials_log;
                chain = chain_id,
                iteration = iteration,
                trial = trial.trial,
                T = T,
                chi2_before = chi2_before_iter,
                rms_before = rms_before_iter,
                chi2_prop = trial.fit.chi2,
                dchi2 = trial.fit.chi2 - chi2_before_iter,
                rms_prop = trial.fit.rms,
                drms = trial.fit.rms - rms_before_iter,
                dE = trial.dE,
                p_acc = trial.p_acc,
                u_acc = trial.u_acc,
                chi2_best = trial.chi2_best,
                rms_best = trial.rms_best,
                acc = trial.accepted ? 1 : 0,
                updated_controls = trial.updated_controls,
                model_rel = trial.model_rel,
                data_rel = "",
            )
        end

        # accepted_trial holds the winning trial index when its metropolis test
        # passed (0 otherwise)
        push!(
            iteration_records,
            MT2DIterationRecord(
                chain = chain_id,
                iteration = iteration,
                accepted_trial = last_accepted_trial,
                temperature = T,
                proposal_chi2 = current_fit.chi2,
                proposal_rms = current_fit.rms,
                current_chi2 = current_fit.chi2,
                current_rms = current_fit.rms,
                best_chi2 = best_fit.chi2,
                best_rms = best_fit.rms,
                delta_chi2 = current_fit.chi2 - chi2_before_iter,
                accepted = n_accepted_iter > 0,
            ),
        )

        _append_iteration_row(
            iteration_log;
            chain = chain_id,
            iteration = iteration,
            T = T,
            chi2_curr = current_fit.chi2,
            dchi2 = current_fit.chi2 - chi2_before_iter,
            rms_curr = current_fit.rms,
            current_rms = current_fit.rms,
            chi2_best = best_fit.chi2,
            rms_best = best_fit.rms,
            nacc = n_accepted_iter,
            model_rel = last_accepted_trial > 0 ?
                joinpath(chain_slug, _vfsa2d_trial_model_filename(iteration, last_accepted_trial)) : "",
        )
        _print_vfsa_progress_bar(chain_id, iteration, config.max_iter, current_fit.rms, best_fit.rms)

        # stop once best rms hits the target
        if best_fit.rms <= config.target_rms
            println()
            @info @sprintf("chain %d: target RMS %.3f reached at iter %d (best RMS %.5f); stopping.",
                           chain_id, config.target_rms, iteration, best_fit.rms)
            break
        end
    end

    MT2DChainResult(
        chain_id = chain_id,
        chain_dir = chain_dir,
        best_model_path = best_model_path,
        best_data_path = best_data_path,
        final_model_path = "",
        final_data_path = "",
        best_chi2 = best_fit.chi2,
        best_rms = best_fit.rms,
        best_resistivity = best_resistivity,
        best_response = best_response,
        best_data = best_data,
        final_chi2 = current_fit.chi2,
        final_rms = current_fit.rms,
        final_resistivity = current_resistivity,
        final_response = current_response,
        final_data = current_data,
        iterations = iteration_records,
        control_y_m = control_y_m,
        control_z_m = control_z_m,
    )
end

"""
    _write_iteration_snapshots(snapshot_dir, chains, mesh, config)

After all chains finish, read each chain's `best_iter_NNNNN.rho` snapshots,
compute the cross-chain mean resistivity at each snapshot iteration, and write
the averaged model to `snapshot_dir/avg_iter_NNNNN.rho`.
"""
function _write_iteration_snapshots(
    snapshot_dir::AbstractString,
    chains::Vector{MT2DChainResult},
    mesh::MT2DMesh,
    config::VFSA2DMTConfig,
)
    config.snapshot_interval <= 0 && return nothing
    mkpath(snapshot_dir)

    iterations = config.snapshot_interval:config.snapshot_interval:config.max_iter
    for iteration in iterations
        snap_filename = @sprintf("best_iter_%05d.rho", iteration)
        # Collect all available chain snapshots for this iteration.
        chain_models = Matrix{Float64}[]
        for chain in chains
            snap_path = joinpath(chain.chain_dir, snap_filename)
            if isfile(snap_path)
                model = load_model2d(snap_path)
                push!(chain_models, model.resistivity)
            end
        end
        isempty(chain_models) && continue
        # Average resistivity in log10 space (more meaningful for resistivity).
        avg_log10 = mean(log10.(m) for m in chain_models)
        avg_resistivity = 10.0 .^ avg_log10
        avg_path = joinpath(snapshot_dir, @sprintf("avg_iter_%05d.rho", iteration))
        write_model2d(avg_path, mesh, avg_resistivity; title = "Chain-mean model @ iter $(iteration)")
    end
    return nothing
end

function run_mt2d_vfsa(
    start_model_path::AbstractString,
    observed_data_path::AbstractString;
    run_dir::Union{Nothing, AbstractString} = nothing,
    true_model_path::Union{Nothing, AbstractString} = nothing,
    config::VFSA2DMTConfig = VFSA2DMTConfig(),
)
    run_info = run_dir === nothing ?
        _create_mt2d_run_dir(config.output_root) :
        (
            run_dir = String(run_dir),
            timestamp_slug = _vfsa2d_timestamp_slug(),
            timestamp_label = _vfsa2d_timestamp_label(),
            workflow = "VFSA2DMT",
        )
    mkpath(run_info.run_dir)
    output_dirs = _vfsa_output_dirs(run_info.run_dir)

    copied_observed_data_path = _ensure_file_in_dir(observed_data_path, output_dirs.data; target_name = _vfsa2d_observed_output_filename())
    _ensure_file_in_dir(start_model_path, output_dirs.models; target_name = "model.start")
    if true_model_path === nothing
        _maybe_copy_mt2d_true_model(observed_data_path, run_info.run_dir)
    elseif !isempty(String(true_model_path)) && isfile(String(true_model_path))
        _ensure_file_in_dir(String(true_model_path), run_info.run_dir; target_name = _vfsa2d_true_model_output_filename())
    end
    observed_data = load_data2d(copied_observed_data_path)
    start_model = load_model2d(start_model_path)
    mesh = build_mesh_from_model2d(
        start_model;
        frequencies = observed_data.frequencies,
        receiver_positions = observed_data.receivers,
    )
    start_resistivity = start_model.resistivity

    chains = MT2DChainResult[]
    for chain_id in 1:config.n_chains
        push!(chains, _run_mt2d_chain(chain_id, mesh, start_resistivity, observed_data, output_dirs, config))
    end
    best_chain = chains[argmin(getfield.(chains, :best_chi2))]
    ensemble = _summarize_best_chain_models(chains, mesh, observed_data)

    # Write per-iteration snapshot models (cross-chain mean of best-so-far).
    snapshot_dir = joinpath(run_info.run_dir, "snapshots")
    _write_iteration_snapshots(snapshot_dir, chains, mesh, config)

    (
        run_info = run_info,
        output_dirs = output_dirs,
        config = config,
        observed_data_path = copied_observed_data_path,
        snapshot_dir = snapshot_dir,
        chains = chains,
        best_chain = best_chain,
        ensemble = ensemble,
    )
end

# Writes the generic 2D VFSA summary for file-driven runs.
function _write_summary_2d_generic(
    path::AbstractString,
    config::VFSA2DMTConfig,
    result,
    output_paths,
)
    mkpath(dirname(path))
    lines = String[
        "# VFSA2DMT inversion summary",
        "",
        "- timestamp: $(result.run_info.timestamp_label)",
        "- n_chains: $(config.n_chains)",
        "- n_trials: $(config.n_trials)",
        "- max_iter: $(config.max_iter)",
        "- n_ctrl: $(config.n_ctrl)",
        "- log_bounds: $(config.log_bounds)",
        "- perturb_depth_m: $(config.perturb_depth_m)",
        "- best_chain: $(result.best_chain.chain_id)",
        "- best_chain_chi_square: $(round(result.best_chain.best_chi2, digits = 4))",
        "- best_chain_rms: $(round(result.best_chain.best_rms, digits = 4))",
        "- ensemble_mean_chi_square: $(round(result.ensemble.mean_fit.chi2, digits = 4))",
        "- ensemble_mean_rms: $(round(result.ensemble.mean_fit.rms, digits = 4))",
        "",
        "## Files",
        "- observed_data: $(basename(result.observed_data_path))",
        "- best_chain_model: $(basename(result.best_chain.best_model_path))",
        "- best_chain_data: $(basename(result.best_chain.best_data_path))",
        "- mean_model: $(basename(output_paths.mean_model_path))",
        "- median_model: $(basename(output_paths.median_model_path))",
        "- uncertainty_table: $(basename(output_paths.uncertainty_table_path))",
        "",
        "## Chains",
    ]

    for chain in result.chains
        acceptance_fraction = isempty(chain.iterations) ? 0.0 : mean(Float64.(getfield.(chain.iterations, :accepted)))
        push!(
            lines,
            "- chain $(chain.chain_id): best_chi_square=$(round(chain.best_chi2, digits = 4)), best_rms=$(round(chain.best_rms, digits = 4)), final_chi_square=$(round(chain.final_chi2, digits = 4)), final_rms=$(round(chain.final_rms, digits = 4)), acceptance_fraction=$(round(acceptance_fraction, digits = 4))",
        )
    end

    open(path, "w") do io
        println(io, join(lines, "\n"))
    end
    path
end

# Runs the file-based 2D VFSA inversion workflow from a starting model and observed data file.
function VFSA2DMT(
    start_model_path::AbstractString,
    observed_data_path::AbstractString;
    run_dir::Union{Nothing, AbstractString} = nothing,
    true_model_path::Union{Nothing, AbstractString} = nothing,
    config::VFSA2DMTConfig = VFSA2DMTConfig(),
)
    result = run_mt2d_vfsa(
        start_model_path,
        observed_data_path;
        run_dir = run_dir,
        true_model_path = true_model_path,
        config = config,
    )
    observed_data = load_data2d(result.observed_data_path)
    start_model = load_model2d(start_model_path)
    mesh = build_mesh_from_model2d(
        start_model;
        frequencies = observed_data.frequencies,
        receiver_positions = observed_data.receivers,
    )

    mean_model_path = write_model2d(
        joinpath(result.output_dirs.models, _vfsa2d_mean_model_filename()),
        mesh,
        result.ensemble.mean_resistivity;
        title = "Best-chain ensemble mean model",
    )
    median_model_path = write_model2d(
        joinpath(result.output_dirs.models, _vfsa2d_median_model_filename()),
        mesh,
        result.ensemble.median_resistivity;
        title = "Best-chain ensemble median model",
    )
    uncertainty_table_path = _write_mt2d_uncertainty_table(
        joinpath(result.output_dirs.models, _vfsa2d_std_model_filename()),
        mesh,
        result.ensemble,
    )
    plot_mt2d_model(mesh, result.ensemble.mean_resistivity; output_path = joinpath(result.output_dirs.plots, _vfsa2d_mean_model_plot_filename()))
    plot_mt2d_model(mesh, result.ensemble.median_resistivity; output_path = joinpath(result.output_dirs.plots, _vfsa2d_median_model_plot_filename()))
    plot_mt2d_data_maps(data_to_response2d(observed_data); output_path = joinpath(result.output_dirs.plots, _vfsa2d_observed_maps_plot_filename()))
    plot_mt2d_data_maps(result.best_chain.best_response; output_path = joinpath(result.output_dirs.plots, _vfsa2d_best_maps_plot_filename()))
    plot_mt2d_site_curves(
        data_to_response2d(observed_data);
        station_index = mt2d_center_station(mesh),
        output_path = joinpath(result.output_dirs.plots, _vfsa2d_observed_curves_plot_filename()),
    )
    plot_mt2d_site_curves(
        result.best_chain.best_response;
        station_index = mt2d_center_station(mesh),
        output_path = joinpath(result.output_dirs.plots, _vfsa2d_best_curves_plot_filename()),
    )
    convergence_plot_path = _plot_mt2d_vfsa_convergence(joinpath(result.output_dirs.plots, _vfsa2d_convergence_plot_filename()), result.chains)
    summary_path = _write_summary_2d_generic(
        joinpath(result.output_dirs.logs, "Summary.md"),
        result.config,
        result,
        (
            mean_model_path = mean_model_path,
            median_model_path = median_model_path,
            uncertainty_table_path = uncertainty_table_path,
        ),
    )

    merge(
        result,
        (
            mean_model_path = mean_model_path,
            median_model_path = median_model_path,
            mean_data_path = "",
            uncertainty_table_path = uncertainty_table_path,
            convergence_plot_path = convergence_plot_path,
            summary_path = summary_path,
        ),
    )
end

function _parse_cli_args(args::AbstractVector{<:AbstractString})
    values = Dict{String, String}()
    index = 1
    while index <= length(args)
        arg = args[index]
        if startswith(arg, "--") && index < length(args)
            values[arg] = args[index + 1]
            index += 2
        else
            index += 1
        end
    end
    values
end

function _parse_cli_int(cli::Dict{String, String}, key::AbstractString, default::Int)
    haskey(cli, key) || return default
    parsed = tryparse(Int, cli[key])
    parsed === nothing && error("expected integer value for $key, got '$(cli[key])'")
    parsed
end

function _parse_cli_float(cli::Dict{String, String}, key::AbstractString, default::Float64)
    haskey(cli, key) || return default
    parsed = tryparse(Float64, cli[key])
    parsed === nothing && error("expected numeric value for $key, got '$(cli[key])'")
    parsed
end

function VFSA2DMT(params::VFSA2DMTParams)
    script_dir = dirname(abspath(params.script_path))
    run_info = _create_mt2d_run_dir(script_dir; run_name = params.run_name)
    VFSA2DMT(
        abspath(params.start_model_path),
        abspath(params.data_path);
        run_dir = run_info.run_dir,
        true_model_path = isempty(params.model_path) ? nothing : abspath(params.model_path),
        config = params.config,
    )
end

# Prints the required and optional inputs for the 2D VFSA workflow.
function _print_usage_2d(io::IO = stdout)
    d = VFSA2DMTConfig()
    println(io, "VFSA2DMT required inputs:")
    println(io, "  1. start_model_path")
    println(io, "  2. observed_data_path")
    println(io, "VFSA2DMT optional inputs with defaults:")
    println(io, "  --n-chains $(d.n_chains)")
    println(io, "  --n-trials $(d.n_trials)")
    println(io, "  --max-iter $(d.max_iter)")
    println(io, "  --n-ctrl $(d.n_ctrl)")
    println(io, "  --log-bounds $(d.log_bounds[1]),$(d.log_bounds[2])")
    println(io, "  --perturb-depth $(d.perturb_depth_m)")
    println(io, "  --output-root $(d.output_root)")
    println(io, "  --seed $(d.seed)")
end

# Parses the log10 resistivity bounds from the CLI.
function _parse_log_bounds_2d(cli::Dict{String, String}, default_value::Tuple{Float64, Float64})
    haskey(cli, "--log-bounds") || return default_value
    parts = split(cli["--log-bounds"], ",")
    length(parts) == 2 || error("expected --log-bounds lo,hi")
    (parse(Float64, strip(parts[1])), parse(Float64, strip(parts[2])))
end

# Runs the CLI entry point for the 2D VFSA workflow.
function main_vfsa2dmt(args::AbstractVector{<:AbstractString} = ARGS)
    length(args) < 2 && return _print_usage_2d()

    d = VFSA2DMTConfig()
    start_model_path = args[1]
    observed_data_path = args[2]
    cli = _parse_cli_args(args[3:end])
    config = VFSA2DMTConfig(
        n_chains = _parse_cli_int(cli, "--n-chains", d.n_chains),
        n_trials = _parse_cli_int(cli, "--n-trials", d.n_trials),
        max_iter = _parse_cli_int(cli, "--max-iter", d.max_iter),
        n_ctrl = _parse_cli_int(cli, "--n-ctrl", d.n_ctrl),
        log_bounds = _parse_log_bounds_2d(cli, d.log_bounds),
        perturb_depth_m = _parse_cli_float(cli, "--perturb-depth", d.perturb_depth_m),
        seed = _parse_cli_int(cli, "--seed", d.seed),
        output_root = get(cli, "--output-root", d.output_root),
    )
    result = VFSA2DMT(start_model_path, observed_data_path; config = config)

    println("VFSA2DMT outputs:")
    println("  run_dir = ", result.run_info.run_dir)
    println("  best_chi_square = ", round(result.best_chain.best_chi2, digits = 4))
    println("  ensemble_rms = ", round(result.ensemble.mean_fit.rms, digits = 4))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_vfsa2dmt()
end

#= 
This is an experimental MCMC code for solving 1D MT inverse problem. Don't read too
    much into this one at the momement. I will add proper documentation later
    when this works in a way I want it to. 
    -- @pankajkmishra 
=# 


using Distributed
addprocs(11)

@everywhere begin
    cd(@__DIR__)
    using Random
    using LinearAlgebra
    using Distributions
    using Turing
    using Statistics
    using LinearSolve
    using MCMCChains
    using AbstractMCMC
    include("MT1D_response_FD.jl")
end

function read_mesh(filename::String)
    lines = readlines(filename)
    [parse(Float64, line) for line in lines[2:end]]
end

function read_model(filename::String)
    lines = readlines(filename)
    [parse(Float64, line) for line in lines[2:end]]
end

@everywhere function project_main_model_to_mesh(main_boundaries, main_resistivities, mesh, main_depth)
    n_cells = length(mesh) - 1
    model_out = similar(main_resistivities, n_cells)
    for i in 1:n_cells
        z_mid = 0.5 * (mesh[i] + mesh[i+1])
        if z_mid <= main_depth
            idx = searchsortedfirst(main_boundaries, z_mid)
            idx = clamp(idx, 1, length(main_resistivities))
            model_out[i] = main_resistivities[idx]
        else
            model_out[i] = main_resistivities[end]
        end
    end
    return model_out
end

@everywhere function generate_synthetic_data(freqs, mesh, true_model, noise_rhoa, noise_phi, main_depth, n_main_layers)
    boundary_log = range(log10(1.0), log10(main_depth+1), length=n_main_layers+1)
    main_boundaries = 10 .^ boundary_log .- 1
    main_boundaries[1] = 0.0
    main_boundaries[end] = main_depth
    fine_model = project_main_model_to_mesh(main_boundaries, true_model, mesh, main_depth)
    pred_rhoa, pred_phase = MT1D_response_FD(freqs, mesh, fine_model)
    obs_rhoa = pred_rhoa .+ noise_rhoa .* randn(length(freqs))
    obs_phase = pred_phase .+ noise_phi .* randn(length(freqs))
    return obs_rhoa, obs_phase
end

@everywhere @model function MT_Nlayer_model(freqs, obs_rhoa, obs_phase, noise_rhoa, noise_phi, mesh, main_depth, n_main_layers, true_model, prior_std)
    res = Vector{Float64}(undef, n_main_layers)
    for j in 1:n_main_layers
        log_true = log(true_model[j])
        res[j] ~ Normal(log_true, prior_std)
    end
    res = exp.(res)
    boundary_log = range(log10(1.0), log10(main_depth+1), length=n_main_layers+1)
    main_boundaries = 10 .^ boundary_log .- 1
    main_boundaries[1] = 0.0
    main_boundaries[end] = main_depth
    fine_model = project_main_model_to_mesh(main_boundaries, res, mesh, main_depth)
    pred_rhoa, pred_phase = MT1D_response_FD(freqs, mesh, fine_model)
    for i in 1:length(freqs)
        obs_rhoa[i] ~ Normal(pred_rhoa[i], noise_rhoa)
        obs_phase[i] ~ Normal(pred_phase[i], noise_phi)
    end
end

function compute_posterior_statistics(chain)
    total_iter, num_params, num_chains = size(chain.value)
    flattened = reshape(chain.value, total_iter*num_chains, num_params)
    phys_samples = exp.(flattened)
    mean_model = mean(phys_samples, dims=1)[:]
    std_model = std(phys_samples, dims=1)[:]
    return mean_model, std_model
end

function save_model_to_ubc(filename::String, model::Vector{Float64})
    if isempty(model)
        error("Model is empty. Cannot write to file.")
    end
    open(filename, "w") do io
        println(io, length(model))
        for res in model
            println(io, res)
        end
    end
end

function run_inversion_Nlayer(n_main_layers; prior_std=0.3)
    mesh = read_mesh("model.mesh")
    true_model = read_model("model.rho")
    if length(true_model) < length(mesh)-1
        main_depth = mesh[length(true_model)+1]
    else
        main_depth = maximum(mesh)
    end
    freqs = 10 .^ range(-3, 3, length=50)
    noise_rhoa = 5.0
    noise_phi = 2.5
    obs_rhoa, obs_phase = generate_synthetic_data(freqs, mesh, true_model, noise_rhoa, noise_phi, main_depth, n_main_layers)
    mymodel = MT_Nlayer_model(freqs, obs_rhoa, obs_phase, noise_rhoa, noise_phi, mesh, main_depth, n_main_layers, true_model, prior_std)
    nchains = 10
    nsamples = 6000
    println("Sampling with $nchains chains, $nsamples samples each (including burn-in).")
    spinner = ["|", "/", "-", "\\"]
    stop_flag = false
    function spinner_thread()
        i = 1
        while !stop_flag
            print("\rSampling... ", spinner[i])
            flush(stdout)
            i = (i % length(spinner)) + 1
            sleep(0.1)
        end
        print("\rSampling... Done!   \n")
    end
    t = Threads.@spawn spinner_thread()
    chain = sample(mymodel, MH(), MCMCDistributed(), nsamples, nchains)
    stop_flag = true
    wait(t)
    println("MCMC diagnostics:")
    display(describe(chain))
    mean_model, std_model = compute_posterior_statistics(chain)
    boundary_log = range(log10(1.0), log10(main_depth+1), length=n_main_layers+1)
    main_boundaries = 10 .^ boundary_log .- 1
    main_boundaries[1] = 0.0
    main_boundaries[end] = main_depth
    mean_model_interp = project_main_model_to_mesh(main_boundaries, mean_model, mesh, main_depth)
    std_model_interp = project_main_model_to_mesh(main_boundaries, std_model, mesh, main_depth)
    save_model_to_ubc("model.mean", mean_model_interp)
    save_model_to_ubc("model.std", std_model_interp)
    println("Interpolated models saved in UBC format: 'model.mean' and 'model.std'.")
    return mean_model_interp, std_model_interp
end

mean_model_result, std_model_result = run_inversion_Nlayer(10)
println("Done!")
rmprocs(workers())
GC.gc()

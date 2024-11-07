# Import necessary packages
using Flux
using Flux: DataLoader, Optimise
using Statistics
using Random
using Plots
using LinearAlgebra
using BSON: @save, @load  # For saving and loading models and data
using GR  # Use the GR backend for headless plotting

GR.inline("png")  # Set GR to save plots directly as images

# 1. Define the MT1D function with complex handling
function MT1D(T, rho, h)
    mu0 = 4 * π * 1e-7
    m = length(rho)
    
    # Initialize arrays as Complex{Float64} explicitly
    k = Complex{Float64}[0.0 + 0.0im for _ in 1:m, _ in 1:length(T)]
    Z = Complex{Float64}[0.0 + 0.0im for _ in 1:length(T)]
    
    for n in 1:m
        k[n, :] = sqrt.(-im * 2 * π * mu0 ./ (T .* rho[n]))
    end
    
    Z .= -(im * mu0 * 2 * π) ./ (T .* k[m, :])  # Ensure Z remains complex
    
    for n in (m-1):-1:1
        A = -(im * mu0 * 2 * π) ./ (T .* k[n, :])
        B = exp.(-2 .* k[n, :] .* h[n])
        Z .= A .* (A .* (1 .- B) .+ Z .* (1 .+ B)) ./ (A .* (1 .+ B) .+ Z .* (1 .- B))
    end
    
    # Apparent resistivity and phase should be extracted as real values
    rho_a = real.((T ./ (mu0 * 2 * π)) .* abs.(Z).^2)  # Extract real part element-wise
    phase = real.(-atan.(imag.(Z) ./ real.(Z)) .* 180 / π)  # Extract real part element-wise
    
    return rho_a, phase
end

# 2. Generate Synthetic MT Data Using MT1D Model
# Function to generate data using MT1D model response
function generate_MT1D_data(num_samples::Int)
    Random.seed!(1234)  # For reproducibility
    T = 10 .^ range(-3, 3, length=100)  # Period range (0.001 to 1000 s)
    
    # Generate random resistivity and thickness values for each sample
    data_X = []
    data_y = []
    
    for _ in 1:num_samples
        rho = 10 .^ (rand(4) .* 3 .+ 1)  # Four layers with resistivity (10-1000 ohm-m)
        h = rand(3) .* 1000  # Three layers with thickness (0-1000 m)
        
        rho_a, phase = MT1D(T, rho, h)
        
        # Store (resistivity, thickness) as input and apparent resistivity and phase as output
        push!(data_X, vcat(rho, h))
        push!(data_y, vcat(rho_a, phase))
    end
    
    return hcat(data_X...), hcat(data_y...)
end

# Generate 1000 synthetic data samples
X, y = generate_MT1D_data(1000)

# Split dataset into training and validation sets
train_size = Int(0.8 * size(X, 2))
train_X, train_y = X[:, 1:train_size], y[:, 1:train_size]
val_X, val_y = X[:, train_size+1:end], y[:, train_size+1:end]

# 3. Define the Neural Network Model for MT Inversion
# A basic feedforward neural network
model = Chain(
    Dense(7, 32, relu),   # Input layer with 7 features (4 resistivities + 3 thicknesses)
    Dense(32, 64, relu),  # Hidden layer
    Dense(64, 64, relu),  # Hidden layer
    Dense(64, 200)        # Output layer with 200 values (100 apparent resistivity + 100 phase)
)

# 4. Train the Neural Network
# Define loss function and optimizer
loss(x, y) = Flux.Losses.mse(model(x), y)
opt = Flux.Adam(0.001)  # Use Flux.Adam if Optimisers.Adam is not available

# Training loop
train_data = DataLoader((train_X, train_y), batchsize=32, shuffle=true)

println("Training started...")
for epoch in 1:10000
    # Train over batches
    Flux.train!(loss, Flux.params(model), train_data, opt)
    
    # Validation loss
    val_loss = loss(val_X, val_y)
    println("Epoch: $epoch, Validation Loss: $val_loss")
end
println("Training completed.")

# Save trained model and training data
@save "mt_inversion_model.bson" model
@save "training_data.bson" train_X train_y val_X val_y

# 5. Evaluate the Model on Test Data
# Generate new synthetic test data for evaluation
test_X, test_y = generate_MT1D_data(200)

# Predict with the trained model
predictions = model(test_X)

# Extract predicted apparent resistivity and phase
predicted_rho_a = predictions[1:100, :]
predicted_phase = predictions[101:200, :]

# Extract true apparent resistivity and phase from test data
test_rho_a = test_y[1:100, :]
test_phase = test_y[101:200, :]

# Plot true vs predicted values for the test set (for one sample)
T = 10 .^ range(-3, 3, length=100)
Plots.plot(T, test_rho_a[:, 1], xscale=:log10, yscale=:log10, label="True Apparent Resistivity", linewidth=2)
Plots.plot!(T, predicted_rho_a[:, 1], xscale=:log10, yscale=:log10, label="Predicted Apparent Resistivity", linewidth=2, linestyle=:dash)
Plots.xlabel!("Period (s)")
Plots.ylabel!("Apparent Resistivity (Ohm-m)")
Plots.title!("MT Inversion: True vs Predicted Apparent Resistivity (Test Set)")
Plots.savefig("MT_Inversion_Test_Set_rhoa.png")

# Plot true vs predicted phase for the test set (for one sample)
Plots.plot(T, test_phase[:, 1], xscale=:log10, label="True Phase", linewidth=2)
Plots.plot!(T, predicted_phase[:, 1], xscale=:log10, label="Predicted Phase", linewidth=2, linestyle=:dash)
Plots.xlabel!("Period (s)")
Plots.ylabel!("Phase (degrees)")
Plots.title!("MT Inversion: True vs Predicted Phase (Test Set)")
Plots.savefig("MT_Inversion_Test_Set_phase.png")

# Plot the MT model used for prediction for the test set (for one sample)
mt_rho = test_X[1:4, 1]  # Extract resistivities from test data
mt_h = test_X[5:7, 1]    # Extract thicknesses from test data
mt_depth = [0; cumsum(mt_h); 2000.0]  # Add an arbitrary depth for the last layer
Plots.plot(mt_depth, [mt_rho; mt_rho[end]], xlabel="Depth (m)", ylabel="Resistivity (Ohm-m)", xflip=true, label="MT Model", linewidth=2)
Plots.title!("MT Model for the Test Sample")
Plots.savefig("MT_Model_Test_Sample.png")

# Save test data and predictions
@save "test_data_and_predictions.bson" test_X test_y predictions

# Tutorials

This section provides practical examples and tutorials for using MTGeophysics.jl.

## Basic Usage

### MT 1D Forward Modeling

The most common use case is computing magnetotelluric responses for 1D layered earth models.

#### Analytical Solution

```julia
using MTGeophysics

# Define a simple 3-layer earth model
frequencies = 10.0 .^ range(-2, 2, length=50)  # 0.01 to 100 Hz
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]  # depths in meters
resistivities = [100.0, 10.0, 1000.0]  # resistivity in Ω⋅m

# Compute MT response using analytical solution
rho_a, phi = MT1D_response(frequencies, layer_boundaries, resistivities)

println("Computed $(length(frequencies)) frequency points")
println("Frequency range: $(minimum(frequencies)) - $(maximum(frequencies)) Hz")
println("Apparent resistivity range: $(round(minimum(rho_a), digits=2)) - $(round(maximum(rho_a), digits=2)) Ω⋅m")
println("Phase range: $(round(minimum(phi), digits=2))° - $(round(maximum(phi), digits=2))°")
```

#### Finite-Difference Solution

For more complex models or when numerical stability is required, use the finite-difference method:

```julia
using MTGeophysics

frequencies = 10.0 .^ range(-2, 2, length=50)
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]
resistivities = [100.0, 10.0, 1000.0]

# Compute MT response using finite-difference diagonal method
rho_a, phi = MT1D_response_FDD(frequencies, layer_boundaries, resistivities)
```

If LinearSolve is available, you can also use:

```julia
# Requires LinearSolve package
rho_a, phi = MT1D_response_FD(frequencies, layer_boundaries, resistivities)
```

### UBC Format File Creation

Create mesh and model files in UBC format for use with other geophysical software:

```julia
using MTGeophysics

# Define layer boundaries and resistivities
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]
resistivities = [100.0, 10.0, 1000.0]

# Generate UBC mesh and model files
# d0: initial cell size, r: growth factor, pad: number of padding cells
UBC_1D(layer_boundaries, resistivities; d0=10.0, r=1.1, pad=15)

# This creates two files: model.mesh and model.rho
if isfile("model.mesh") && isfile("model.rho")
    println("✓ Created UBC format files: model.mesh and model.rho")
end
```

### Working with ModEM Data

#### Creating Empty Data Structures

```julia
using MTGeophysics

# Create empty ModEM data structure
data = make_nan_data()
println("Created ModEMData structure")
println("Initial sites: $(data.ns)")
println("Initial frequencies: $(data.nf)")
```

#### Loading ModEM Data Files

```julia
using MTGeophysics

# Load observed data
data_obs = load_data_modem("path/to/observed_data.dat")
println("Loaded $(data_obs.ns) sites with $(data_obs.nf) frequencies")

# Load predicted data
data_pred = load_data_modem("path/to/predicted_data.dat")
```

#### Calculating Apparent Resistivity and Phase

ModEM data files store impedance tensor elements. You can calculate apparent resistivity and phase:

```julia
# The data structure automatically computes rho and pha fields
# Access them directly:
println("Apparent resistivity: ", data_obs.rho)
println("Phase: ", data_obs.pha)
```

### Inversion Quality Metrics

Calculate chi-squared and RMS values to assess inversion quality:

```julia
using MTGeophysics

# Calculate chi-squared and RMS between observed and predicted data
result = chi2_and_rms("observed_data.dat", "predicted_data.dat";
                      use_impedance=true, 
                      use_tipper=false, 
                      components=["ZXY", "ZYX"])

println("Chi-squared: $(result.chi2)")
println("RMS: $(result.rms)")
```

You can control which data components are included:

```julia
# Include all impedance components and tipper
result = chi2_and_rms("observed_data.dat", "predicted_data.dat";
                      use_impedance=true,
                      use_tipper=true)

# Only specific impedance components
result = chi2_and_rms("observed_data.dat", "predicted_data.dat";
                      use_impedance=true,
                      components=["ZXX", "ZXY", "ZYX", "ZYY"])
```

### 3D Model Visualization

!!! note
    Model visualization requires the GLMakie package to be installed.

#### Loading and Viewing Models

```julia
using MTGeophysics
using GLMakie

# Load a ModEM model
model = load_model_modem("path/to/model.rho")

# Create interactive 3D viewer
fig, parts = gl_modem_viewer(model)

# Display the figure
display(fig)
```

#### Customizing Visualization

```julia
# Use logarithmic scale (default)
fig, parts = gl_modem_viewer(model; 
                              log10scale=true,
                              resistivity_range=(1.0, 4.0))

# Use linear scale
fig, parts = gl_modem_viewer(model;
                              log10scale=false,
                              resistivity_range=(10, 1e4))

# Remove padding cells and limit depth
fig, parts = gl_modem_viewer(model;
                              withPadding=false,
                              max_depth=50000.0,
                              resistivity_range=[1.5, 3.5])
```

The viewer provides:
- Three orthogonal slice planes (YZ, XZ, XY)
- Interactive sliders to move through the model
- Toggles to show/hide each slice
- Color bar showing resistivity values
- Support for both log10 and linear scales

## Advanced Examples

### Comparing Different Forward Modeling Methods

```julia
using MTGeophysics

frequencies = 10.0 .^ range(-2, 2, length=50)
layer_boundaries = [0.0, 500.0, 2000.0, 10000.0]
resistivities = [100.0, 10.0, 1000.0]

# Analytical solution
rho_analytical, phi_analytical = MT1D_response(frequencies, layer_boundaries, resistivities)

# Finite-difference diagonal solution
rho_fdd, phi_fdd = MT1D_response_FDD(frequencies, layer_boundaries, resistivities)

# Compare results
max_rho_diff = maximum(abs.(rho_analytical .- rho_fdd))
max_phi_diff = maximum(abs.(phi_analytical .- phi_fdd))

println("Maximum apparent resistivity difference: $(max_rho_diff) Ω⋅m")
println("Maximum phase difference: $(max_phi_diff)°")
```

### Batch Processing Multiple Models

```julia
using MTGeophysics

frequencies = 10.0 .^ range(-2, 2, length=50)

# Define multiple models
models = [
    ([0.0, 500.0, 2000.0], [100.0, 10.0]),
    ([0.0, 1000.0, 3000.0], [50.0, 500.0]),
    ([0.0, 300.0, 1000.0, 5000.0], [200.0, 20.0, 2000.0])
]

# Process each model
results = []
for (i, (boundaries, resistivities)) in enumerate(models)
    rho_a, phi = MT1D_response(frequencies, boundaries, resistivities)
    push!(results, (model=i, rho_a=rho_a, phi=phi))
    println("Processed model $i")
end
```

## Tips and Best Practices

### Frequency Range Selection

- For shallow structures: Use higher frequencies (> 1 Hz)
- For deep structures: Use lower frequencies (< 0.1 Hz)
- Typical MT surveys: 0.001 Hz to 1000 Hz

### Layer Boundaries

- First boundary should always be 0.0 (surface)
- Last boundary should be deep enough to simulate a half-space
- Use logarithmic spacing for better resolution at depth

### Resistivity Values

- Typical sediments: 1-100 Ω⋅m
- Typical crystalline rocks: 100-10,000 Ω⋅m
- Conductive layers: < 10 Ω⋅m
- Resistive layers: > 1,000 Ω⋅m

### Visualization Performance

- For large 3D models, consider using `withPadding=false` to focus on the core region
- Use `max_depth` parameter to limit visualization to depths of interest
- Adjust `resistivity_range` to enhance contrast in your area of interest

## Next Steps

- Explore the [API Reference](@ref) for detailed function documentation
- Check the [GitHub repository](https://github.com/JuliaGeophysics/MTGeophysics.jl) for the latest updates
- Report issues or request features on [GitHub Issues](https://github.com/JuliaGeophysics/MTGeophysics.jl/issues)

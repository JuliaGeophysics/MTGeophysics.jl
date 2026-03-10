# 1D Forward Modelling

Compute the 1D magnetotelluric forward response for a layered earth model using analytical and finite-difference solvers.

## Running the example

```bash
julia --project=. Examples/response_1d.jl
```

This reads the default benchmark model (`Examples/0Layered1D/Layered1D.true`) and data specification (`Examples/0Layered1D/Layered1D.ref`), computes the predicted response, and writes a comparison plot.

Custom input files:

```bash
julia --project=. Examples/response_1d.jl path/to/model path/to/dataspec
```

## From Julia

```julia
using MTGeophysics

# File-based: solve and write predicted data
pred_path = ForwardSolve1D("Examples/0Layered1D/Layered1D.true",
                           "Examples/0Layered1D/Layered1D.ref")

# Plot observed vs predicted
PlotData1D("Examples/0Layered1D/Layered1D.ref", pred_path;
           output_path = "DataPlot1D.png")

# Plot the model structure
PlotModel1D("Examples/0Layered1D/Layered1D.true";
            output_path = "ModelPlot1D.png")
```

## Building models programmatically

```julia
using MTGeophysics

# Define layers
thicknesses   = [120.0, 280.0, 650.0, 1400.0]     # metres
resistivities = [100.0, 20.0, 350.0, 40.0, 800.0]  # Ω·m (last = basement)

# Build mesh
mesh = BuildMesh1D(thicknesses, resistivities)

# Compute responses
frequencies = 10 .^ range(-3, 3, length=50)
r_analytical = solve_mt1d_analytical(frequencies, resistivities, thicknesses)
r_fd         = solve_mt1d_fd(frequencies, mesh)

# Plot
plot_mt1d_data(Dict(:Analytical => r_analytical, :FD => r_fd);
               output_path = "comparison.png")
plot_mt1d_model(mesh; output_path = "model.png")
```

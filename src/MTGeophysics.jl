"""
# MTGeophysics.jl

Julia package for magnetotelluric (MT) geophysics research and applications.

## Features
- MT 1D forward modeling with analytical and finite-difference methods
- ModEM data loading and processing
- Model visualization and plotting
- Chi-squared and RMS calculations for inversions
- UBC format mesh and model handling

## Main exports
- `load_data_modem`: Load ModEM data files
- `MT1D_response`: Analytical MT 1D forward modeling  
- `MT1D_response_FDD`: Finite-difference MT 1D forward modeling
- `UBC_1D`: Create UBC format mesh and model files
- `chi2_and_rms`: Calculate chi-squared and RMS for inversion analysis
- `read_mackie3d_model`, `load_model_modem`: Load 3D models
- `gl_modem_viewer`: 3D model visualization
"""
module MTGeophysics

using LinearAlgebra
using Statistics
using DelimitedFiles

# Include submodules
include("data/ModEMData.jl")
include("model/ModEMModel.jl") 
include("forward/MT1D_response.jl")
include("forward/MT1D_response_FDD.jl")
include("forward/UBC_1D.jl")
include("inversion/Chi2RMS.jl")

# Conditionally include LinearSolve-dependent functionality
has_linearsolve = false
try
    using LinearSolve
    include("forward/MT1D_response_FD.jl")
    global has_linearsolve = true
catch LoadError
    @warn "LinearSolve not available, some finite-difference functionality disabled"
end

# Conditionally include visualization if GLMakie is available
has_visualization = false
try
    using GLMakie
    include("visualisation/PlotModel.jl")
    global has_visualization = true
    # Export visualization function
    export gl_modem_viewer
catch LoadError
    @warn "GLMakie not available, visualization functionality disabled"
end

# Export main data structures
export ModEMData, ModEMModel

# Export data loading functions  
export load_data_modem, make_nan_data, calc_rho_pha

# Export model loading functions
export read_mackie3d_model, load_model_modem

# Export forward modeling functions
export MT1D_response, MT1D_response_FDD

# Conditionally export LinearSolve-dependent functions
if has_linearsolve
    export MT1D_response_FD
end

# Export mesh/model creation functions  
export UBC_1D

# Export inversion analysis functions
export chi2_and_rms

end # module
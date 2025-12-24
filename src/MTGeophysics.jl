"""
# MTGeophysics.jl

Julia package for magnetotelluric (MT) geophysics research and applications.

## Features
- MT 1D forward modeling with analytical and finite-difference methods
- ModEM data loading and processing
- Model visualization and plotting
- Chi-squared and RMS calculations for inversions
"""
module MTGeophysics

using LinearAlgebra
using Statistics

# Include submodules
include("io/ModEMData.jl")
include("io/ModEMModel.jl")
include("fwd/MT1D.jl")
include("inv/Chi2RMS.jl")

# Conditionally include visualization if GLMakie is available
has_visualization = false
try
    using GLMakie
    include("viz/PlotModel.jl")
    global has_visualization = true
    export gl_modem_viewer
catch LoadError
    @warn "GLMakie not available, interactive visualization functionality disabled"
end

# Export main data structures
export ModEMData, ModEMModel

# Export data loading functions  
export load_data_modem, make_nan_data, calc_rho_pha

# Export model loading functions
export read_mackie3d_model, load_model_modem

# Export forward modeling functions
export MT1D, MT1DMethod, Analytical, FiniteDifference

# Export inversion analysis functions
export chi2_and_rms

end # module
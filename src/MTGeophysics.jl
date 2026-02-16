# MTGeophysics package module entry.
# Author: @pankajkmishra
# This file includes all core source scripts and defines the exported public API.
# It also conditionally enables visualization helpers when GLMakie is available.

module MTGeophysics

using LinearAlgebra
using Statistics

include("Data.jl")
include("Model.jl")
include("MT1D.jl")
include("Chi2RMS.jl")

has_visualization = false
try
    using GLMakie
    include("PlotModel.jl")
    global has_visualization = true
    export edges_from_centers, core_indices, z_indices_for_max_depth
    export compute_colorrange, prepare_model_arrays
catch LoadError
    @warn "GLMakie not available, interactive visualization functionality disabled"
end

export Data, Model, ModEMData, ModEMModel

export load_data_modem, make_nan_data, calc_rho_pha

export read_mackie3d_model, load_model_modem, write_model_modem

export MT1D, MT1DMethod, Analytical, FiniteDifference

export chi2_and_rms

end

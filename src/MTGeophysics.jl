# MTGeophysics package module entry.
# Author: @pankajkmishra
# This file includes all core source scripts and defines the exported public API.
# It also conditionally enables visualization helpers when GLMakie is available.

module MTGeophysics

using LinearAlgebra
using Statistics

#----- 3-D ModEM I/O and misfit (existing code) ---------------------------#

include("Data.jl")
include("Model.jl")
include("Chi2RMS.jl")

#----- Headless core/padding utilities (always available) ------------------#

include("CoreUtils3D.jl")

#----- WS3D format model I/O (log10 internal) -----------------------------#

include("WS3DModel.jl")

#----- Visualization (optional, requires GLMakie) -------------------------#

has_visualization = false
try
    using GLMakie
    include("PlotModel.jl")
    global has_visualization = true
    export compute_colorrange, prepare_model_arrays
catch LoadError
    @warn "GLMakie not available, interactive visualization functionality disabled"
end

#----- 1-D / 2-D forward, I/O, plotting, and inversion (patch) -----------#

include("MTGeophysics1D.jl")
include("MTGeophysics2D.jl")
include("VFSA2DMT.jl")

#----- 3-D VFSA inversion and ensemble statistics -------------------------#

include("VFSA3DMT.jl")

#----- Exports: 3-D ModEM ------------------------------------------------#

export Data, Model, ModEMData, ModEMModel
export load_data_modem, make_nan_data, calc_rho_pha
export read_mackie3d_model, load_model_modem, write_model_modem
export chi2_and_rms

#----- Exports: Core utilities (always available) -------------------------#

export edges_from_centers, core_indices, z_indices_for_max_depth
export lateral_core_ranges, core_view

#----- Exports: WS3D model I/O -------------------------------------------#

export WS3DModel
export load_ws3d_model, read_ws3d_model, write_ws3d_model

#----- Exports: 3-D VFSA inversion ---------------------------------------#

export VFSA3DMTConfig
export VFSA3DMT
export AnalyseEnsemble3D
export core_statistics
export RBFMap, build_rbf_map, apply_rbf_map!

#----- Exports: 1-D / 2-D ------------------------------------------------#

export MT1DMesh
export MT2DMesh
export MT1DDataSpec
export MT1DResponse
export MT2DResponse
export ModelFile2D
export DataFile2D
export FitSummary2D
export VFSA2DMTConfig
export VFSA2DMTParams

export BuildMesh1D
export BuildMesh2D
export MakeMesh1D
export MakeMesh2D

export Forward1D
export Forward2D
export ForwardSolve1D
export ForwardSolve2D

export load_mt1d_model
export write_mt1d_model
export load_mt1d_data_spec
export write_mt1d_data_template
export load_mt1d_observed_data
export write_mt1d_observed_data
export mt1d_layered_model
export solve_mt1d_analytical
export solve_mt1d_fd

export load_model2d
export write_model2d
export build_mesh_from_model2d
export build_default_mt2d_mesh
export build_mt2d_halfspace_model
export build_mt2d_data_template
export write_mt2d_data_template
export load_data2d
export write_data2d
export data_from_response2d
export data_to_response2d
export chi2_rms2d
export run_mt2d_forward

export plot_mt1d_data
export plot_mt1d_model
export plot_mt2d_data_maps
export plot_mt2d_site_curves
export plot_mt2d_model
export PlotData1D
export PlotModel1D
export PlotData2D
export PlotModel2D

export AnalyseEnsemble2D
export VFSA2DMT

end

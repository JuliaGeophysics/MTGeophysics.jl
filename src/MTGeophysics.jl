# MTGeophysics package module entry.
# Author: @pankajkmishra
# This file includes all core source scripts and defines the exported public API.
# It also conditionally enables visualization helpers when GLMakie is available.

module MTGeophysics

using LinearAlgebra
using Statistics
using Dates
using Printf
using Shapefile
using GeoInterface
using Proj

#----- 3-D ModEM I/O and misfit (existing code) ---------------------------#

include("Data.jl")
include("Rotate.jl")
include("Model.jl")
include("Chi2RMS.jl")
include("Distortion.jl")

#----- Headless core/padding utilities (always available) ------------------#

include("CoreUtils3D.jl")

#----- Headless shapefile overlay utilities --------------------------------#

include("ShapefileOverlay.jl")
include("GeoRef3D.jl")

#----- Phase tensors, induction vectors and their GIS export ---------------#

include("PhaseTensor.jl")

#----- WS3D format model I/O (log10 internal) -----------------------------#

include("WS3DModel.jl")

#----- Topography / bathymetry extraction and air / water masks -----------#

include("Mask3D.jl")
include("MeshToMesh.jl")
include("MakeMesh3D.jl")

#----- Visualization (optional, requires GLMakie) -------------------------#

try
    using GLMakie
    include("PlotModel.jl")
    include("PlotModel3D.jl")
    include("EditModel3D.jl")
    include("MakeMesh3DGUI.jl")
    include("PlotData3D.jl")
    export compute_colorrange, prepare_model_arrays
    export PlotModelXY, PlotModelXZ, PlotModelYZ, PlotModelXYZ
    export EditModelByLayers, EditModelByDrawing
    export PlotPTIVMap
catch LoadError
    @warn "GLMakie not available, interactive visualization functionality disabled"
end

#----- 1-D / 2-D forward, I/O, plotting, and inversion ---------------------#

include("Control2D.jl")
include("Mesh2D.jl")
include("Fwd2D.jl")
include("Strike2D.jl")
include("Topo2D.jl")
include("Fwd1D.jl")
include("PlotModel2D.jl")
include("PlotData2D.jl")
include("Inv2D.jl")
include("Inv2D_GN.jl")
include("Inv2D_NLCG.jl")
include("Inv2D_VFSA.jl")
include("Inv1D.jl")
include("Mask2D.jl")
include("MakeMesh2D.jl")
include("MakeMesh2DGUI.jl")

#----- 3-D VFSA inversion and ensemble statistics -------------------------#

include("VFSA3DMT.jl")

#----- Exports: 3-D ModEM ------------------------------------------------#

export Data, Model, ModEMData, ModEMModel
export load_data_modem, write_data_modem, make_nan_data, calc_rho_pha
export rotate_data, RotationStep
export read_mackie3d_model, load_model_modem, write_model_modem
export chi2_and_rms
export chi2_and_rms_distorted, DistortionFit, write_distortion_file

#----- Exports: Core utilities (always available) -------------------------#

export edges_from_centers, core_indices, z_indices_for_max_depth
export lateral_core_ranges, core_view

#----- Exports: Shapefile overlay utilities --------------------------------#

export detect_shapefile_crs, shapefile_coord_transform
export load_shapefile_geometries, prepare_shapefiles
export pick_shapefile

#----- Exports: Phase tensors and induction vectors ------------------------#

export phase_tensor, induction_vector, has_tipper_data
export phase_tensors_from_data, induction_vectors_from_data
export write_ptiv_gis

#----- Exports: mesh design and mesh-to-mesh projection -------------------#

export MakeMesh3D
export MeshToMesh

#----- Exports: WS3D model I/O -------------------------------------------#

export WS3DModel
export load_ws3d_model, read_ws3d_model, write_ws3d_model

#----- Exports: Bathymetry / water masks ----------------------------------#

export extract_bathymetry, write_bathymetry, read_bathymetry
export water_mask_from_bathymetry, water_mask_from_model
export extract_topography, write_topography, read_topography
export air_mask_from_topography, air_mask_from_model

#----- Exports: 3-D VFSA inversion ---------------------------------------#

export VFSA3DMTConfig
export VFSA3DMT
export AnalyseEnsemble3D
export core_statistics
export RBFMap, build_rbf_map, apply_rbf_map!

#----- Exports: 1-D / 2-D ------------------------------------------------#

export MT2DMesh
export MT2DResponse
export ModelFile2D
export DataFile2D
export FitSummary2D
export VFSA2DConfig
export VFSA2D, mt2d_ensemble

export BuildMesh2D

export Forward2D
export ForwardSolve1D, MakeMesh1D, mt1d_layers, mt1d_skin_depth, mt1d_impedance, mt1d_frechet, mt1d_site_data, WriteFrechet1D
export Invert1D, InvCtrl1D, ReadInvCtrl1D, WriteInvCtrl1D, PlotInversion1D
export Mask2D, mt2d_ground, MakeMesh2D
export ForwardSolve2D
export EstimateStrike2D, RotateData2D, StrikeData2D, RotateToStrike2D


export build_default_mt2d_mesh
export mt2d_skin_depth, mt2d_skin_depth_layers
export build_mt2d_halfspace_model
export build_mt2d_data_template
export write_mt2d_data_template
export load_data2d
export write_data2d
export data_from_response2d
export chi2_rms2d
export FrechetDerivative2D, ApplyFrechet2D, ApplyFrechetTranspose2D
export Invert2D, Inv2DOptions, Inv2DResult, AbstractInversion2D
export inv2d_frechet, inv2d_gradient
export inv2d_tag, inv2d_init, inv2d_prepare!, inv2d_direction, inv2d_reject!, inv2d_accept!, inv2d_info, inv2d_validate
export GaussNewton2D, GaussNewton2DConfig, GaussNewton2DResult
export NLCG2D, NLCG2DConfig, NLCG2DResult
export run_mt2d_forward
export FwdCtrl2D, InvCtrl2D, VFSACtrl2D, Cov2D
export ReadFwdCtrl2D, WriteFwdCtrl2D, ReadInvCtrl2D, WriteInvCtrl2D, ReadVFSACtrl2D, WriteVFSACtrl2D
export ReadCov2D, WriteCov2D, ReadMask2D, WriteMask2D
export ReadModel2D, WriteModel2D, Mesh2DFromInputs, mt2d_air_layers, mt2d_geometric_layers
export WriteFrechet2D
export Topo2D, ReadTopo2D, WriteTopo2D, Topography2D, mt2d_profile_topography
export mt2d_air_mask, mt2d_station_offsets, mt2d_receiver_depths, mt2d_receiver_columns, mt2d_topo_air

export plot_mt1d_model, plot_mt1d_convergence
export plot_mt2d_site_curves
export plot_mt2d_model
export plot_mt2d_mesh
export plot_mt2d_data_fit
export plot_inv2d_convergence, plot_vfsa2d_convergence
export PlotInversion2D
export PlotModel1D
export PlotData2D
export PlotModel2D

export AnalyseEnsemble2D

end

# API Reference

This page documents all public functions and types in MTGeophysics.jl.

!!! note
    Full API documentation with docstrings will be available once the package dependencies are installed.

## Data Structures

- `ModEMData` - Data structure for ModEM format data
- `ModEMModel` - Data structure for ModEM format models

## Data Loading and Creation

- `load_data_modem(path)` - Load ModEM data files
- `make_nan_data()` - Create empty ModEM data structure
- `calc_rho_pha` - Calculate apparent resistivity and phase

## Model Loading

- `read_mackie3d_model(path)` - Read Mackie 3D model files
- `load_model_modem(path)` - Load ModEM model files

## Forward Modeling

### 1D MT Response

- `MT1D_response(frequencies, boundaries, resistivities)` - Analytical 1D MT forward modeling
- `MT1D_response_FDD(frequencies, boundaries, resistivities)` - Finite-difference diagonal 1D MT forward modeling  
- `MT1D_response_FD(frequencies, boundaries, resistivities)` - Finite-difference 1D MT forward modeling (requires LinearSolve)

### File Generation

- `UBC_1D(boundaries, resistivities; kwargs...)` - Create UBC format mesh and model files

## Inversion Analysis

- `chi2_and_rms(obs_path, pred_path; kwargs...)` - Calculate chi-squared and RMS metrics

## Visualization

- `gl_modem_viewer(model; kwargs...)` - Interactive 3D model visualization (requires GLMakie)

```@raw html
---
layout: home

hero:
  name: "MTGeophysics.jl"
  tagline: A software repository for magnetotelluric geophysics research and applications.
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: Forward modelling
      link: /forward/2d
    - theme: alt
      text: View on GitHub
      link: https://github.com/JuliaGeophysics/MTGeophysics.jl

features:
  - title: Forward modelling
    details: Exact 1D layered-earth responses and a 2D TE/TM finite-difference solver with topography and lakes. 3D runs through ModEM.
    link: /forward/2d
  - title: Deterministic inversion
    details: Gauss–Newton and NLCG in 1D and 2D, set up from a handful of small text files.
    link: /inversion/2d_deterministic
  - title: VFSA inversion
    details: Very fast simulated annealing in 1D, 2D and 3D. Several chains give an ensemble of models and an estimate of uncertainty.
    link: /inversion/2d_vfsa
  - title: Data and meshes
    details: Read and rotate ModEM data files, and build 2D and 3D meshes with topography and water from a data file.
    link: /data/files
  - title: Interactive 3D viewers
    details: Depth slices, cross-sections, phase tensor maps and model editors, with shapefile overlays and GIS export.
    link: /visualisation/3d_models
---
```

## What is MTGeophysics.jl?

MTGeophysics.jl is part of the [JuliaGeophysics ecosystem](https://github.com/JuliaGeophysics) and is intended for both research and real-world applications. It provides reusable forward-modelling, inversion, and visualization components so you can prototype new machine-learning methods, inversion strategies, and data-analysis ideas quickly without rebuilding core MT tooling from scratch. More broadly, JuliaGeophysics aims to build a tightly integrated yet modular ecosystem for multiphysics workflows, multisource data integration, and uncertainty quantification.

## Install

```julia
pkg> add MTGeophysics
```

Then see [Getting Started](getting_started.md) to run your first forward model and inversion.

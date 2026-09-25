```@raw html
---
layout: home

hero:
  name: "MTGeophysics.jl"
  text: "Magnetotellurics in Julia"
  tagline: Forward modelling, inversion and visualisation of MT data in 1D, 2D and 3D, using the same files as ModEM.
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

MTGeophysics.jl is a Julia package for magnetotelluric (MT) research and applications. It gives you
the building blocks of an MT workflow (reading data, building meshes, forward modelling, inversion
and plotting) so that you can test new ideas without rebuilding the basic tools first.

The package reads and writes the same data and model files as
[ModEM](https://github.com/magnetotellurics/ModEM), so you can move models and data back and forth
between the two. It is part of the [JuliaGeophysics](https://github.com/JuliaGeophysics) ecosystem.

## Install

```julia
pkg> add MTGeophysics
```

Then see [Getting Started](getting_started.md) to run your first forward model and inversion.

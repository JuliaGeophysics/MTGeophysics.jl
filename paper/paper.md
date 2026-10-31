---
title: 'MTGeophysics.jl: A software repository for magnetotelluric research and applications'
tags:
  - Julia
  - magnetotellurics
  - geophysics
  - forward modelling
  - inversion
  - uncertainty quantification
authors:
  - name: Pankaj K Mishra
    orcid: 0000-0003-4907-4724
    affiliation: 1
    corresponding: true
affiliations:
  - name: Geological Survey of Finland (GTK), Espoo, Finland
    index: 1
date: 10 March 2026
bibliography: paper.bib
---

# Summary

The magnetotelluric (MT) method uses natural electromagnetic field variations
to image subsurface electrical resistivity at depths ranging from tens of
metres to hundreds of kilometres. It is widely used in geophysical
exploration for natural resources, including minerals, geothermal energy, and
groundwater, as well as in the characterisation of sites for CO$_2$ storage and
in studies of crustal and lithospheric structure. The MT inverse problem is,
however, inherently non-unique: many distinct resistivity models reproduce the
observations equally well. Conventional
practice nonetheless reports a single regularised model, with no quantitative
measure of its uncertainty.

MTGeophysics.jl is a Julia package developed to address this limitation. Its
principal aim is the quantification of uncertainty in MT model building,
advancing interpretation from a single deterministic model, through an
ensemble of plausible models, towards the full posterior distribution
(\autoref{fig:vision}). The package also provides forward modelling, inversion
and visualisation tools for MT research and applications.

![The vision behind MTGeophysics.jl. Deterministic inversion yields a single
regularised model $\hat{\mathbf{m}}$; an ensemble of plausible models
$\{\mathbf{m}_i\}$ and, ultimately, the full posterior distribution
$p(\mathbf{m}\mid\mathbf{d}^{\mathrm{obs}})$ characterise model uncertainty
with increasing rigour and computational cost. The package aims to make
ensemble-based uncertainty quantification routine and to provide a practical
path towards the full posterior.\label{fig:vision}](vision.png){ width=98% }

# Statement of need

Sampling many models consistent with the data, rather than a single
regularised model, is computationally demanding in 3D, because a regional
mesh carries millions of free parameters and every candidate model requires a
full electromagnetic forward solve.

MTGeophysics.jl addresses this gap with a stochastic
inversion workflow based on Very Fast Simulated Annealing (VFSA)
[@SenStoffa2013] that produces an *ensemble* of plausible 3D models,
summarised by its mean, median, and standard deviation. The workflow reuses the community-standard [ModEM](https://github.com/magnetotellurics/ModEM) forward
solver [@Egbert2012; @Kelbert2014] by reading and writing its native file
formats, so it slots into established MT practice without asking users to
change their models, data, or solver.

The package targets MT researchers and students who want an open-source
toolkit that builds on Julia's [@Bezanson2017] strengths in numerical
computing, composability, and interactive graphics. It is designed as a
research repository and as a core component of a broader JuliaGeophysics
ecosystem: forward solvers, data structures, and inversion routines are meant
to be reused and recombined, so that new ideas can be prototyped without
rebuilding core MT tooling.

# State of the field

Several open-source codes support MT research. ModEM
[@Egbert2012; @Kelbert2014] is widely used for deterministic 3D MT inversion
on structured grids. [FEMTIC](https://github.com/yoshiya-usui/femtic) [@Usui2015; @Usui2017] performs 3D MT inversion on
unstructured meshes, incorporating topography, bathymetry and galvanic
distortion, and [FEMTIC-DABIC](https://github.com/hsong-28/FEMTIC-DABIC) [@Song2025] extends it with statistical selection
of the regularisation parameter. [jif3D](https://sourceforge.net/projects/jif3d/) [@Moorkamp2011] is a multi-physics
joint inversion framework that combines 3D MT with gravity, magnetic and
seismic data. [MARE2DEM](https://mare2dem.bitbucket.io/master/) [@Key2016]
provides adaptive finite-element 2D modelling and inversion of MT and
controlled-source electromagnetic data, and [custEM](https://gitlab.com/Rochlitz.R/custEM) [@Rochlitz2019] offers 3D
finite-element electromagnetic forward modelling. [MTpy](https://github.com/MTgeophysics/mtpy-v2) [@Kirkby2019] provides
Python utilities for MT data handling and visualisation, while [pyGIMLi](https://www.pygimli.org)
[@Ruecker2017] and [SimPEG](https://simpeg.xyz) [@Cockett2015] are general geophysical inversion
frameworks that include MT modules. [HiQGA.jl](https://github.com/GeoscienceAustralia/HiQGA.jl) [@Ray2021] provides deterministic
and probabilistic (trans-dimensional Bayesian) inversion for several
geophysical methods.

MTGeophysics.jl complements these tools, with a focus on stochastic inversion
and ensemble-based uncertainty quantification that is practical in 3D. Because
many 3D MT models and datasets are held in ModEM formats, the package reads and
writes those files directly and uses ModEM as its 3D forward solver, so it fits
into existing workflows without converting meshes or data. Building instead on
a larger framework would have required that conversion, which motivated a
standalone package. The package also includes its own pure-Julia 1D and 2D
forward solvers and inversion routines. It is intended to work alongside the existing
open-source ecosystem in magnetotelluric research.

# Software design

The central design choice in MTGeophysics.jl is reduced model
parameterisation, the strategy that makes stochastic 3D inversion tractable.
Stochastic inversion becomes intractable if every mesh cell is a free
parameter: the search space is too large and each proposal requires an
expensive forward solve. Rather than perturbing the full grid, the VFSA
workflow perturbs a small set of radial-basis-function (RBF) control points and
maps them back to the mesh with a compactly supported Gaussian RBF (truncated
at $3\sigma$), as shown for one iteration in \autoref{fig:loop}.

![One VFSA iteration as a sparse-model update loop. The padded modelling mesh
(1) is reduced to its core (2), sampled at $M$ random control points, the
*sparse model* (3), perturbed by the VFSA rule (4), then mapped back to the
full padded mesh by compactly supported Gaussian-RBF interpolation (5) for the
ModEM forward solve, advancing $\mathbf{m}_k \to
\mathbf{m}_{k+1}$.\label{fig:loop}](vfsa_loop.png){ width=90% }

The trade-off is resolution: structure finer than the control-point spacing
and RBF width cannot be represented, so both are user-set parameters. The
resulting ensemble consists of independently optimised models rather than
samples from a posterior distribution; its spread indicates where the data
weakly constrain the model rather than providing calibrated uncertainty. A
low-dimensional representation is also what would make fully Bayesian sampling
affordable in 3D.

Around this core, the package is organised as a single Julia module with
separate layers: ModEM 3D data and model I/O; a 1D module that solves the
layered-earth impedance recursion exactly, with Fréchet derivatives by
automatic differentiation; a 2D module that assembles sparse finite-difference
operators for the TE and TM modes; and inversion modules for 1D, 2D and 3D:
deterministic inversion, which in 3D is left to ModEM, and stochastic inversion
by VFSA, with adaptive simulated annealing (ASA) as a further 3D option. Below
3D the package uses its own forward solvers and runs chains in parallel
threads; in 3D it
calls ModEM as a subprocess, runs each chain as an independent job and computes
ensemble statistics across the completed runs. Because ModEM is reached only
through its files, the forward solver can be swapped. Visualisation uses
[GLMakie](https://github.com/MakieOrg/Makie.jl) [@Danisch2021] for interactive 3D slice viewers, with coordinate reprojection via
[Proj.jl](https://github.com/JuliaGeo/Proj.jl) and shapefile overlays, and [CairoMakie](https://github.com/MakieOrg/Makie.jl) for static figures. All model,
data and control files are plain text, and argument-driven example scripts make
every workflow reproducible from the command line. Detailed usage is described
in the [package documentation](https://juliageophysics.github.io/MTGeophysics.jl/dev/).

Development is continuing in two directions: portable, GPU-accelerated forward
solvers written in pure Julia, which would remove the dependence on an
external solver, and neural surrogates of the forward problem, which together
aim to make full Bayesian inversion computationally practical. These
components are experimental and not part of the current release. Implicit
neural representations, which have been applied to 3D gravity inversion
[@Mishra2026INR], offer a further route to a compact, continuous model
parameterisation.

MTGeophysics.jl is intended to support the latest Julia release and the one
before it, currently Julia 1.13 and 1.12. Continuous integration runs the test
suite on both versions under a virtual display, so the visualisation layer is
also exercised. The tests cover the 3D data structures and ModEM I/O, analytic
checks of the 1D solver, physical-consistency and finite-difference Fréchet
checks of the 2D solver, the COMMEMI benchmark generators, and end-to-end 1D
and 2D inversions.

# Benchmarks and validation

The package is validated at two levels. The COMMEMI 2D benchmark models
[@Zhdanov1997] are generated natively by `helpers/benchmarks_2D.jl`, so the 2D
forward and inversion workflows can be exercised from the repository without
external downloads. At regional scale, the 3D VFSA workflow has been applied
to USArray MT data from the Cascadia subduction zone [@Patro2008] and compared
with a deterministic ModEM inversion [@Mishra2026]: the ensemble mean recovers
the major conductive structures, while the ensemble standard deviation
indicates where the data poorly constrain the model. The comparison is
reproduced in the [package documentation](https://juliageophysics.github.io/MTGeophysics.jl/dev/).

# Research impact statement

MTGeophysics.jl supports ongoing magnetotelluric research at the Geological
Survey of Finland (GTK) and opens uncertainty-aware MT interpretation to the
wider community. The package has been in development on GitHub since August
2024 and has attracted external users, who have reported issues, requested
features and contributed fixes.

The package and the research built on it have been presented to the
geophysics and scientific-computing communities by the author and
collaborators, including at the EGU General Assembly 2026 [@Mishra2026], the 27th
Electromagnetic Induction Workshop (EMIW 2026), and JuliaCon 2026.
MTGeophysics.jl was also used to participate in the
[MT3DINV-4 workshop](https://mt3dinv4.mtnet.info/Real_data.html)
(Memorial University of Newfoundland, 2026), a community 3D MT inversion
benchmarking exercise in which many groups inverted a common real field
dataset acquired over the Raglan mining district in northern Quebec, Canada.

# AI usage disclosure

Different versions of Claude (Anthropic) and Codex (OpenAI) were used over the
duration of development of this package to implement ideas, reorganise the
repository, write documentation and tests, and edit this paper. The GitHub
Copilot coding agent contributed a small number of commits, and is also used
to check for continuous integration (CI) errors on GitHub; Dependabot is used
to keep dependency versions current. All AI-generated code and text were
reviewed, tested, and verified by the author for correctness.

# Acknowledgements

This work was supported by the Research Council of Finland (project 359261).
The author wishes to acknowledge CSC – IT Center for Science, Finland, for
computational resources. The author also thanks the electromagnetic geophysics team at the Geological Survey of Finland for discussions on magnetotelluric geophysics.

# References

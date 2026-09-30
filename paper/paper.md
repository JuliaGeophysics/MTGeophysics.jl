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
(\autoref{fig:vision}). In addition, the package provides modern capabilities
for forward modelling, inversion, and visualization in support of
state-of-the-art MT research and applications.

![The vision behind MTGeophysics.jl. Deterministic inversion yields a single
regularised model $\hat{\mathbf{m}}$; an ensemble of plausible models
$\{\mathbf{m}_i\}$ and, ultimately, the full posterior distribution
$p(\mathbf{m}\mid\mathbf{d}^{\mathrm{obs}})$ characterise model uncertainty
with increasing rigour and computational cost. The package aims to make
ensemble-based uncertainty quantification routine and to provide a practical
path towards the full posterior.\label{fig:vision}](vision.png){ width=98% }

# Statement of need

Sampling many models consistent with the data, rather than a single
regularised model, has been considered impractical in 3D, because a regional
mesh carries millions of free parameters and every candidate model requires a
full electromagnetic forward solve.

MTGeophysics.jl closes this gap in practice. It provides a stochastic
inversion workflow based on Very Fast Simulated Annealing (VFSA)
[@SenStoffa2013] that produces an *ensemble* of plausible 3D models,
summarised by its mean, median, and standard deviation. The workflow reuses the community-standard ModEM forward
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

Several open-source tools address parts of the MT workflow. ModEM
[@Egbert2012; @Kelbert2014] is the community standard for deterministic 3D MT
inversion but provides no built-in visualization or uncertainty
quantification. MARE2DEM [@Key2016] provides adaptive finite-element 2D
modelling and inversion of MT and controlled-source electromagnetic data, but
does not address 3D problems. MTpy [@Kirkby2019]
provides comprehensive Python utilities for MT data handling and
visualization but no forward solvers or stochastic inversion. pyGIMLi
[@Ruecker2017] and SimPEG [@Cockett2015] are general Python inversion
frameworks with MT modules, but their MT functionality is embedded in much
larger codebases. We are not aware of one that integrates stochastic inversion
with ensemble uncertainty quantification in both 2D and 3D. Implementing such a
workflow within those frameworks would mean adopting their own meshes and
forward solvers, whereas 3D MT models and data in practice are largely held in
ModEM formats and inverted with ModEM. A standalone package operating directly
on those files, with its own lightweight 1D and 2D solvers for rapid
experimentation, was therefore the more practical route.

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

The same principle extends beyond VFSA: searching a low-dimensional
representation is a prerequisite for practical probabilistic inversion. Fully
Bayesian
approaches such as Markov chain Monte Carlo mirror the multi-chain structure
already used here, but become affordable in 3D only under a similarly compact
model. The roadmap follows the same logic: differentiable open-source forward
solvers, neural surrogates for the expensive 3D solve, and implicit neural
representations that recast the reduced parameterisation as a learned
continuous field. That combination has already
been demonstrated for 3D gravity inversion [@Mishra2026INR], and Julia's
scientific-machine-learning ecosystem makes MT the natural next target.

Around this core, the package is organised as a single Julia module with
clearly separated layers: ModEM 3D data/model I/O; a 1D module that solves the
layered-earth impedance recursion exactly, with Fréchet derivatives by
automatic differentiation; a 2D module that assembles sparse
finite-difference operators for the TE/TM Maxwell equations and solves them
with LinearSolve.jl; and the inversion modules, with Gauss–Newton in 1D and
2D, NLCG in 2D, and VFSA in 1D, 2D and 3D. In 1D and 2D the package uses its
own forward engines and runs VFSA chains in parallel threads; in 3D it wraps
the external ModEM solver, each chain runs as an independent job, and ensemble
statistics are computed across the completed runs. Because the driver reaches
ModEM only through its files and a subprocess call, the forward engine is
swappable. Extending this interface to open-source 3D MT solvers is planned,
which would remove the one component users must currently obtain separately.
Interactive visualization is handled through
GLMakie (GPU-accelerated 3D slice viewers with coordinate reprojection via
Proj.jl and optional shapefile overlays), while CairoMakie produces
publication-quality static plots. Both are regular dependencies; the
interactive viewers need OpenGL and are disabled with a warning if GLMakie
cannot be loaded, while the static plots need no display. All model, data and
control files are plain text for reproducibility and version control, and
argument-driven example scripts make every workflow reproducible from the
command line. In future, this scriptable, file-based interface is intended to
support the development of agentic AI systems, in which AI agents could set
up, run, and interpret MT modelling and inversion workflows. Detailed usage (the
interactive slice viewers and the polygon-based model editor) is documented in
the package documentation [@MTGeophysicsDocs].

Correctness is guarded by continuous integration on GitHub Actions. Every push
and pull request runs the full test suite on Linux against the two most recent
Julia releases (at the time of writing, 1.13 and 1.12), so the package tracks
the latest toolchain while the previous release remains supported.
Because GLMakie cannot precompile without a display, the workflow provisions a
virtual framebuffer (Xvfb) so the visualization layer is exercised headlessly
in a clean environment. The suite combines checks of the core 3D data
structures and ModEM I/O with analytic checks of the 1D solver, physical
consistency and finite-difference Fréchet checks of the 2D solver, the COMMEMI
benchmark generators, and end-to-end 1D and 2D inversions (Gauss–Newton, NLCG
and VFSA), so solver and inversion behaviour is re-verified on every change. Companion workflows build and deploy the documentation and keep
dependency bounds current.

# Benchmarks and validation

The package is validated at two levels. The COMMEMI 2D benchmark models
[@Zhdanov1997] are generated natively by `helpers/benchmarks_2D.jl`, so the 2D
forward and inversion workflows can be exercised from the repository without
external downloads. At regional scale, the 3D VFSA workflow has been applied
to USArray MT data from the Cascadia subduction zone [@Patro2008] and compared
with a deterministic ModEM inversion [@Mishra2026]: the ensemble mean recovers
the major conductive structures, while the ensemble standard deviation
indicates where the data poorly constrain the model. The comparison is
reproduced in the package documentation [@MTGeophysicsDocs].

# Research impact statement

MTGeophysics.jl supports ongoing magnetotelluric research at the Geological
Survey of Finland (GTK) and opens uncertainty-aware MT interpretation to the
wider community.

The package and the research built on it have been presented to the
geophysics and scientific-computing communities by Pankaj K Mishra and
Cedric Patzer, including at the EGU General Assembly 2026, the 27th
Electromagnetic Induction Workshop (EMIW 2026), and JuliaCon 2026.
MTGeophysics.jl was also used to participate in the
[MT3DINV-4 workshop](https://mt3dinv4.mtnet.info/Real_data.html)
(Memorial University of Newfoundland, 2025), a community 3D MT inversion
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
computational resources. Author wishes to thank Jochen Kamm, Cedric Patzer, Uula Autio for discussion on various aspects of magnetotelluric geophysics.

# References

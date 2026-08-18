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
to image electrical resistivity tens to hundreds of kilometres into the Earth.
Resistivity is diagnostic of the constituents that matter most for resources
and hazards: saline fluids, partial melt, graphite, and metallic sulphides. MT
therefore underpins exploration for critical minerals and geothermal energy,
characterisation of aquifers, CO$_2$ storage sites and repository host rocks,
and studies of volcanic and tectonic systems. That sensitivity, however, is
ambiguous: a conductor at depth may be brine, melt, or ore. Conductance also
trades off against depth, so many different subsurface configurations explain
the same observations equally well. Standard practice nonetheless reports a
single "best" model, with no measure of which imaged features the data require
and which merely survive the modeller's choice of regularisation.

MTGeophysics.jl is a Julia package for MT modelling, inversion, and
interactive visualization built around this problem. It provides the usual
building blocks (1D and 2D forward solvers, ModEM 3D I/O [@Kelbert2014], and
interactive 3D viewers), but is organised around a single guiding idea: making
stochastic inversion and uncertainty quantification practical for MT. The
package is designed to move MT interpretation along the spectrum from one
deterministic model, to an ensemble of plausible models, toward the full
posterior distribution (\autoref{fig:vision}).

![The vision behind MTGeophysics.jl. Deterministic inversion returns a single
model $\hat{\mathbf{m}}$ chosen by a regularisation weight $\lambda$ that the
user picks; the regulariser expresses a *preference*, not information, so the
"best" model shifts as $\lambda$ changes. Greater rigour in uncertainty
quantification (a VFSA ensemble $\{\mathbf{m}_i\}$, ultimately the full
posterior $p(\mathbf{m}\mid\mathbf{d}^{\mathrm{obs}})$) costs more
computation. The package is built to make the ensemble regime routine at
regional scale and to open a practical path toward the full
posterior.\label{fig:vision}](vision.png){ width=98% }

# Statement of need

Standard 3D MT inversion answers non-uniqueness by regularisation: it returns
the single model minimising data misfit plus a smoothness penalty whose weight
the user picks. The alternative is to sample many models consistent with the
data. In 3D this has been considered impractical, because a regional mesh
carries millions of free parameters and every candidate model requires a full
electromagnetic forward solve.

MTGeophysics.jl closes this gap in practice. It provides a stochastic
inversion workflow based on Very Fast
Simulated Annealing (VFSA) [@SenStoffa2013] that produces an *ensemble* of
plausible 3D models, from which ensemble statistics (mean, median, standard
deviation) give a per-voxel picture of what the data do and do not constrain.
The workflow reuses the community-standard ModEM forward
solver [@Egbert2012; @Kelbert2014] by reading and writing its native file
formats, so it slots into established MT practice without asking users to
change their models, data, or solver.

The package targets MT researchers and students who want an open-source
toolkit that builds on Julia's [@Bezanson2017] strengths in numerical
computing, composability, and interactive graphics. It is designed as a
research repository and as a core component of a broader JuliaGeophysics
ecosystem: forward solvers, data structures, and inversion routines are meant
to be reused and recombined, so that new ideas can be prototyped without
rebuilding core MT tooling. That composability also positions the package for
the scientific-machine-learning direction set out below, and its
plain-text formats and argument-driven example scripts keep every workflow
reproducible from the command line and drivable by external tooling.

# State of the field

Several open-source tools address parts of the MT workflow. ModEM
[@Egbert2012; @Kelbert2014] is the community standard for 3D MT inversion but
provides no built-in visualization or 1D/2D forward capabilities. MARE2DEM
[@Key2016] focuses on 2D and 3D marine electromagnetic modelling but is
Fortran-based and limited to controlled-source methods. MTpy [@Kirkby2019]
provides comprehensive Python utilities for MT data handling and
visualization but no forward solvers or stochastic inversion. pyGIMLi
[@Ruecker2017] and SimPEG [@Cockett2015] are general Python inversion
frameworks with MT modules, but their MT functionality is embedded in much
larger codebases. We are not aware of one that integrates stochastic inversion
with ensemble uncertainty quantification in both 2D and 3D. MTGeophysics.jl is
distinguished by the abstraction its components are organised around. A
reduced model parameterisation makes ensemble inversion affordable at regional
scale, and the package reuses the forward solver the MT community already
trusts.

# Software design

The central design choice in MTGeophysics.jl is reduced model
parameterisation, the strategy that makes stochastic 3D inversion tractable.
Stochastic inversion becomes intractable if every mesh cell is a free
parameter: the search space is too large and each proposal requires an
expensive forward solve. Rather than perturbing the full grid, the VFSA
workflow perturbs a small set of radial-basis-function (RBF) control points and
maps them back to the mesh with a compactly supported Gaussian RBF (truncated
at $3\sigma$). Each iteration reduces the padded modelling mesh to its core,
samples $M$ random control points, perturbs them by the VFSA rule, and
interpolates back to the full mesh for the ModEM forward solve
(\autoref{fig:loop}). Inverting for a handful of coefficients is what makes
ensemble exploration feasible at regional scale.

![One VFSA iteration as a sparse-model update loop. The padded modelling mesh
(1) is reduced to its core (2), sampled at $M$ random control points, the
*sparse model* (3), perturbed by the VFSA rule (4), then mapped back to the
full padded mesh by compactly supported Gaussian-RBF interpolation (5) for the
ModEM forward solve, advancing $\mathbf{m}_k \to \mathbf{m}_{k+1}$. This
reduced parameterisation is the strategy that makes stochastic 3D inversion
tractable.\label{fig:loop}](vfsa_loop.png){ width=90% }

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
clearly separated layers: ModEM 3D data/model I/O; a 1D module with analytical
and finite-difference solvers; a 2D module that assembles sparse
finite-difference operators for the TE/TM Maxwell equations and solves them
with LinearSolve.jl; and the VFSA inversion module. In 2D the package includes
its own finite-difference forward engine; in 3D it wraps the external ModEM
solver, running multiple independent Markov chains in parallel and computing
ensemble statistics across them. The inversion driver reaches that solver only
through plain-text files and a subprocess call, so the forward engine is
swappable. Extending this interface to open-source 3D MT solvers is planned,
which would remove the one component users must currently obtain separately.
Interactive visualization is handled through
GLMakie (GPU-accelerated 3D slice viewers with coordinate reprojection via
Proj.jl and optional shapefile overlays), while CairoMakie produces
publication-quality static plots. The visualization layer is conditionally
loaded, so the core package runs without OpenGL dependencies. All file formats
are plain text for reproducibility and version control. Detailed usage (the
interactive slice viewers and the polygon-based model editor) is documented in
the package documentation [@MTGeophysicsDocs].

Correctness is guarded by continuous integration on GitHub Actions. Every push
and pull request runs the full test suite on Linux across two Julia versions:
the minimum release declared in `Project.toml` (1.10) and a pinned current
release, so both the supported floor and the latest toolchain are exercised.
Because GLMakie cannot precompile without a display, the workflow provisions a
virtual framebuffer (Xvfb) so the visualization layer is exercised headlessly
in a clean environment. The suite combines smoke tests of the core 3D data structures and
ModEM I/O with regression checks of the 1D and 2D forward solvers, including
the COMMEMI benchmark generators, so forward-solver results are re-verified on
every change. Companion workflows build and deploy the documentation and keep
dependency bounds current.

# Benchmarks and validation

The package is validated at two levels. The 2D forward solver and VFSA
inversion reproduce the community-standard COMMEMI 2D benchmark suite
[@Zhdanov1997], generated natively by `helpers/benchmarks_2D.jl` and
reproducible from the repository without external downloads. At regional
scale, the 3D VFSA inversion is benchmarked against USArray MT data from the
Cascadia subduction zone [@Patro2008], a dataset extensively studied with
the deterministic ModEM NLCG inversion [@Mishra2026]. The VFSA ensemble mean
recovers the major conductive structures found by deterministic inversion,
while the ensemble standard deviation provides a per-voxel uncertainty
estimate that no single deterministic inversion can offer; high-variance
regions flag where the data poorly constrain the model. The benchmark is a
substantial computation: each non-greedy chain performs roughly 12,000 ModEM
forward solves at about 13.5 s each, close to two days per chain, with five
chains run in parallel. The slice-by-slice
comparison of deterministic model, ensemble mean, and ensemble standard
deviation is reproduced in the package documentation [@MTGeophysicsDocs].

# Research impact statement

MTGeophysics.jl supports ongoing magnetotelluric research at the Geological
Survey of Finland (GTK) and opens uncertainty-aware MT interpretation to the
wider community. Its components are reusable building blocks for the broader
JuliaGeophysics ecosystem.

The package and the research built on it have been presented to the
geophysics and scientific-computing communities by Pankaj K. Mishra and
Cédric Patzer, including at the EGU General Assembly 2026, the 27th
Electromagnetic Induction Workshop (EMIW 2026), and JuliaCon 2026.
MTGeophysics.jl was also used to participate in the
[MT3DINV-4 workshop](https://mt3dinv4.mtnet.info/Real_data.html)
(Memorial University of Newfoundland, 2025), a community 3D MT inversion
benchmarking exercise in which many groups inverted a common real field
dataset acquired over the Raglan mining district in northern Quebec, Canada.

# AI usage disclosure

GitHub Copilot was used as a coding assistant during development of this
software and in drafting this paper. All AI-generated code and text were
reviewed, tested, and verified by the author for correctness.

# Acknowledgements

This work was supported by the Research Council of Finland (project 359261).
The author wishes to acknowledge CSC – IT Center for Science, Finland, for
computational resources. The COMMEMI benchmark models used for validation were
defined by @Zhdanov1997.

# References

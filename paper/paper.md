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

MTGeophysics.jl is a Julia package for magnetotelluric (MT) geophysical
modelling, inversion, and interactive visualization. The magnetotelluric
method uses natural electromagnetic signals to image the electrical
resistivity structure of the Earth's subsurface, with applications ranging
from mineral exploration to studies of tectonic structure. Beyond providing
the usual building blocks (1D and 2D forward solvers, ModEM 3D I/O
[@Kelbert2014], and interactive 3D viewers), the package is organised around
a single guiding idea: making stochastic inversion and uncertainty
quantification practical for MT. Because many different subsurface
configurations can explain the same observations, a single "best" model is
rarely enough; MTGeophysics.jl is designed to move MT interpretation along
the spectrum from one deterministic model, to an ensemble of plausible
models, toward the full posterior distribution (\autoref{fig:vision}).

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

Standard 3D MT inversion returns a single resistivity model, typically the
minimiser of a data misfit plus a smoothness regulariser whose weight is
chosen by hand. This is computationally efficient but hides the
non-uniqueness that is intrinsic to MT: the recovered model depends on the
regularisation choice, and no uncertainty is reported. Quantifying that
uncertainty rigorously, by sampling many models consistent with the data,
is the natural remedy, but it has been considered impractical in 3D because
each forward solve is expensive and the model space is enormous.

MTGeophysics.jl exists to close this gap in practice rather than only in
principle. It provides a stochastic inversion workflow based on Very Fast
Simulated Annealing (VFSA) [@SenStoffa2013] that produces an *ensemble* of
plausible 3D models, from which ensemble statistics (mean, median, standard
deviation) give a per-voxel picture of what the data do and do not constrain.
The workflow reuses the community-standard ModEM forward
solver [@Egbert2012; @Kelbert2014] by reading and writing its native file
formats, so it slots into established MT practice without asking users to
change their models, data, or solver.

The package targets MT researchers and students who want a unified,
open-source toolkit that builds on Julia's [@Bezanson2017] strengths in
numerical computing, automatic differentiation, composability, and
interactive graphics. It is designed as a research repository and as a core
component of a broader JuliaGeophysics ecosystem: forward solvers, data
structures, and inversion routines are meant to be reused and recombined, so
that new ideas (alternative parameterisations, hybrid inversion strategies,
or novel uncertainty quantification schemes) can be prototyped without
rebuilding core MT tooling from scratch.

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
larger codebases. None of these offer integrated stochastic inversion with
ensemble uncertainty quantification in both 2D and 3D. MTGeophysics.jl
contributes a Julia-native package that ties together 1D/2D forward
modelling, 2D/3D VFSA inversion with ensemble uncertainty quantification, and
interactive 3D visualization with GIS overlay support.

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
(\autoref{fig:loop}). Inverting for a handful of coefficients instead of
millions of voxels is what makes ensemble exploration feasible at regional
scale.

![One VFSA iteration as a sparse-model update loop. The padded modelling mesh
(1) is reduced to its core (2), sampled at $M$ random control points, the
*sparse model* (3), perturbed by the VFSA rule (4), then mapped back to the
full padded mesh by compactly supported Gaussian-RBF interpolation (5) for the
ModEM forward solve, advancing $\mathbf{m}_k \to \mathbf{m}_{k+1}$. This
reduced parameterisation is the strategy that makes stochastic 3D inversion
tractable.\label{fig:loop}](vfsa_loop.png){ width=90% }

This reduced parameterisation is more than a VFSA implementation detail. The
same principle, searching a low-dimensional representation rather than the
full voxel space, is a prerequisite for practical probabilistic inversion.
Fully Bayesian approaches such as Markov chain Monte Carlo mirror the
multi-chain structure already used here, but they only become affordable in 3D
when the model is expressed through a similarly compact parameterisation. By
building the package around reduced parameterisation from the outset,
MTGeophysics.jl is positioned to extend naturally from ensemble-based VFSA
toward full posterior sampling.

Around this core, the package is organised as a single Julia module with
clearly separated layers: ModEM 3D data/model I/O; a 1D module with analytical
and finite-difference solvers; a 2D module that assembles sparse
finite-difference operators for the TE/TM Maxwell equations and solves them
with LinearSolve.jl; and the VFSA inversion module. In 2D the package includes
its own finite-difference forward engine; in 3D it wraps the external ModEM
solver, running multiple independent Markov chains in parallel and computing
ensemble statistics across them. Interactive visualization is handled through
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
Because GLMakie is a hard dependency that cannot precompile without a display,
the workflow provisions a virtual framebuffer (Xvfb) and runs the tests
headlessly, ensuring the visualization layer builds in a clean, display-free
environment. The suite combines smoke tests of the core 3D data structures and
ModEM I/O with regression checks of the 1D and 2D forward solvers, including
the COMMEMI benchmark generators, so numerical results are re-verified on every
change. Companion workflows build and deploy the documentation, keep dependency
bounds current, and, together with the plain-text file formats, keep the
repository reproducible and contribution-friendly.

# Benchmarks and validation

The package is validated at two levels, both fully reproducible from the
repository. The 2D forward solver and VFSA inversion reproduce the
community-standard COMMEMI 2D benchmark suite [@Zhdanov1997], generated
natively by `helpers/benchmarks_2d.jl` without external downloads. At regional
scale, the 3D VFSA inversion is benchmarked against USArray MT data from the
Cascadia subduction zone [@Patro2008], a dataset extensively studied with
the deterministic ModEM NLCG inversion [@Mishra2026]. The VFSA ensemble mean
recovers the major conductive structures found by deterministic inversion,
while the ensemble standard deviation provides a per-voxel uncertainty
estimate that no single deterministic inversion can offer; high-variance
regions flag where the data poorly constrain the model. Because the workflow
returns a *distribution* of plausible models rather than one solution, it
directly quantifies the non-uniqueness inherent in MT inversion
(\autoref{fig:cascadia}).

![Cascadia field benchmark. Top row: two map slices from the deterministic
inversion (ModEM NLCG). Middle row: the ensemble mean of five independent VFSA
chains, which recovers the main conductive structures of the deterministic
model. Bottom row: the ensemble standard deviation, a per-voxel uncertainty
map showing where the data constrain resistivity and where they do not. Each
non-greedy VFSA chain performs roughly 12,000 ModEM forward solves (about
13.5 s each, near two days per chain), with the five chains run in
parallel.\label{fig:cascadia}](VFSA3DBenchmark.png){ width=95% }

# Research impact statement

MTGeophysics.jl is developed to support ongoing magnetotelluric research at
the Geological Survey of Finland (GTK) and is intended to lower the barrier to
uncertainty-aware MT interpretation for the wider community. By making
ensemble-based stochastic inversion practical on top of the established ModEM
solver, the package lets researchers report where a resistivity model is well
constrained and where it is not, information no single deterministic
inversion provides. Its reduced-parameterisation design is a deliberate step
toward fully probabilistic MT inversion, and its reusable Julia components are
meant to serve as a foundation for the broader JuliaGeophysics ecosystem.

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

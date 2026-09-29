# Data dashboard

The dashboard pages through the sites of a ModEM data file, one site at a time. It is the place to
look at field data before an inversion, and at the fit after one.

## Open it

```julia
using MTGeophysics

DataDashboard("data.dat")                              # observed only
DataDashboard("data.dat"; predicted = "dpred.dat")     # observed against a ModEM response
```

EDI files are converted to a ModEM file first (see [EDI files](../data/files.md#EDI-files)):

```julia
DataDashboard(EDIToModEM("EDI/COPROD2"))               # writes EDI/COPROD2.dat
```

From the command line, with a directory of EDIs or a ModEM file:

```bash
julia --project=. examples/DataDashboard.jl <data.dat | edi_dir> [predicted.dat]
```

The window needs GLMakie and a display. On a cluster, `--png <dir>` (or
`interactive = false, snapshot_dir = "<dir>"`) writes one PNG per site with CairoMakie instead.

## What it shows

The top row is the apparent resistivity and phase of the four impedance components. Observed data are
circles with error bars, predicted data are lines, and masked points are hollow. The yx phase is
turned by 180° so both off-diagonal phases sit in the first quadrant. The **Diagonals** toggle shows or
hides Zxx and Zyy.

Below are four panels, each with its own menu:

| Panel | Content |
|:------|:--------|
| Tipper | Re and Im of Tzx and Tzy |
| Phase tensors & induction arrows | PT ellipses along period, coloured by Φmin, with the real (black) and imaginary (grey) induction arrows; one row each for observed and predicted |
| PT skew β | with the ±3° band |
| Ellipticity, Swift & Bahr skew | the PT ellipticity, Swift κ and Bahr η, with η = 0.3 marked |
| Strike | PT azimuth and Swift strike, clockwise from north, modulo 90° |
| Normalised residuals | (predicted − observed) / error for the real and imaginary parts, with the site RMS |
| Relative errors | \|δZ\|/\|Z\| per component and the tipper errors, with the error floor marked |
| Niblett–Bostick | Bostick resistivity against depth for xy, yx and the determinant |

The map on the left shows the sites, coloured by their RMS when a response is loaded, and the phase
tensor ellipses at one period (the slider underneath). Click a site on the map to go to it.

Residuals and RMS use the errors raised to `z_floor`·√|Zxy·Zyx| (default 5 %) and `t_floor` (0.03),
as an inversion would. The error bars show the recorded errors.

## Keys and editing

| Action | Effect |
|:-------|:-------|
| ← / → | previous / next site |
| ↑ / ↓ | map period |
| click a site on the map | go to that site |
| **Edit mask** on, click a ρ, φ or tipper point | mask it, or restore it |
| **Drop site** / **Reset site** | mask or restore the whole site |
| **Export ModEM** | write the kept data, with floored errors, to `<name>Edited.dat` (or `export_path`) |
| **Save PNG** | write the window as a PNG |

## Options

```julia
DataDashboard("data.dat";
    predicted = "dpred.dat",
    z_floor = 0.05, t_floor = 0.03,                   # errors for residuals, RMS and export
    panels = [:tipper, :ptstrip, :resid, :beta],      # the four panels at start
    iv_convention = :parkinson,                       # :parkinson towards conductors, :wiese away
    export_path = "edited.dat")
```

A response from a rotated mesh is turned back into the frame of the observed data before it is
compared. Sites are matched by name, periods within 2 %.

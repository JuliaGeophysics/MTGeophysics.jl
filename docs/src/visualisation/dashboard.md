# Data dashboard

The dashboard pages through the sites of a ModEM data file, one site at a time, and masks data. It is
the place to clean field data before an inversion, and to look at the fit after one. The mask is saved
as a text file that applies to the EDIs or to any ModEM file of the same survey.

## Open it

```julia
using MTGeophysics

DataDashboard("data.dat")                              # observed only
DataDashboard("data.dat"; predicted = "dpred.dat")     # observed against a ModEM response
DataDashboard("data.dat"; mask = "data_mask.txt")      # continue from a saved mask
```

EDI files are converted to a ModEM file first (see [EDI files](../data/files.md#EDI-files)):

```julia
DataDashboard(EDIToModEM("EDI/COPROD2"))               # writes EDI/COPROD2.dat
```

From the command line, with a directory of EDIs or a ModEM file:

```bash
julia --project=. examples/edit_data.jl <data.dat | edi_dir> [predicted.dat] [--mask mask.txt]
```

The window needs GLMakie and a display. On a cluster, `--png <dir>` (or
`interactive = false, snapshot_dir = "<dir>"`) writes one PNG per site with CairoMakie instead.

## What it shows

| Panel | Content |
|:------|:--------|
| Apparent resistivity | Zxy (red) and Zyx (blue), with error bars; the range fits the kept points with half a decade to spare |
| Phase | the same components as recorded, on a fixed −200° to 200°: Zxy in the first quadrant, Zyx in the third |
| Tipper | Re (red) and Im (blue) of Tzx above Tzy, on the right |
| Map | the sites in longitude and latitude (true distances), on the left; click one to go there |

**Full tensor** adds Zxx (pale green) and Zyy (pale purple) to the same two panels; it only changes the
view, and masking always takes all four impedances. The map and the
tipper collapse with the buttons at the left and right edges; the tipper starts collapsed when the survey
has none. With all three open the map, ρa/φ and tipper columns share the width 1 : 2 : 2; with one side
closed the other two are equal. Every site shares the survey's period axis. Observed data are markers, a
predicted response is lines, and masked points are hollow.

## Masking

No panel zooms or pans, and masking is by drag only; a click on a data panel does nothing:

| Action | Effect |
|:-------|:-------|
| drag across ρa or φ | select a band of periods over both panels and mask all four impedances in it |
| drag across a tipper panel | mask Tzx and Tzy in the band |
| drag with Shift held | restore the band instead |
| **Mask site** / **Restore site** | mask or restore the whole site |
| **Save mask** | write the mask to `mask_path` (default: the `mask` it started from, else `<name>_mask.txt`) |
| **Write ModEM** | write the kept data to `<name>_<date_time>.dat` |
| **Write EDI** | write the kept data as one EDI per site to `<name>_EDI_<date_time>/` |
| **Save PNG** | write the window as a PNG |
| ← / → | previous / next site |

Both writes keep the errors as they were read; the editor applies no error floor. A new time-stamped
name each time means no write replaces an earlier one.

A band on ρa or φ takes the full impedance tensor at its periods, and one on the tipper takes Re and Im of both components.
Axis ranges follow the kept points, so masking an outlier zooms in on the rest; a masked point that
falls outside comes back with **Restore site** or a restoring band.

## Mask files

A mask file lists what is left out, one line per site, period and component. `*` stands for every period
or every component:

```text
# MTGeophysics data mask, 2026-09-29 14:02:11, from Quantec2017.dat
# one masked datum per line; * is every period or every component; components Zxx Zxy Zyx Zyy Tzx Tzy
# zrot 0.00
# site  period_s  component
MT1247  8.192021e-02  *
MT1247  2.925700e+01  Zxy
MT1250  *  *
```

Apply it to the survey in either form; the originals are not touched:

```julia
apply_data_mask("data.dat", "data_mask.txt")   # -> data_masked_<date_time>.dat
apply_data_mask("EDI/", "data_mask.txt")       # -> EDI_masked_<date_time>/ and EDI_masked_<date_time>.dat
```

For a ModEM file the masked lines are dropped and everything else is kept as written. EDIs are read as
`EDIToModEM` reads them, written again without the masked values (fully masked sites are left out), and
converted to the ModEM file beside them. Sites match by name and periods within 1 % (`rtol`). The mask
records the frame it was drawn in (`zrot`), and applying it to data in another frame gives a warning,
since Zxy in one frame is not Zxy in another. In code, `apply_data_mask!(d, mask)` masks a `Data` in
place and `mask_keep(d, mask)` returns what a mask keeps.

## Options

```julia
DataDashboard("data.dat";
    predicted = "dpred.dat",
    mask = "data_mask.txt", mask_path = "data_mask.txt",
    show_map = true, show_tipper = nothing,           # side panels at start; nothing: open if there is a tipper
    full_tensor = false,
    export_dir = "exports")                           # where Write ModEM / Write EDI go (default: beside the data)
```

A response from a rotated mesh is turned back into the frame of the observed data before it is
compared. Sites are matched by name, periods within 2 %.

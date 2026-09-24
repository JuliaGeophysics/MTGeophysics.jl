# Data Rotation

`rotate_data` turns the impedance tensor and tipper of a ModEM data file and records why. The algebra is the same
for 1D, 2D and 3D: turning the axes by θ clockwise gives

```text
Z' = R Z Rᵀ,   T' = T Rᵀ,   R = [cos θ  sin θ; −sin θ  cos θ]
```

and the errors are propagated as independent variances. The 2D strike code (`fwd.ctrl` `Strike (deg)`) uses it.

```julia
rotate_data("data.obs", 30.0)                                  # mesh turn, writes data-r.obs
rotate_data("data-r.obs", 9.5; kind = :declination)            # rewrites data-r.obs
rotate_data("data.obs", [9.4, 9.5, 9.6]; kind = :declination)  # one declination per site
d = rotate_data(load_data_modem("data.obs"), 30.0)             # in memory
```

## The header angle

The rotation line of the ModEM header (`> 30.00`) is the azimuth of the data's x axis, degrees clockwise from
geographic north, summed over all rotations. ModEM reads and writes it back but computes in the grid frame, so it
is a label: the data must already be in the mesh's frame.

| `kind` | Fields | Station x, y | Header angle | Use |
|:-------|:-------|:-------------|:-------------|:----|
| `:mesh` (default) | turn by `angle` | turn by `angle` | + `angle` | match a rotated mesh |
| `:strike` | turn by `angle` | turn by `angle` | + `angle` | 2D strike frame |
| `:declination` | turn by −`angle` | unchanged | unchanged | data recorded with x on magnetic north |

A declination correction assumes the file was recorded with x on magnetic north but labelled geographic. Turning
the fields by −D (the declination, east positive) makes that label true, so the header angle and the geographic
station positions stay. The angle may differ per site; mesh and strike turns take one angle, since the frame
turns as a whole. Lat/lon never change.

## History

Every call appends a step to the first `#` line of each data block, after the description:

```text
# Rotated by MTGeophysics.jl | rotated: mesh +30.00; declination [JK01 +9.40, JK02 +9.50, JK03 +9.60]
```

`load_data_modem` and `load_data2d` read it back as `rotations` (`RotationStep`s), and every file the package
writes from those data carries it on. ModEM keeps only the first 200 characters of that line in its own output,
so a long per-site list may be cut there; the header angle is unaffected. Correcting twice for declination warns.

## Output and incomplete data

The rotated file is `<stem>-r<ext>` next to the input (`data.obs` gives `data-r.obs`); a file already named
`-r` is rewritten in place, since the history is in it. `output_path` names another file. The file keeps its
sign convention and units.

A site-period without the full tensor (or both tipper components) cannot be rotated and is dropped with a
warning. Data with no complete tensor at all, such as the ZXY/ZYX-only 2D synthetics, are an error.

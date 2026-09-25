# Rotation and strike

MT data are recorded in one coordinate frame, but a mesh or a 2D model often needs another. In this
section we rotate data files with `rotate_data`, and see how the 2D workflows find and use the
geoelectric strike.

## Rotate a data file

Turning the axes by an angle θ (clockwise) changes the impedance tensor Z and the tipper T as

```text
Z' = R Z Rᵀ,   T' = T Rᵀ,   R = [cos θ  sin θ; −sin θ  cos θ]
```

and the errors are carried along. This is the same for 1D, 2D and 3D data.

```julia
rotate_data("data.obs", 30.0)                                  # writes data-r.obs
rotate_data("data-r.obs", 9.5; kind = :declination)            # updates data-r.obs
rotate_data("data.obs", [9.4, 9.5, 9.6]; kind = :declination)  # one angle per site
d = rotate_data(load_data_modem("data.obs"), 30.0)             # in memory
```

The rotated file is written next to the input as `<stem>-r<ext>`, so `data.obs` becomes
`data-r.obs`. A file that already ends in `-r` is updated in place.

## Three kinds of rotation

| `kind` | What it does | Use it when |
|:-------|:-------------|:------------|
| `:mesh` (default) | turns the data and the station positions | your mesh is rotated |
| `:strike` | the same, for the 2D strike frame | a 2D workflow does this for you |
| `:declination` | turns the data by −D, keeps the station positions | the data were recorded with x on magnetic north |

The rotation line of the ModEM header (`> 30.00`) always holds the total angle of the data's x axis
from geographic north. Latitude and longitude never change.

!!! warning
    A site and period without the full impedance tensor (or both tipper components) cannot be
    rotated. It is dropped with a warning.

## Every rotation is recorded

Each call adds a step to the first `#` line of the data block, for example:

```text
# Rotated by MTGeophysics.jl | rotated: mesh +30.00; declination [JK01 +9.40, JK02 +9.50, JK03 +9.60]
```

The package reads this history back with the data and carries it into every file it writes, so you
can always see what was done to a data set.

## The 2D strike

The 2D solver assumes that x runs along the geological strike and y along the profile. Every 2D
workflow therefore first rotates the data to the strike given in `fwd.ctrl`:

```text
Strike (deg)             : auto
```

- `auto` (the default) estimates the strike from the phase tensors of the data.
- A number, such as `32.5`, sets it in degrees clockwise from north.

You can also estimate the strike or rotate to it yourself:

```julia
s = EstimateStrike2D(load_data2d("data.obs"))   # strike, with consistency and skew
RotateToStrike2D("data.obs")                    # writes data-r.obs
```

Inversions write the rotated data to their run folder and the strike to `Summary.txt`.

!!! note
    The mesh depends on the station positions in the strike frame. If you change the strike, run
    the [mesh tool](mesh2d.md) again.

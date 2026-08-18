using MTGeophysics

source_model = length(ARGS) >= 1 ? ARGS[1] : ""
source_data  = length(ARGS) >= 2 ? ARGS[2] : ""
target_mesh  = length(ARGS) >= 3 ? ARGS[3] : ""

outdir      = ""
stem        = ""
method      = :nearest
water_log10 = NaN
snap_sites  = true
smooth_xy   = 0.3
smooth_z    = 0.3
napply      = 1

MeshToMesh(source_model, source_data, target_mesh;
    outdir      = outdir,
    stem        = stem,
    method      = method,
    water_log10 = water_log10,
    snap_sites  = snap_sites,
    smooth_xy   = smooth_xy,
    smooth_z    = smooth_z,
    napply      = napply,
)

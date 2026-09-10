using MTGeophysics

flags      = filter(a -> startswith(a, "--"), ARGS)
positional = filter(a -> !startswith(a, "--"), ARGS)

data_file  = length(positional) >= 1 ? positional[1] : ""
target_crs = length(positional) >= 2 ? positional[2] : "EPSG:4326"

# extra overlays can also be imported from the viewer with "Add Shapefile"
shapefiles = []

shapefile_color      = :grey30
shapefile_alpha      = 0.9
shapefile_line_width = 1.2
shapefile_point_size = 6

pt_fill         = :beta        # :beta, :phimin or :phi2
pt_as_angle     = true
pt_scale        = 0.66
pt_scale_step   = 1.25
pt_colormap     = :Spectral
pt_range        = nothing
pt_stroke       = :black
pt_strokewidth  = 1.1
skip_beta_above = nothing

iv_convention    = :parkinson  # :parkinson points towards conductors, :wiese away
iv_scale         = 1.2
iv_scale_step    = 1.25
iv_real_color    = :black
iv_imag_color    = :grey45
iv_linewidth     = 1.4
iv_head_frac     = 0.30
iv_head_width    = 0.42
iv_max_magnitude = 1.0

show_phase_tensor = true
show_iv_real      = true
show_iv_imag      = false
show_sites        = true
show_site_labels  = false
site_color        = :grey20
site_markersize   = 4

map_pad        = 0.06
viewer_figsize = (1250, 950)
export_dpi     = 3
export_figsize = (1150, 950)

gis_output_dir = ""            # default <data>-PTIV-GIS in the working directory

# --gis-only writes the shapefiles and skips the viewer; it needs no display
gis_only = ("--gis-only" in flags) || ("--no-interactive" in flags)

if gis_only
    write_ptiv_gis(data_file;
        crs              = target_crs,
        output_dir       = gis_output_dir,
        as_angle         = pt_as_angle,
        pt_scale         = pt_scale,
        iv_scale         = iv_scale,
        iv_convention    = iv_convention,
        iv_max_magnitude = iv_max_magnitude,
        iv_head_frac     = iv_head_frac,
        iv_head_width    = iv_head_width,
        skip_beta_above  = skip_beta_above,
    )
    exit(0)
end

PlotPTIVMap(data_file;
    crs                  = target_crs,
    shapefiles           = shapefiles,
    shapefile_color      = shapefile_color,
    shapefile_alpha      = shapefile_alpha,
    shapefile_line_width = shapefile_line_width,
    shapefile_point_size = shapefile_point_size,
    pt_fill              = pt_fill,
    pt_as_angle          = pt_as_angle,
    pt_scale             = pt_scale,
    pt_scale_step        = pt_scale_step,
    pt_colormap          = pt_colormap,
    pt_range             = pt_range,
    pt_stroke            = pt_stroke,
    pt_strokewidth       = pt_strokewidth,
    skip_beta_above      = skip_beta_above,
    iv_convention        = iv_convention,
    iv_scale             = iv_scale,
    iv_scale_step        = iv_scale_step,
    iv_real_color        = iv_real_color,
    iv_imag_color        = iv_imag_color,
    iv_linewidth         = iv_linewidth,
    iv_head_frac         = iv_head_frac,
    iv_head_width        = iv_head_width,
    iv_max_magnitude     = iv_max_magnitude,
    show_phase_tensor    = show_phase_tensor,
    show_iv_real         = show_iv_real,
    show_iv_imag         = show_iv_imag,
    show_sites           = show_sites,
    show_site_labels     = show_site_labels,
    site_color           = site_color,
    site_markersize      = site_markersize,
    map_pad              = map_pad,
    viewer_figsize       = viewer_figsize,
    export_dpi           = export_dpi,
    export_figsize       = export_figsize,
    gis_output_dir       = gis_output_dir,
)

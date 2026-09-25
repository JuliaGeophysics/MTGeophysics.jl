using Documenter, DocumenterVitepress

makedocs(;
    sitename = "MTGeophysics.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = DocumenterVitepress.MarkdownVitepress(;
        repo       = "github.com/JuliaGeophysics/MTGeophysics.jl",
        devbranch  = "main",
        devurl     = "dev",
        # full URL with https://, otherwise the host is taken as part of the base path
        deploy_url = "https://juliageophysics.com/MTGeophysics.jl",
        description = "Magnetotelluric forward modelling, inversion and visualisation in Julia",
    ),
    # nested lists become the dropdown menus of the top navigation bar
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Data & Meshes" => [
            "Data and model files" => "data/files.md",
            "Rotation and strike"  => "data/rotation.md",
            "1D and 2D meshes"     => "data/mesh2d.md",
            "3D meshes"            => "data/mesh3d.md",
        ],
        "Forward" => [
            "1D" => "forward/1d.md",
            "2D" => "forward/2d.md",
            "3D" => "forward/3d.md",
        ],
        "Inversion" => [
            "1D"               => "inversion/1d.md",
            "2D deterministic" => "inversion/2d_deterministic.md",
            "2D VFSA"          => "inversion/2d_vfsa.md",
            "3D VFSA"          => "inversion/3d_vfsa.md",
        ],
        "Visualisation" => [
            "1D"              => "visualisation/1d.md",
            "2D"              => "visualisation/2d.md",
            "3D models"       => "visualisation/3d_models.md",
            "3D data maps"    => "visualisation/3d_data.md",
            "3D model editing" => "visualisation/editing.md",
        ],
        "API" => "api.md",
        "For Developers" => [
            "Code layout"          => "developers/index.md",
            "Control file reference" => "developers/control_files.md",
            "2D solver notes"      => "developers/solver2d.md",
            "2D inversion theory"  => "developers/inversion2d.md",
            "Adding an algorithm"  => "developers/new_algorithm.md",
        ],
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaGeophysics/MTGeophysics.jl.git",
    devbranch = "main",
    push_preview = true,
)

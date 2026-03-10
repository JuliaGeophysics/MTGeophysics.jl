using Documenter

makedocs(;
    sitename = "MTGeophysics.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://pankajkmishra.github.io/MTGeophysics.jl",
        assets     = ["assets/custom.css"],
        collapselevel = 1,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started"       => "getting_started.md",
        "1D Forward Modelling"  => "forward1d.md",
        "2D Forward Modelling"  => "forward2d.md",
        "2D VFSA Inversion"     => "inversion2d.md",
        "3D Visualization"      => "visualization3d.md",
        "API Reference"         => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/pankajkmishra/MTGeophysics.jl.git",
    devbranch = "main",
    push_preview = true,
)

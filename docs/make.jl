using Documenter

makedocs(;
    sitename = "MTGeophysics.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://juliageophysics.github.io/MTGeophysics.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started"       => "getting_started.md",
        "1D Forward Modelling"  => "forward1d.md",
        "2D Forward Modelling"  => "forward2d.md",
        "2D VFSA Inversion"     => "inversion2d.md",
        "3D VFSA Inversion"     => "inversion3d.md",
        "3D Visualization"      => "visualization3d.md",
        "Model Editing"         => "editing.md",
        "API Reference"         => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaGeophysics/MTGeophysics.jl.git",
    devbranch = "main",
    push_preview = true,
)

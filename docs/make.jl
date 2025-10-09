using Documenter

# Try to load MTGeophysics, but don't fail if dependencies are missing  
modules_list = Module[]
try
    using Pkg
    Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
    using MTGeophysics
    global modules_list = [MTGeophysics]
catch e
    @warn "Could not load MTGeophysics, documentation will be built without docstrings" exception=e
end

makedocs(;
    modules=modules_list,
    authors="JuliaGeophysics community, Pankaj K Mishra and contributors",
    repo="https://github.com/JuliaGeophysics/MTGeophysics.jl/blob/{commit}{path}#{line}",
    sitename="MTGeophysics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGeophysics.github.io/MTGeophysics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Tutorials" => "tutorials.md",
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references, :eval_block],
)

deploydocs(;
    repo="github.com/JuliaGeophysics/MTGeophysics.jl",
    devbranch="main",
    push_preview=true,
)

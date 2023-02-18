using SubsidenceChron
using Documenter

DocMeta.setdocmeta!(SubsidenceChron, :DocTestSetup, :(using SubsidenceChron); recursive=true)

makedocs(;
    modules=[SubsidenceChron],
    authors="C. Brenhin Keller",
    repo="https://github.com/JuliaGeochronology/SubsidenceChron.jl/blob/{commit}{path}#{line}",
    sitename="SubsidenceChron.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliageochronology.github.io/SubsidenceChron.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeochronology/SubsidenceChron.jl",
    devbranch = "main",
)

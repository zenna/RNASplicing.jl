using Documenter, RNASplicing

makedocs(;
    modules=[RNASplicing],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/zenna/RNASplicing.jl/blob/{commit}{path}#L{line}",
    sitename="RNASplicing.jl",
    authors="Zenna Tavares",
    assets=String[],
)

deploydocs(;
    repo="github.com/zenna/RNASplicing.jl",
)

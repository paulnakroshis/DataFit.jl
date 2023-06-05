using DataFit
using Documenter

DocMeta.setdocmeta!(DataFit, :DocTestSetup, :(using DataFit); recursive=true)

makedocs(;
    modules=[DataFit],
    authors="Paul Nakroshis <pauln@maine.edu>",
    repo="https://github.com/paulnakroshis/DataFit.jl/blob/{commit}{path}#{line}",
    sitename="DataFit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://paulnakroshis.github.io/DataFit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/paulnakroshis/DataFit.jl",
    devbranch="main",
)

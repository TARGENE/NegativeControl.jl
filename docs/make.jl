using NegativeControl
using Documenter

DocMeta.setdocmeta!(NegativeControl, :DocTestSetup, :(using NegativeControl); recursive=true)

makedocs(;
    modules=[NegativeControl],
    authors="Olivier Labayle",
    repo="https://github.com/olivierlabayle/NegativeControl.jl/blob/{commit}{path}#{line}",
    sitename="NegativeControl.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olivierlabayle.github.io/NegativeControl.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olivierlabayle/NegativeControl.jl",
    devbranch="main",
)

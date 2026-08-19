using Documenter, SymplecticPauli

DocMeta.setdocmeta!(SymplecticPauli, :DocTestSetup, :(using SymplecticPauli); recursive=true)

makedocs(
    sitename="SymplecticPauli.jl",
    modules=[SymplecticPauli],
    authors="Omar Alsheikh",
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Pauli strings" => "strings.md",
            "Operators" => "operators.md",
            "Rotations" => "rotations.md",
        ],
        "API reference" => "api.md",
    ],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://ooalshei.github.io/SymplecticPauli.jl",
        repolink="https://github.com/ooalshei/SymplecticPauli.jl",
    ),
)

deploydocs(repo="https://github.com/ooalshei/SymplecticPauli.jl", devbranch="main")

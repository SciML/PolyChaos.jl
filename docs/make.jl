# using Pkg
# Pkg.add("Documenter")
using Documenter, PolyChaos

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(
    sitename = "PolyChaos.jl",
    clean = true,
    doctest = true,
    checkdocs = :exports,
    linkcheck = true,
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/PolyChaos/stable/"
    ),
    modules = [PolyChaos],
    authors = "tillmann.muehlpfordt@kit.edu",
    pages = pages
)

deploydocs(
    repo = "github.com/SciML/PolyChaos.jl.git";
    push_preview = true
)

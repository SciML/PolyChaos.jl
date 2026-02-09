# using Pkg
# Pkg.add("Documenter")
using Documenter, PolyChaos

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(
    sitename = "PolyChaos.jl",
    clean = true, doctest = false, linkcheck = true,
    linkcheck_ignore = [
        "https://www.sciencedirect.com/science/article/pii/S235246771830105X",
        "https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html",
    ],
    warnonly = [:docs_block, :missing_docs],
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

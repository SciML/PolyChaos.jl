# using Pkg
# Pkg.add("Documenter")
using Documenter, PolyChaos

include("pages.jl")

makedocs(
    sitename = "PolyChaos.jl",
    format = Documenter.HTML(analytics = "UA-90474609-3",
                         assets = ["assets/favicon.ico"],
                         canonical="https://polychaos.sciml.ai/stable/"),
    modules = [PolyChaos],
    authors = "tillmann.muehlpfordt@kit.edu",
    doctest = false,
    pages = pages
)


deploydocs(
    repo = "github.com/SciML/PolyChaos.jl.git",
    target = "build",
)

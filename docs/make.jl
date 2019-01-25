using Documenter
# push!(LOAD_PATH,"C:\\Users\\la5373\\Documents\\Code\\pcejl\\src")
# push!(LOAD_PATH,"..\\src\\")
using PolyChaos
using DocumenterLaTeX

makedocs(
    sitename = "PolyChaos",
    format = Documenter.HTML(),
    assets = ["assets/myfont.css"],
    # format = LaTeX(),
    modules = [PolyChaos],
    authors = "Tillmann Muehlpfordt",
    repo = "https://github.com/timueh/PolyChaos",
    doctest = true,
    pages = Any[
        "index.md",
        "TypeHierarchy.md",
        "Usage" => [
            "Numerical Integration"=>"NumericalIntegration.md",
            "Monic Orthogonal Polynomials"=>"OrthogonalPolynomials_canonical.md",
            "Scalar Products" => "ScalarProducts.md",
            "Polynomial Chaos" => [ "Basic Usage" => "PCEtutorial.md",
                                    "Chi Squared, One DOF" => "ChiSquared_k1.md",
                                    "Chi Squared, Several DOFs" => "ChiSquared_kGreater1.md",
                                    "Random ODE" => "RandomODE.md"]
            ],
    "math.md",
    "functions.md"
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

deploydocs(
    repo = "https://github.com/timueh/PolyChaos.jl.git",
)

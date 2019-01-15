using Documenter
push!(LOAD_PATH,"C:\\Users\\la5373\\Documents\\Code\\pcejl\\src")
push!(LOAD_PATH,"..\\src\\")
using PolyChaos

makedocs(
    sitename = "PolyChaos",
    Documenter.HTML(prettyurls=false),
    modules = [PolyChaos],
    pages = Any[
        "index.md",
        "Tutorials" => [
            "Numerical Integration"=>"NumericalIntegration.md",
            "Monic Orthogonal Polynomials"=>"OrthogonalPolynomials_canonical.md",
            "Polynomial Chaos" => [ "Basic Usage" => "PCEtutorial.md",
                                    "Chi Squared, One DOF" => "ChiSquared_k1.md",
                                    "Chi Squared, Several DOFs" => "ChiSquared_kGreater1.md",
                                    "Random ODE" => "RandomODE.md"]
            ],
    "math.md",
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

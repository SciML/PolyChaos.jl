# using Pkg
# Pkg.add("Documenter")
using Documenter, PolyChaos

makedocs(
    sitename = "PolyChaos.jl",
    format = Documenter.HTML(),
    assets = ["assets/myfont.css"],
    # format = LaTeX(),
    modules = [PolyChaos.jl],
    authors = "tillmann.muehlpfordt@kit.edu",
    # repo = "github.com/timueh/PolyChaos.jl.git",
    doctest = true,
    pages = Any[
        "index.md",
        "type_hierarchy.md",
        "Usage" => [
            "Numerical Integration"=>"numerical_integration.md",
            "Monic Orthogonal Polynomials"=>"orthogonal_polynomials_canonical.md",
            "Scalar Products" => "scalar_products.md",
            "Polynomial Chaos" => [ "Basic Usage" => "pce_tutorial.md",
                                    "Chi Squared, One DOF" => "chi_squared_k1.md",
                                    "Chi Squared, Several DOFs" => "chi_squared_k_greater1.md",
                                    "Random ODE" => "random_ode.md"]
            ],
    "math.md",
    "functions.md"
    ]
)


deploydocs(
    repo = "github.com/timueh/PolyChaos.jl.git",
    target = "build",
)

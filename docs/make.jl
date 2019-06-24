# using Pkg
# Pkg.add("Documenter")
using Documenter, PolyChaos

makedocs(
    sitename = "PolyChaos.jl",
    format = Documenter.HTML(assets = ["assets/myfont.css"]),
    modules = [PolyChaos],
    authors = "tillmann.muehlpfordt@kit.edu",
    doctest = true,
    pages = Any[
        "index.md",
        "type_hierarchy.md",
        "Usage" => [
            "numerical_integration.md",
            "quadrature_rules.md",
            "orthogonal_polynomials_canonical.md",
            "multiple_discretization.md",
            "scalar_products.md",
            "Polynomial Chaos" => [ "Basic Usage" => "pce_tutorial.md",
                                    "Chi Squared, One DOF" => "chi_squared_k1.md",
                                    "Chi Squared, Several DOFs" => "chi_squared_k_greater1.md",
                                    "Random ODE" => "random_ode.md"
                                  ],
            "Optimal Power Flow" => "DCsOPF.md"
            ],
    "math.md",
    "functions.md"
    ]
)


deploydocs(
    repo = "github.com/timueh/PolyChaos.jl.git",
    target = "build",
)

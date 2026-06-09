using Pkg

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    include("constructors.jl")
    include("quadrature.jl")
    include("recurrence_coefficients.jl")
    include("discretization.jl")
    include("logistic.jl")
    include("comparison_gaussquadrature.jl")
    include("quadrature_rules.jl")
    include("show.jl")
    include("types.jl")
    include("auxfuns.jl")
    include("densities.jl")
    include("polynomial_chaos.jl")
    include("inverse_transform_sampling.jl")
    include("multi_indices.jl")
    include("scalarproducts.jl")
    include("evaluate.jl")
end

if GROUP == "All" || GROUP == "QA"
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(path = joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    include("qa.jl")
end

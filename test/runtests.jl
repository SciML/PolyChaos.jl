using SciMLTesting

run_tests(;
    core = () -> begin
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
    end,
    # qa= runs under both "QA" and "All" (matches original GROUP == "All" || GROUP == "QA")
    qa = (; env = joinpath(@__DIR__, "qa"), body = joinpath(@__DIR__, "qa.jl")),
)

using Pkg
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    @safetestset "Constructors" begin include("constructors.jl") end
    @safetestset "Quadrature" begin include("quadrature.jl") end
    @safetestset "Recurrence Coefficients" begin include("recurrence_coefficients.jl") end
    @safetestset "Discretization" begin include("discretization.jl") end
    @safetestset "Logistic" begin include("logistic.jl") end
    @safetestset "Comparison Gauss Quadrature" begin include("comparison_gaussquadrature.jl") end
    @safetestset "Quadrature Rules" begin include("quadrature_rules.jl") end
    @safetestset "Show" begin include("show.jl") end
    @safetestset "Types" begin include("types.jl") end
    @safetestset "Auxiliary Functions" begin include("auxfuns.jl") end
    @safetestset "Densities" begin include("densities.jl") end
    @safetestset "Polynomial Chaos" begin include("polynomial_chaos.jl") end
    @safetestset "Inverse Transform Sampling" begin include("inverse_transform_sampling.jl") end
    @safetestset "Multi Indices" begin include("multi_indices.jl") end
    @safetestset "Scalar Products" begin include("scalarproducts.jl") end
    @safetestset "Evaluate" begin include("evaluate.jl") end
end

if GROUP == "All" || GROUP == "QA"
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(path = joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    include("qa.jl")
end

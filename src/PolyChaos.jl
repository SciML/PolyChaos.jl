__precompile__()

"""
    PolyChaos

Orthogonal-polynomial, quadrature, and polynomial-chaos tools for uncertainty
quantification. Construct a canonical basis such as [`LegendreOrthoPoly`](@ref),
evaluate it with [`evaluate`](@ref), and use its quadrature rule through
[`integrate`](@ref).
"""
module PolyChaos

    using SpecialFunctions: beta, gamma
    using SparseArrays: SparseVector, spzeros
    using Distributions: Beta, Continuous, Distribution, Gamma, Logistic, Normal,
        Uniform, Univariate
    import LinearAlgebra: I, dot
    import FFTW: ifft
    import Combinatorics: with_replacement_combinations
    import Base: show
    import AdaptiveRejectionSampling: RejectionSampler, run_sampler!
    import Statistics
    import GaussQuadrature: special_eigenproblem!

    include("typesAbstract.jl")
    include("typesMeasures.jl")
    include("typesQuad.jl")
    include("typesOrthoPolys.jl")
    include("typesTensor.jl")
    include("multi_indices.jl")
    include("evaluate.jl")
    include("quadrature_rules.jl")
    include("recurrence_coefficients_monic.jl")
    include("stieltjes.jl")
    include("densities.jl")
    include("auxfuns.jl")
    include("scalar_product.jl")
    include("tensor.jl")
    include("show.jl")
    include("sample.jl")
    include("polynomial_chaos.jl")

end

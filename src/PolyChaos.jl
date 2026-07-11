__precompile__()

module PolyChaos

    using SpecialFunctions: beta, gamma
    using SparseArrays: SparseVector, spzeros
    using Distributions: Beta, Continuous, Distribution, Gamma, Logistic, Normal,
        Uniform, Univariate
    import LinearAlgebra: I, dot, issymmetric
    import FFTW: ifft
    import Combinatorics: with_replacement_combinations
    import Base: show
    import AdaptiveRejectionSampling: RejectionSampler, run_sampler!
    import Statistics: mean, std, var
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

    @doc """
        AbstractMeasure

    Abstract supertype for measures used to define orthogonality in PolyChaos.
    Concrete subtypes provide a weight function, domain, and symmetry information.
    """ AbstractMeasure

    @doc """
        AbstractCanonicalMeasure

    Abstract supertype for built-in canonical measures with known recurrence
    coefficients or quadrature rules.
    """ AbstractCanonicalMeasure

    @doc """
        AbstractQuad{T}

    Abstract supertype for quadrature rules with nodes and weights of real element
    type `T`.
    """ AbstractQuad

    @doc """
        AbstractOrthoPoly{M, Q}

    Abstract supertype for orthogonal polynomial bases associated with measure type
    `M` and quadrature type `Q`.
    """ AbstractOrthoPoly

    @doc """
        AbstractCanonicalOrthoPoly

    Abstract supertype for orthogonal polynomial bases generated from canonical
    measures.
    """ AbstractCanonicalOrthoPoly

    @doc """
        AbstractTensor{T}

    Abstract supertype for tensors of scalar products for orthogonal polynomial basis
    type `T`.
    """ AbstractTensor

    @doc """
        Measure(name, w, dom, symmetric, pars = Dict())

    Create a custom measure with weight function `w`, support `dom`, symmetry flag
    `symmetric`, and optional parameter dictionary `pars`.
    """ Measure

    @doc """
        ProductMeasure(w, measures)

    Product measure built from univariate component measures and a product weight
    function.
    """ ProductMeasure

    @doc """
        LegendreMeasure()

    Canonical Legendre measure on `(-1, 1)` with unit weight.
    """ LegendreMeasure

    @doc """
        JacobiMeasure(a, b)

    Canonical Jacobi measure on `(-1, 1)` with shape parameters `a` and `b`.
    """ JacobiMeasure

    @doc """
        LaguerreMeasure()

    Canonical Laguerre measure on `(0, Inf)` with exponential weight.
    """ LaguerreMeasure

    @doc """
        genLaguerreMeasure(shape)

    Generalized Laguerre measure on `(0, Inf)` with shape parameter `shape`.
    """ genLaguerreMeasure

    @doc """
        HermiteMeasure()

    Canonical Hermite measure on the real line.
    """ HermiteMeasure

    @doc """
        genHermiteMeasure(mu)

    Generalized Hermite measure on the real line with parameter `mu`.
    """ genHermiteMeasure

    @doc """
        MeixnerPollaczekMeasure(lambda, phi)

    Meixner-Pollaczek measure on the real line with parameters `lambda` and `phi`.
    """ MeixnerPollaczekMeasure

    @doc """
        GaussMeasure()

    Standard Gaussian measure on the real line.
    """ GaussMeasure

    @doc """
        Uniform01Measure()

    Uniform measure on `(0, 1)`.
    """ Uniform01Measure

    @doc """
        Uniform_11Measure()

    Uniform measure on `(-1, 1)`.
    """ Uniform_11Measure

    @doc """
        Beta01Measure(a, b)

    Beta measure on `(0, 1)` with shape parameters `a` and `b`.
    """ Beta01Measure

    @doc """
        GammaMeasure(shape)

    Gamma measure on `(0, Inf)` with shape parameter `shape`.
    """ GammaMeasure

    @doc """
        LogisticMeasure()

    Logistic probability measure on the real line.
    """ LogisticMeasure

    @doc """
        Quad(name, N, nodes, weights)

    Quadrature rule named `name` with `N` nodes and corresponding weights.
    """ Quad

    @doc """
        EmptyQuad()

    Empty quadrature placeholder for orthogonal polynomial bases without an attached
    quadrature rule.
    """ EmptyQuad

    @doc """
        InconsistencyError(message)

    Exception used when basis dimensions, coefficients, or quadrature inputs are
    internally inconsistent.
    """ InconsistencyError

    @doc """
        OrthoPoly(args...; kwargs...)

    Construct a univariate orthogonal polynomial basis from recurrence coefficients,
    a measure, or a custom weight function.
    """ OrthoPoly

    @doc """
        LegendreOrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the Legendre measure up to degree `deg`.
    """ LegendreOrthoPoly

    @doc """
        JacobiOrthoPoly(deg, a, b; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the Jacobi measure with shape parameters `a` and
    `b`.
    """ JacobiOrthoPoly

    @doc """
        LaguerreOrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the Laguerre measure up to degree `deg`.
    """ LaguerreOrthoPoly

    @doc """
        genLaguerreOrthoPoly(deg, shape; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the generalized Laguerre measure.
    """ genLaguerreOrthoPoly

    @doc """
        HermiteOrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the Hermite measure up to degree `deg`.
    """ HermiteOrthoPoly

    @doc """
        genHermiteOrthoPoly(deg, mu; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the generalized Hermite measure.
    """ genHermiteOrthoPoly

    @doc """
        MeixnerPollaczekOrthoPoly(deg, lambda, phi; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the Meixner-Pollaczek measure.
    """ MeixnerPollaczekOrthoPoly

    @doc """
        GaussOrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the standard Gaussian measure.
    """ GaussOrthoPoly

    @doc """
        Uniform01OrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the uniform measure on `(0, 1)`.
    """ Uniform01OrthoPoly

    @doc """
        Uniform_11OrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the uniform measure on `(-1, 1)`.
    """ Uniform_11OrthoPoly

    @doc """
        Beta01OrthoPoly(deg, a, b; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the beta measure on `(0, 1)`.
    """ Beta01OrthoPoly

    @doc """
        GammaOrthoPoly(deg, shape; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the gamma measure.
    """ GammaOrthoPoly

    @doc """
        LogisticOrthoPoly(deg; Nrec = deg + 1, addQuadrature = true)

    Orthogonal polynomial basis for the logistic measure.
    """ LogisticOrthoPoly

    @doc """
        MultiOrthoPoly(args...; kwargs...)

    Multivariate orthogonal polynomial basis assembled from univariate bases and a
    multi-index set.
    """ MultiOrthoPoly

    @doc """
        Tensor(dim, basis)

    Tensor of scalar products of order `dim` for an orthogonal polynomial basis.
    """ Tensor

    @doc """
        dim(op)

    Return the number of basis polynomials represented by an orthogonal polynomial
    basis.
    """ dim

    @doc """
        deg(op)

    Return the maximum polynomial degree represented by an orthogonal polynomial
    basis.
    """ deg

    @doc """
        calculateMultiIndices(d, n)

    Construct the total-degree multi-index matrix for `d` variables up to degree `n`.
    """ calculateMultiIndices

    @doc """
        assign2multi(x, i, ind)

    Embed univariate coefficients `x` for variable `i` into the multivariate
    coefficient vector described by multi-index matrix `ind`.
    """ assign2multi

    @doc """
        multi2uni(a, ind)

    Convert scalar-product basis indices `a` into the corresponding univariate
    multi-indices from `ind`.
    """ multi2uni

    @doc """
        getentry(a, T, ind, dim)

    Return the sparse tensor entry indexed by basis tuple `a`.
    """ getentry

    @doc """
        computeTensorizedSP(m, args...)

    Compute the tensorized scalar products of order `m` for a univariate or
    multivariate orthogonal polynomial basis.
    """ computeTensorizedSP

    @doc """
        golubwelsch(alpha, beta, maxiter = 30)

    Compute Gaussian quadrature nodes and weights from recurrence coefficients using
    the Golub-Welsch eigenproblem.
    """ golubwelsch

    @doc """
        rm_chebyshev1(N)

    Return monic recurrence coefficients for first-kind Chebyshev polynomials.
    """ rm_chebyshev1

    @doc """
        w_legendre(t)

    Legendre weight function on `(-1, 1)`.
    """ w_legendre

    @doc """
        build_w_jacobi(a, b)

    Return a Jacobi weight function closure with parameters `a` and `b`.
    """ build_w_jacobi

    @doc """
        w_jacobi(t, a, b)

    Jacobi weight function on `(-1, 1)` with parameters `a` and `b`.
    """ w_jacobi

    @doc """
        w_hermite(t)

    Hermite weight function on the real line.
    """ w_hermite

    @doc """
        build_w_genhermite(mu)

    Return a generalized Hermite weight function closure with parameter `mu`.
    """ build_w_genhermite

    @doc """
        w_genhermite(t, mu)

    Generalized Hermite weight function with parameter `mu`.
    """ w_genhermite

    @doc """
        build_w_genlaguerre(a)

    Return a generalized Laguerre weight function closure with parameter `a`.
    """ build_w_genlaguerre

    @doc """
        w_laguerre(t)

    Laguerre weight function on `(0, Inf)`.
    """ w_laguerre

    @doc """
        w_meixner_pollaczek(t, lambda, phi)

    Meixner-Pollaczek weight function with parameters `lambda` and `phi`.
    """ w_meixner_pollaczek

    @doc """
        build_w_meixner_pollaczek(lambda, phi)

    Return a Meixner-Pollaczek weight function closure.
    """ build_w_meixner_pollaczek

    @doc """
        w_gaussian(t)

    Standard Gaussian probability density.
    """ w_gaussian

    @doc """
        build_w_beta(alpha, beta)

    Return a beta weight function closure with shape parameters `alpha` and `beta`.
    """ build_w_beta

    @doc """
        build_w_gamma(alpha)

    Return a gamma weight function closure with shape parameter `alpha`.
    """ build_w_gamma

    @doc """
        w_uniform01(t)

    Uniform density on `(0, 1)`.
    """ w_uniform01

    @doc """
        w_uniform_11(t)

    Uniform density on `(-1, 1)`.
    """ w_uniform_11

    @doc """
        w_logistic(t)

    Logistic probability density.
    """ w_logistic

end

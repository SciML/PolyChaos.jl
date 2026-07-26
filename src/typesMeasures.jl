export Measure,
    ProductMeasure,
    LegendreMeasure,
    JacobiMeasure,
    LaguerreMeasure,
    genLaguerreMeasure,
    HermiteMeasure,
    genHermiteMeasure,
    GaussMeasure,
    Beta01Measure,
    GammaMeasure,
    Uniform01Measure,
    Uniform_11Measure,
    LogisticMeasure,
    MeixnerPollaczekMeasure

Type_for_domain = Tuple{Float64, Float64}
Type_for_function = Tuple{<:Real, <:Real}

"""
    Measure(name, w, dom, symmetric, pars = Dict())

Construct a user-defined univariate measure.

# Arguments

- `name`: descriptive measure name; it is stored in lowercase.
- `w`: nonnegative weight function on `dom`.
- `dom`: ordered finite or infinite support bounds.
- `symmetric`: whether the measure is symmetric about zero.
- `pars`: optional parameter dictionary retained with the measure.

# Fields

- `name`: lowercase measure name.
- `w`: weight function.
- `dom`: support bounds.
- `symmetric`: symmetry flag used by scalar-product algorithms.
- `pars`: user-supplied metadata.

# Examples

```jldoctest
julia> using PolyChaos

julia> Measure("uniform", x -> 1.0, (0.0, 1.0), false).name
"uniform"
```
"""
struct Measure <: AbstractMeasure
    name::String
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    pars::Dict
    function Measure(
            name::String, w::Function, dom::Type_for_function, symm::Bool,
            d::Dict = Dict()
        )
        !(dom[1] < dom[2]) && throw(DomainError(dom, "invalid domain bounds specified"))
        return new(lowercase(name), w, dom, symm, d)
    end
end

"""
    ProductMeasure(w, measures)

Represent a product measure assembled from univariate component measures.

# Arguments

- `w`: product weight function accepting one point per component.
- `measures`: component [`AbstractMeasure`](@ref)s in coordinate order.

# Fields

- `w`: product weight function.
- `measures`: univariate component measures.
"""
struct ProductMeasure <: AbstractMeasure
    w::Function
    measures::Vector{<:AbstractMeasure}
end

# constructor for classic distributions
"""
    LegendreMeasure()

Canonical Legendre measure with unit weight on `(-1, 1)`.

# Fields

- `w`: the Legendre weight function.
- `dom`: `(-1.0, 1.0)`.
- `symmetric`: `true`.

# Examples

```jldoctest
julia> using PolyChaos

julia> LegendreMeasure().dom
(-1.0, 1.0)
```
"""
struct LegendreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    function LegendreMeasure()
        return new(w_legendre, (-1.0, 1.0), true)
    end
end

"""
    JacobiMeasure(shape_a, shape_b)

Canonical Jacobi measure on `(-1, 1)`.

# Arguments

- `shape_a`: exponent of `(1 - t)`; must exceed `-1`.
- `shape_b`: exponent of `(1 + t)`; must exceed `-1`.

# Fields

- `w`, `dom`, `symmetric`: measure data.
- `ashapeParameter`, `bshapeParameter`: validated shape parameters.
"""
struct JacobiMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    ashapeParameter::Float64
    bshapeParameter::Float64

    function JacobiMeasure(shape_a::Real, shape_b::Real)
        any
        shape_a <= -1 && throw(DomainError(shape_a, "shape parameter a must be > -1"))
        shape_b <= -1 && throw(DomainError(shape_b, "shape parameter b must be > -1"))
        return new(
            build_w_jacobi(shape_a, shape_b), (-1, 1), isapprox(shape_a, shape_b), shape_a,
            shape_b
        )
    end
end

"""
    LaguerreMeasure()

Canonical Laguerre measure with weight `exp(-t)` on `(0, Inf)`.

# Fields

- `w`: Laguerre weight.
- `dom`: `(0, Inf)`.
- `symmetric`: `false`.
"""
struct LaguerreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool

    function LaguerreMeasure()
        return new(w_laguerre, (0, Inf), false)
    end
end

"""
    genLaguerreMeasure(shape)

Generalized Laguerre measure with weight `t^shape * exp(-t)` on `(0, Inf)`.

# Arguments

- `shape`: exponent of `t`; must exceed `-1`.

# Fields

- `w`, `dom`, `symmetric`: measure data.
- `shapeParameter`: validated generalized-Laguerre shape.
"""
struct genLaguerreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    shapeParameter::Float64

    function genLaguerreMeasure(shape::Real)
        shape <= -1 && throw(DomainError(shape, "invalid shape parameter"))
        return new(build_w_genlaguerre(shape), (0, Inf), false, shape)
    end
end

"""
    HermiteMeasure()

Canonical Hermite measure with weight `exp(-t^2)` on the real line.

# Fields

- `w`: Hermite weight.
- `dom`: `(-Inf, Inf)`.
- `symmetric`: `true`.
"""
struct HermiteMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool

    function HermiteMeasure()
        return new(w_hermite, (-Inf, Inf), true)
    end
end

"""
    genHermiteMeasure(mu)

Generalized Hermite measure with weight `abs(t)^(2mu) * exp(-t^2)`.

# Arguments

- `mu`: generalized-Hermite parameter; must exceed `-0.5`.

# Fields

- `w`, `dom`, `symmetric`: measure data.
- `muParameter`: validated parameter.
"""
struct genHermiteMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    muParameter::Float64

    function genHermiteMeasure(mu::Real)
        mu <= -0.5 && throw(DomainError(mu, "invalid parameter value (must be > - 0.5)"))
        return new(build_w_genhermite(mu), (-Inf, Inf), true, mu)
    end
end

"""
    MeixnerPollaczekMeasure(lambda, phi)

Meixner-Pollaczek measure on the real line.

# Arguments

- `lambda`: positive shape parameter.
- `phi`: angle parameter in `(0, pi)`.

# Fields

- `w`, `dom`, `symmetric`: measure data.
- `λParameter`, `ϕParameter`: validated family parameters.
"""
struct MeixnerPollaczekMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    λParameter::Float64
    ϕParameter::Float64

    function MeixnerPollaczekMeasure(λ::Real, ϕ::Real)
        λ <= 0 && throw(DomainError(λ, "λ has to be positive"))
        !(0 < ϕ < pi) && throw(DomainError(ϕ, "ϕ has to be between 0 and pi"))
        return new(build_w_meixner_pollaczek(λ, ϕ), (-Inf, Inf), false, λ, ϕ)
    end
end

"""
    GaussMeasure()

Standard Gaussian probability measure on the real line.

# Fields

- `w`: standard-normal density.
- `dom`: `(-Inf, Inf)`.
- `symmetric`: `true`.
"""
struct GaussMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool

    function GaussMeasure()
        return new(w_gaussian, (-Inf, Inf), true)
    end
end

"""
    Uniform01Measure()

Uniform probability measure on `(0, 1)`.

# Fields

- `w`: unit density on `(0, 1)`.
- `dom`: `(0, 1)`.
- `symmetric`: `true` for the centered polynomial family.
"""
struct Uniform01Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool

    function Uniform01Measure()
        return new(w_uniform01, (0, 1), true)
    end
end

"""
    Uniform_11Measure()

Uniform probability measure on `(-1, 1)`.

# Fields

- `w`: density `1 / 2`.
- `dom`: `(-1, 1)`.
- `symmetric`: `true`.
"""
struct Uniform_11Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool

    function Uniform_11Measure()
        return new(w_uniform_11, (-1, 1), true)
    end
end

"""
    Beta01Measure(a, b)

Beta probability measure on `(0, 1)`.

# Arguments

- `a`: first positive beta shape parameter.
- `b`: second positive beta shape parameter.

# Fields

- `w`, `dom`, `symmetric`: measure data.
- `ashapeParameter`, `bshapeParameter`: validated beta parameters.
"""
struct Beta01Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    ashapeParameter::Float64
    bshapeParameter::Float64

    function Beta01Measure(a, b)
        a <= 0 && throw(DomainError(a, "shape parameter a must be positive"))
        b <= 0 && throw(DomainError(b, "shape parameter b must be positive"))
        return new(build_w_beta(a, b), (0, 1), isapprox(a, b), a, b)
    end
end

"""
    GammaMeasure(shape, rate = 1)

Gamma probability measure on `(0, Inf)` with unit rate.

# Arguments

- `shape`: positive gamma shape.
- `rate`: currently required to equal `1`.

# Fields

- `w`, `dom`, `symmetric`: measure data.
- `shapeParameter`, `rateParameter`: validated gamma parameters.
"""
struct GammaMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool
    shapeParameter::Float64
    rateParameter::Float64

    function GammaMeasure(shape::Real, rate::Real)
        shape <= 0 && throw(DomainError(shape, "shape parameter needs to be positive"))
        rate != 1 && throw(DomainError(rate, "rate must be unity (currently!)"))
        return new(build_w_gamma(shape), (0, Inf), false, shape, rate)
    end
end

"""
    LogisticMeasure()

Standard logistic probability measure on the real line.

# Fields

- `w`: logistic density.
- `dom`: `(-Inf, Inf)`.
- `symmetric`: `true`.
"""
struct LogisticMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Type_for_domain
    symmetric::Bool

    function LogisticMeasure()
        return new(w_logistic, (-Inf, Inf), true)
    end
end

export convert2affinePCE,
       calculateAffinePCE,
       assign2multi,
       evaluatePCE,
       sampleMeasure,
       samplePCE,
       mean,
       var,
       std

# auxiliary functions
function _createMethodVector(n::Int, word::String = "adaptiverejection")
    [word for i in 1:n]
end

function _createMethodVector(m::ProductMeasure, word::String = "adaptiverejection")
    _createMethodVector(length(m.measures), word)
end

function _createMethodVector(mop::MultiOrthoPoly, word::String = "adaptiverejection")
    _createMethodVector(mop.measure, word)
end

function _checkStandardDevation(σ::Real)
    σ < 0 && throw(DomainError(σ, "σ has to be non-negative"))
    σ < 1e-4 && @warn "σ is close to zero (σ = $σ)"
end

function _checkKind(kind::String)
    lowercase(kind) ∉ ["lbub", "μσ"] &&
        throw(DomainError(kind, "this kind is not supported"))
    lowercase(kind)
end

function _checkBounds(lb::Real, ub::Real)
    lb >= ub && throw(DomainError((lb, ub), "inconsistent bounds"))
end

function _checkNumberOfSamples(n::Int)
    n < 1 && throw(DomainError(n, "invalid number of samples"))
end
"""
    calculateAffinePCE(α::AbstractVector{<:Real})

Computes the affine PCE coefficients ``x_0`` and ``x_1`` from recurrence coefficients ``\alpha``.
"""
function calculateAffinePCE(α::AbstractVector{<:Real})
    length(α) < 1 && throw(DomainError(length(α), "not enough recursion coefficients"))
    return [α[1], 1.0]
end

calculateAffinePCE(op::AbstractOrthoPoly) = calculateAffinePCE(op.α)

function calculateAffinePCE(i::Int, mop::MultiOrthoPoly)
    p = length(mop.name)
    i > p &&
        throw(InconsistencyError("mop is $p-variate, but the PC for the $i-variate entry was requested"))
    calculateAffinePCE(mop.uni[i])
end

"""
```
convert2affinePCE(mu::Real, sigma::Real, op::AbstractCanonicalOrthoPoly; kind::String)
```

Computes the affine PCE coefficients ``x_0`` and ``x_1`` from

```math
X = a_1 + a_2 \\Xi = x_0 + x_1 \\phi_1(\\Xi),
```

where ``\\phi_1(t) = t-\\alpha_0`` is the first-order monic basis polynomial.

Works for subtypes of AbstractCanonicalOrthoPoly. The keyword `kind in ["lbub", "μσ"]`
specifies whether `p1` and `p2` have the meaning of lower/upper bounds or mean/standard deviation.
"""
function convert2affinePCE(a1::Real, a2::Real, α0::Real)
    [a1 + α0 * a2; a2]
end

function convert2affinePCE(mu::Real, sigma::Real, op::GaussOrthoPoly)
    _checkStandardDevation(sigma)
    convert2affinePCE(mu, sigma, first(op.α))
end

function convert2affinePCE(par1::Real, par2::Real, op::Uniform01OrthoPoly;
                           kind::String = "lbub")
    kind = _checkKind(kind)
    a1, a2 = if kind == "lbub"
        _checkBounds(par1, par2)
        par1, par2 - par1
    elseif kind == "μσ"
        _checkStandardDevation(par2)
        par1 - sqrt(3) * par2, 2 * sqrt(3) * par2
    end
    convert2affinePCE(a1, a2, first(op.α))
end

function convert2affinePCE(par1::Real, par2::Real, op::Uniform_11OrthoPoly;
                           kind::String = "lbub")
    kind = _checkKind(kind)
    a1, a2 = if kind == "lbub"
        _checkBounds(par1, par2)
        0.5 * (par1 + par2), 0.5 * (par2 - par1)
    elseif kind == "μσ"
        _checkStandardDevation(par2)
        par1, sqrt(3) * par2
    end
    convert2affinePCE(a1, a2, first(op.α))
end

function convert2affinePCE(p1::Real, p2::Real, op::Beta01OrthoPoly; kind::String = "lbub")
    kind = _checkKind(kind)
    α, β = op.measure.ashapeParameter, op.measure.bshapeParameter
    a1, a2 = if kind == "lbub"
        _checkBounds(p1, p2)
        a1, a2 = p1, p2 - p1
    elseif kind == "μσ"
        _checkStandardDevation(p2)
        a1, a2 = p1 - sqrt(α / β) * sqrt(1 + α + β) * p2,
                 (α + β) * sqrt((α + β + 1) / (α * β)) * p2
    end
    convert2affinePCE(a1, a2, first(op.α))
end

function convert2affinePCE(p1::Real, p2::Real, op::GammaOrthoPoly)
    throw(error("convert2affine not yet implemented for $(typeof(op))"))
end

function convert2affinePCE(p1::Real, p2::Real, op::LogisticOrthoPoly)
    _checkStandardDevation(p2)
    convert2affinePCE(p1, p2, first(op.α))
end

function assign2multi(x::AbstractVector{<:Real}, i::Int, ind::AbstractMatrix{<:Int})
    l, p = size(ind)
    nx, deg = length(x), ind[end, end]
    nx > deg + 1 &&
        throw(DomainError(nx, "inconsistent number of coefficients ($nx vs $(deg+1))"))
    i > p && throw(DomainError((i, p), "basis is $p-variate, you requested $i-variate"))
    myind = findUnivariateIndices(i, ind)[1:nx]
    y = spzeros(Float64, l)
    y[myind] = x
    return y
end

"""
__Univariate__

```
sampleMeasure(n::Int,name::String,w::Function,dom::Tuple{<:Real,<:Real},symm::Bool,d::Dict;method::String="adaptiverejection")
sampleMeasure(n::Int,m::Measure;method::String="adaptiverejection")
sampleMeasure(n::Int,op::AbstractOrthoPoly;method::String="adaptiverejection")
```

Draw `n` samples from the measure `m` described by its

  - `name`
  - weight function `w`,
  - domain `dom`,
  - symmetry property `symm`,
  - and, if applicable, parameters stored in the dictionary `d`.
    By default, an adaptive rejection sampling method is used (from [AdaptiveRejectionSampling.jl](https://github.com/mauriciogtec/AdaptiveRejectionSampling.jl)),
    unless it is a common random variable for which [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) is used.

The function is dispatched to accept `AbstractOrthoPoly`.

__Multivariate__

```
sampleMeasure(n::Int,m::ProductMeasure;method::Vector{String}=["adaptiverejection" for i=1:length(m.name)])
sampleMeasure(n::Int,mop::MultiOrthoPoly;method::Vector{String}=["adaptiverejection" for i=1:length(mop.meas.name)])
```

Multivariate extension, which provides an array of samples with `n` rows and
as many columns as the multimeasure has univariate measures.
"""
function sampleMeasure(n::Int, w::Function, dom::Tuple{<:Real, <:Real};
                       method::String = "adaptiverejection")
    _checkNumberOfSamples(n)
    method = lowercase(method)

    if method == "adaptiverejection"
        # works only if w is log-concave
        sampler = RejectionSampler(w, dom)
        return run_sampler!(sampler, n)
    elseif method == "rejection"
        # all purpose method but needs a solid envelope PDF
        throw(error("method $method not yet implemented"))
    elseif method == "inversecdf"
        throw(error("method $method not yet implemented"))
        # requires CDF and reliable root-finding method
    else
        throw(error("method $method not implemented"))
    end
end

function sampleMeasure(n::Int, meas::AbstractMeasure; method::String = "adaptiverejection")
    sampleMeasure(n, meas.w, meas.dom; method = method)
end

function sampleMeasure(n::Int, op::AbstractOrthoPoly; method::String = "adaptiverejection")
    sampleMeasure(n, op.measure; method = method)
end
function sampleMeasure(n::Int, op::AbstractCanonicalOrthoPoly)
    sampleMeasure(n, op.measure::AbstractCanonicalMeasure)
end

function sampleMeasure(n::Int, dist::Distribution{Univariate, Continuous})
    _checkNumberOfSamples(n)
    rand(dist, n)
end

sampleMeasure(n::Int, meas::GaussMeasure) = sampleMeasure(n, Normal())
sampleMeasure(n::Int, meas::Uniform01Measure) = sampleMeasure(n, Uniform())
function sampleMeasure(n::Int, meas::Beta01Measure)
    sampleMeasure(n, Beta(meas.ashapeParameter, meas.bshapeParameter))
end
function sampleMeasure(n::Int, meas::GammaMeasure)
    sampleMeasure(n, Gamma(meas.shapeParameter, 1 / meas.rateParameter))
end
sampleMeasure(n::Int, meas::LogisticMeasure) = sampleMeasure(n, Logistic())

function sampleMeasure(n::Int, meas::AbstractCanonicalMeasure;
                       method::String = "adaptiverejection")
    @warn "ignoring keyword method; sampling from Distributions.jl instead"
    sampleMeasure(n, meas)
end

function sampleMeasure(n::Int, measure::ProductMeasure;
                       method::Vector{String} = _createMethodVector(measure))
    samples = Matrix{Float64}(undef, n, 0)
    for (k, unimeasure) in enumerate(measure.measures)
        samples = hcat(samples, sampleMeasure(n, unimeasure; method = method[k]))
    end
    samples
end

function sampleMeasure(n::Int, mop::MultiOrthoPoly;
                       method::Vector{String} = _createMethodVector(mop))
    sampleMeasure(n, mop.measure; method = method)
end

"""
    evaluatePCE(x::AbstractVector{<:Real},ξ::AbstractVector{<:Real},α::AbstractVector{<:Real},β::AbstractVector{<:Real})

Evaluation of polynomial chaos expansion

```math
\\mathsf{x} = \\sum_{i=0}^{L} x_i \\phi_i{\\xi_j},
```

where `L+1 = length(x)` and ``x_j`` is the ``j``th sample where ``j=1,\\dots,m``
with `m = length(ξ)`.
"""
function evaluatePCE(x::AbstractVector{<:Real}, ξ::AbstractVector{<:Real},
                     α::AbstractVector{<:Real}, β::AbstractVector{<:Real})
    length(α) != length(β) &&
        throw(InconsistencyError("inconsistent number of recurrence coefficients"))
    Nsmpl = length(ξ)
    _checkNumberOfSamples(Nsmpl)
    Nx, Nrec = length(x), length(α)
    Nx > Nrec && throw(InconsistencyError("not enough recursion coefficients"))
    if Nrec > Nx
        α, β = α[1:Nx], β[1:Nx]
    end
    ϕ = zeros(Float64, Nsmpl, Nx)
    for n in 1:Nx
        ϕ[:, n] = evaluate(n - 1, ξ, α, β)
    end
    ϕ * x
end
function evaluatePCE(x::AbstractVector{<:Real}, ξ::Real, α::AbstractVector{<:Real},
                     β::AbstractVector{<:Real})
    evaluatePCE(x, [ξ], α, β)
end
function evaluatePCE(x::AbstractVector{<:Real}, ξ::AbstractVector{<:Real},
                     op::AbstractOrthoPoly)
    evaluatePCE(x, ξ, op.α, op.β)
end

function evaluatePCE(x::AbstractVector{<:Real}, ξ::AbstractMatrix{<:Real},
                     α::AbstractVector{<:AbstractVector{<:Real}},
                     β::AbstractVector{<:AbstractVector{<:Real}}, ind::AbstractMatrix{Int})
    Nsmpl = size(ξ, 1)
    _checkNumberOfSamples(Nsmpl)
    !(length(α) == length(β) == size(ξ, 2) == size(ind, 2)) &&
        throw(InconsistencyError("inconsistent number of coefficients"))
    Nx = length(x)
    Nx > size(ind, 1) &&
        throw(InconsistencyError("too few pc coefficients (resp: too small basis)"))
    ϕ = zeros(Float64, Nsmpl, Nx)
    for n in 1:Nx
        ϕ[:, n] = evaluate(ind[n, :], ξ, α, β)
    end
    ϕ * x
end
function evaluatePCE(x::AbstractVector{<:Real}, ξ::AbstractMatrix{<:Real},
                     mOP::MultiOrthoPoly)
    a, b = coeffs(mOP)
    evaluatePCE(x, ξ, a, b, mOP.ind)
end

"""
__Univariate__

```
samplePCE(n::Int,x::AbstractVector{<:Real},op::AbstractOrthoPoly;method::String="adaptiverejection")
```

Combines [`sampleMeasure`](@ref) and [`evaluatePCE`](@ref), i.e. it first draws `n` samples
from the measure, then evaluates the PCE for those samples.

__Multivariate__

```
samplePCE(n::Int,x::AbstractVector{<:Real},mop::MultiOrthoPoly;method::Vector{String}=["adaptiverejection" for i=1:length(mop.meas.name)])
```
"""
function samplePCE(n::Int, x::AbstractVector{<:Real}, op::AbstractOrthoPoly;
                   method::String = "adaptiverejection")
    ξ = sampleMeasure(n, op; method = method)
    evaluatePCE(x, ξ, op)
end

function samplePCE(n::Int, x::AbstractVector{<:Real}, op::AbstractCanonicalOrthoPoly;
                   method::String = "adaptiverejection")
    ξ = sampleMeasure(n, op)
    evaluatePCE(x, ξ, op)
end

function samplePCE(n::Int, x::AbstractVector{<:Real}, mop::MultiOrthoPoly;
                   method::Vector{String} = _createMethodVector(mop))
    ξ = sampleMeasure(n, mop; method = method)
    evaluatePCE(x, ξ, mop)
end

"""
__Univariate__

```
mean(x::AbstractVector,op::AbstractOrthoPoly)
```

__Multivariate__

```
mean(x::AbstractVector,mop::MultiOrthoPoly)
```

compute mean of random variable with PCE `x`
"""
mean(x::AbstractVector, op::AbstractOrthoPoly) = x[1] * computeSP2(0, op.β)

function mean(x::AbstractVector, mop::MultiOrthoPoly)
    nunc = length(mop.uni)
    x[1] * computeSP(zeros(Int64, nunc), mop)
end

"""
__Univariate__

```
var(x::AbstractVector,op::AbstractOrthoPoly)
var(x::AbstractVector,t2::Tensor)
```

__Multivariate__

```
var(x::AbstractVector,mop::MultiOrthoPoly)
var(x::AbstractVector,t2::Tensor)
```

compute variance of random variable with PCE `x`
"""
function var(x::AbstractVector, op::AbstractOrthoPoly)
    t = computeSP2(op)
    # length(t2) > length(x) && throw(InconsistencyError("cannot compute variance; too many PCE coefficients"))
    sum(x[i]^2 * t[i] for i in 2:length(x))
end

function var(x::AbstractVector, t2::AbstractTensor)
    sum(x[i]^2 * t2.get([i - 1, i - 1]) for i in 2:length(x))
end

function var(x::AbstractVector, mop::MultiOrthoPoly)
    length(x) > size(mop.ind, 1) &&
        throw(InconsistencyError("cannot compute variance; too many PCE coefficients"))
    var(x, Tensor(2, mop))
end

"""
__Univariate__

```
std(x::AbstractVector,op::AbstractOrthoPoly)
```

__Multivariate__

```
std(x::AbstractVector,mop::MultiOrthoPoly)
```

compute standard deviation of random variable with PCE `x`
"""
std(x::AbstractVector, op::AbstractOrthoPoly) = sqrt(var(x, op))
std(x::AbstractVector, mop::MultiOrthoPoly) = sqrt(var(x, mop))

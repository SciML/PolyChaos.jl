export OrthoPoly,
       MultiOrthoPoly,
       LegendreOrthoPoly,
       JacobiOrthoPoly,
       LaguerreOrthoPoly,
       genLaguerreOrthoPoly,
       HermiteOrthoPoly,
       genHermiteOrthoPoly,
       GaussOrthoPoly,
       Beta01OrthoPoly,
       GammaOrthoPoly,
       LogisticOrthoPoly,
       Uniform01OrthoPoly,
       Uniform_11OrthoPoly,
       MeixnerPollaczekOrthoPoly,
       InconsistencyError,
       Quad

struct InconsistencyError <: Exception
    var::String
end

Base.showerror(io::IO, err::InconsistencyError) = print(io, err.var)

function _checkConsistency(deg::Int, Nrec::Int)
    deg < 0 && throw(DomainError(deg, "degree has to be non-negative"))
    Nrec < deg + 1 && throw(DomainError(Nrec,
        "not enough recurrence coefficients specified (need >= $(deg + 1))"))
end

struct OrthoPoly{V <: AbstractVector{<:Real}, M, Q} <: AbstractOrthoPoly{M, Q}
    name::String
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function OrthoPoly(name::String, deg::Int, α::AbstractVector{<:Real},
        β::AbstractVector{<:Real}, measure::AbstractMeasure;
        addQuadrature::Bool = true)
    deg < 0 && throw(DomainError(deg, "degree has to be non-negative"))
    !(length(α) == length(β)) && throw(InconsistencyError("Inconsistent lengths"))
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    # @show types = promote_type(eltype(α), eltype(β)), promote_type(typeof(α), typeof(β)), typeof(measure), typeof(quadrature)
    OrthoPoly{promote_type(typeof(α), typeof(β)), typeof(measure), typeof(quadrature)}(
        lowercase(name),
        deg,
        α,
        β,
        measure,
        quadrature)
end

# constructor for known Measure
function OrthoPoly(name::String, deg::Int, measure::AbstractMeasure; Nrec = deg + 1,
        Nquad = 10 * Nrec, quadrature::Function = clenshaw_curtis,
        discretization::Function = stieltjes, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    name = lowercase(name)
    α,
    β = rm_compute(measure, Nrec, Nquad, quadrature = quadrature,
        discretization = discretization)
    OrthoPoly(name, deg, α, β, measure, addQuadrature = addQuadrature)
end

function OrthoPoly(μ::Measure, deg::Int; Nrec = deg + 1, Nquad = 10 * Nrec,
        quadrature::Function = clenshaw_curtis,
        discretization::Function = stieltjes, addQuadrature::Bool = true)
    OrthoPoly(μ.name, deg, μ; Nrec = Nrec, Nquad = Nquad, quadrature = quadrature,
        discretization = discretization, addQuadrature = addQuadrature)
end

# general constructor
function OrthoPoly(name::String, deg::Int, w::Function, s::Tuple{<:Real, <:Real},
        symm::Bool, d::Dict = Dict(); Nrec = deg + 1, Nquad = 10 * Nrec,
        quadrature::Function = clenshaw_curtis,
        discretization::Function = stieltjes, addQuadrature::Bool = true)
    name = lowercase(name)
    measure = Measure(name, w, s, symm, d)
    OrthoPoly(name, deg, measure; Nrec = Nrec, Nquad = Nquad, quadrature = quadrature,
        discretization = discretization, addQuadrature = addQuadrature)
end

##############################################

struct LegendreOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int # maximum degree
    α::V # recurrence coefficients
    β::V # recurrence coefficients
    measure::M
    quad::Q
end

function LegendreOrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_legendre(Nrec)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    LegendreOrthoPoly{promote_type(typeof(α), typeof(β)), LegendreMeasure,
        typeof(quadrature)}(deg, α,
        β,
        LegendreMeasure(),
        quadrature)
end

function OrthoPoly(::LegendreMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    LegendreOrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct JacobiOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int            # maximum degree
    α::V # recurrence coefficients
    β::V # recurrence coefficients
    measure::M
    quad::Q
end

function JacobiOrthoPoly(deg::Int, shape_a::Real, shape_b::Real; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_jacobi(Nrec, shape_a, shape_b)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    JacobiOrthoPoly{promote_type(typeof(α), typeof(β)), JacobiMeasure, typeof(quadrature)}(
        deg,
        α,
        β,
        JacobiMeasure(shape_a,
            shape_b),
        quadrature)
end

function OrthoPoly(μ::JacobiMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    JacobiOrthoPoly(deg, μ.ashapeParameter, μ.bshapeParameter; Nrec = Nrec,
        addQuadrature = addQuadrature)
end

struct LaguerreOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V # recurrence coefficients
    β::V # recurrence coefficients
    measure::M
    quad::Q
end

function LaguerreOrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_laguerre(Nrec)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    LaguerreOrthoPoly{promote_type(typeof(α), typeof(β)), LaguerreMeasure,
        typeof(quadrature)}(deg, α,
        β,
        LaguerreMeasure(),
        quadrature)
end

function OrthoPoly(::LaguerreMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    LaguerreOrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct genLaguerreOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function genLaguerreOrthoPoly(deg::Int, shape::Real; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_laguerre(Nrec, shape)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    genLaguerreOrthoPoly{promote_type(typeof(α), typeof(β)), genLaguerreMeasure,
        typeof(quadrature)}(deg,
        α,
        β,
        genLaguerreMeasure(shape),
        quadrature)
end

function OrthoPoly(μ::genLaguerreMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    genLaguerreOrthoPoly(deg, μ.shapeParameter; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct HermiteOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function HermiteOrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_hermite(Nrec)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    HermiteOrthoPoly{promote_type(typeof(α), typeof(β)), HermiteMeasure, typeof(quadrature)
    }(deg, α,
        β,
        HermiteMeasure(),
        quadrature)
end

function OrthoPoly(::HermiteMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    HermiteOrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct genHermiteOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function genHermiteOrthoPoly(deg::Int, mu::Real; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_hermite(Nrec, mu)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    genHermiteOrthoPoly{promote_type(typeof(α), typeof(β)), genHermiteMeasure,
        typeof(quadrature)}(deg,
        α, β,
        genHermiteMeasure(mu),
        quadrature)
end

function OrthoPoly(μ::genHermiteMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    genHermiteOrthoPoly(deg, μ.muParameter; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct MeixnerPollaczekOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V              # recurrence coefficients
    β::V              # recurrence coefficients
    measure::M
    quad::Q
end

function MeixnerPollaczekOrthoPoly(deg::Int, λ::Real, ϕ::Real; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = rm_meixner_pollaczek(Nrec, λ, ϕ)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    MeixnerPollaczekOrthoPoly{promote_type(typeof(α), typeof(β)), MeixnerPollaczekMeasure,
        typeof(quadrature)}(deg,
        α,
        β,
        MeixnerPollaczekMeasure(λ, ϕ),
        quadrature)
end

function OrthoPoly(μ::MeixnerPollaczekMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    MeixnerPollaczekOrthoPoly(deg, μ.λParameter, μ.ϕParameter; Nrec = Nrec,
        addQuadrature = addQuadrature)
end

struct GaussOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function GaussOrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = r_scale(1 / sqrt(2pi), rm_hermite_prob(Nrec)...)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    GaussOrthoPoly{promote_type(typeof(α), typeof(β)), GaussMeasure, typeof(quadrature)}(
        deg,
        α,
        β,
        GaussMeasure(),
        quadrature)
end

function OrthoPoly(::GaussMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    GaussOrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct Uniform01OrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function Uniform01OrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = r_scale(1.0, rm_legendre01(Nrec)...)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    Uniform01OrthoPoly{promote_type(typeof(α), typeof(β)), Uniform01Measure,
        typeof(quadrature)}(deg,
        α, β,
        Uniform01Measure(),
        quadrature)
end

function OrthoPoly(::Uniform01Measure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    Uniform01OrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct Uniform_11OrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V  # recurrence coefficients
    β::V  # recurrence coefficients
    measure::M
    quad::Q
end

function Uniform_11OrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = r_scale(0.5, rm_legendre(Nrec)...)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    Uniform_11OrthoPoly{promote_type(typeof(α), typeof(β)), Uniform_11Measure,
        typeof(quadrature)}(deg,
        α, β,
        Uniform_11Measure(),
        quadrature)
end

function OrthoPoly(::Uniform_11Measure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    Uniform_11OrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

struct Beta01OrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V # recurrence coefficients
    β::V # recurrence coefficients
    measure::M
    quad::Q
end

function Beta01OrthoPoly(deg::Int, shape_a::Real, shape_b::Real; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α,
    β = r_scale(1 / beta(shape_a, shape_b),
        rm_jacobi01(Nrec, shape_b - 1.0, shape_a - 1.0)...)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    Beta01OrthoPoly{promote_type(typeof(α), typeof(β)), Beta01Measure, typeof(quadrature)}(
        deg,
        α,
        β,
        Beta01Measure(shape_a,
            shape_b),
        quadrature)
end

function OrthoPoly(μ::Beta01Measure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    Beta01OrthoPoly(deg, μ.ashapeParameter, μ.bshapeParameter; Nrec = Nrec,
        addQuadrature = addQuadrature)
end

struct GammaOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V # recurrence coefficients
    β::V # recurrence coefficients
    measure::M
    quad::Q
end

function GammaOrthoPoly(deg::Int, shape::Real, rate::Real; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = r_scale((rate^shape) / gamma(shape), rm_laguerre(Nrec, shape - 1.0)...)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    GammaOrthoPoly{promote_type(typeof(α), typeof(β)), GammaMeasure, typeof(quadrature)}(
        deg,
        α,
        β,
        GammaMeasure(shape,
            rate),
        quadrature)
end

function OrthoPoly(μ::GammaMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    GammaOrthoPoly(deg, μ.shapeParameter, μ.rateParameter; Nrec = Nrec,
        addQuadrature = addQuadrature)
end

struct LogisticOrthoPoly{V, M, Q} <: AbstractCanonicalOrthoPoly{V, M, Q}
    deg::Int          # maximum degree
    α::V # recurrence coefficients
    β::V # recurrence coefficients
    measure::M
    quad::Q
end

function LogisticOrthoPoly(deg::Int; Nrec::Int = deg + 1, addQuadrature::Bool = true)
    _checkConsistency(deg, Nrec)
    α, β = r_scale(1.0, rm_logistic(Nrec)...)
    quadrature = addQuadrature ? Quad(length(α) - 1, α, β) : EmptyQuad()
    LogisticOrthoPoly{promote_type(typeof(α), typeof(β)), LogisticMeasure,
        typeof(quadrature)}(deg, α,
        β,
        LogisticMeasure(),
        quadrature)
end

function OrthoPoly(::LogisticMeasure, deg::Int; Nrec::Int = deg + 1,
        addQuadrature::Bool = true)
    LogisticOrthoPoly(deg; Nrec = Nrec, addQuadrature = addQuadrature)
end

# #####################################################
# #####################################################
# #####################################################

struct MultiOrthoPoly{M, Q, V <: AbstractVector} <: AbstractOrthoPoly{M, Q}
    name::Vector{String}
    deg::Int
    dim::Int
    ind::Matrix{Int} # multi-index
    measure::ProductMeasure
    uni::V
    function MultiOrthoPoly(uniOrthoPolys, deg::Int)
        degs = [op.deg for op in uniOrthoPolys]
        deg > minimum(degs) && throw(DomainError(deg,
            "Requested degree $deg is greater than smallest univariate degree $(minimum(degs))."))

        w(t) = prod([op.measure.w(t) for op in uniOrthoPolys])
        measures = [op.measure for op in uniOrthoPolys]
        measure = ProductMeasure(w, measures)

        names = [hasfield(typeof(op), :name) ? op.name : string(typeof(op))
                 for op in uniOrthoPolys]
        ind = calculateMultiIndices(length(uniOrthoPolys), deg)
        dim = size(ind, 1)

        new{typeof(measure), typeof(first(uniOrthoPolys).quad), typeof(uniOrthoPolys)}(
            names,
            deg,
            dim,
            ind,
            measure,
            uniOrthoPolys)
    end
end

function _hasfield(op::AbstractOrthoPoly, name::Symbol)
    name in fieldnames(op)
end

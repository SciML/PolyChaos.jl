export          OrthoPoly,
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
                MeixnerPollaczekOrthoPoly,
                InconsistencyError,
                Quad,
                Tensor

struct InconsistencyError <: Exception
    var::String
end

Base.showerror(io::IO, err::InconsistencyError) = print(io, err.var)


function _checkConsistency(deg::Int, Nrec::Int)
    deg < 0 && throw(DomainError(deg, "degree has to be non-negative"))
    Nrec < deg + 1 && throw(DomainError(Nrec, "not enough recurrence coefficients specified (need >= $(deg + 1))"))
end

struct LegendreOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int            # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::LegendreMeasure
    quad::AbstractQuad

    # inner constructor
    function LegendreOrthoPoly(deg::Int; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_legendre(Nrec)
        new(deg, α, β, LegendreMeasure(), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct JacobiOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int            # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::JacobiMeasure
    quad::AbstractQuad

    # inner constructor
    function JacobiOrthoPoly(deg::Int, shape_a::Real, shape_b::Real; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_jacobi(Nrec, shape_a, shape_b)
        new(deg, α, β, JacobiMeasure(shape_a, shape_b), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct LaguerreOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::LaguerreMeasure
    quad::AbstractQuad

    # inner constructor
    function LaguerreOrthoPoly(deg::Int;Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_laguerre(Nrec)
        new(deg, α, β, LaguerreMeasure(), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct genLaguerreOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::genLaguerreMeasure
    quad::AbstractQuad

    # inner constructor
    function genLaguerreOrthoPoly(deg::Int, shape::Real ;Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_laguerre(Nrec, shape)
        new(deg, α, β, genLaguerreMeasure(shape), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct HermiteOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::HermiteMeasure
    quad::AbstractQuad

    # inner constructor
    function HermiteOrthoPoly(deg::Int; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_hermite(Nrec)
        new(deg, α, β, HermiteMeasure(), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct genHermiteOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::genHermiteMeasure
    quad::AbstractQuad

    # inner constructor
    function genHermiteOrthoPoly(deg::Int, mu::Real;Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_hermite(Nrec, mu)
        new(deg, α, β, genHermiteMeasure(mu), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct MeixnerPollaczekOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::MeixnerPollaczekMeasure
    quad::AbstractQuad

    # inner constructor
    function MeixnerPollaczekOrthoPoly(deg::Int, λ::Real, ϕ::Real; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = rm_meixner_pollaczek(Nrec, λ, ϕ)
        new(deg, α, β, MeixnerPollaczekMeasure(λ,ϕ), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct GaussOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::GaussMeasure
    quad::AbstractQuad

    # inner constructor
    function GaussOrthoPoly(deg::Int;Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = r_scale(1/sqrt(2pi), rm_hermite_prob(Nrec)...)
        new(deg, α, β, GaussMeasure(), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct Uniform01OrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::Uniform01Measure
    quad::AbstractQuad

    # inner constructor
    function Uniform01OrthoPoly(deg::Int;Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = r_scale(1., rm_legendre01(Nrec)...)
        new(deg, α, β, Uniform01Measure(), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct Beta01OrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::Beta01Measure
    quad::AbstractQuad

    # inner constructor
    function Beta01OrthoPoly(deg::Int, shape_a::Real, shape_b::Real; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = r_scale(1/beta(shape_a, shape_b), rm_jacobi01(Nrec, shape_b-1., shape_a-1.)...)
        new(deg, α, β, Beta01Measure(shape_a, shape_b), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct GammaOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::GammaMeasure
    quad::AbstractQuad

    # inner constructor
    function GammaOrthoPoly(deg::Int, shape::Real, rate::Real; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = r_scale((rate^shape)/gamma(shape), rm_laguerre(Nrec, shape-1.)...)
        new(deg, α, β, GammaMeasure(shape, rate), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct LogisticOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::LogisticMeasure
    quad::AbstractQuad

    # inner constructor
    function LogisticOrthoPoly(deg::Int; Nrec::Int=deg+1, addQuadrature::Bool = true)
        _checkConsistency(deg, Nrec)
        α, β = r_scale(1., rm_logistic(Nrec)...)
        new(deg, α, β, LogisticMeasure(), addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

struct OrthoPoly <: AbstractOrthoPoly
    name::String
    deg::Int          # maximum degree
    α::Vector{<:Real}  # recurrence coefficients
    β::Vector{<:Real}  # recurrence coefficients
    measure::AbstractMeasure
    quad::AbstractQuad
    # inner constructor
    function OrthoPoly(name::String, deg::Int, α::Vector{<:Real}, β::Vector{<:Real}, measure::AbstractMeasure; addQuadrature::Bool = true)
        deg < 0 && throw(DomainError(deg, "degree has to be non-negative"))
        !(length(α) == length(β)) && throw(InconsistencyError("Inconsistent lengths"))
        new(lowercase(name), deg, α, β, measure, addQuadrature ?  Quad(length(α)-1,α,β) : EmptyQuad())
    end
end

# constructor for known Measure
function OrthoPoly(name::String, deg::Int, measure::AbstractMeasure; Nrec=deg+1, Nquad=10*Nrec, quadrature::Function=clenshaw_curtis, discretization::Function=stieltjes)
    _checkConsistency(deg, Nrec)
    name = lowercase(name)
    α, β = rm_compute(measure, Nrec, Nquad, quadrature=quadrature, discretization=discretization)
    return OrthoPoly(name, deg, α, β, measure)
end

# general constructor
function OrthoPoly(name::String, deg::Int, w::Function, s::Tuple{<:Real,<:Real}, symm::Bool, d::Dict=Dict(); Nrec=deg+1, Nquad=10*Nrec, quadrature::Function=clenshaw_curtis, discretization::Function=stieltjes)
  name = lowercase(name)
  measure = Measure(name, w, s, symm, d)
  OrthoPoly(name, deg, measure; Nrec=Nrec, Nquad=Nquad, quadrature=quadrature, discretization=discretization)
end

#####################################################
#####################################################
#####################################################

struct MultiOrthoPoly <: AbstractOrthoPoly
name::Vector{String}
deg::Int
dim::Int
ind::Matrix{<:Int} # multi-index
measure::ProductMeasure
uni::Vector{<:AbstractOrthoPoly}
    function MultiOrthoPoly(uniOrthoPolys::Vector{<:AbstractOrthoPoly}, deg::Int)
      degs = [ op.deg for op in uniOrthoPolys ]
      deg > minimum(degs) && throw(DomainError(deg, "Requested degree $deg is greater than smallest univariate degree $(minimum(degs))."))

      w(t) = prod([op.measure.w(t) for op in uniOrthoPolys])
      measures = [ op.measure for op in uniOrthoPolys ]
      measure = ProductMeasure(w, measures)

      names = [ supertype(typeof(op)) == AbstractCanonicalOrthoPoly ? string(typeof(op)) : op.name for op in uniOrthoPolys ]
      ind = calculateMultiIndices(length(uniOrthoPolys), deg)
      dim = size(ind,1)

      new(names, deg, dim, ind, measure, uniOrthoPolys)
    end
end

struct Tensor
dim::Int          # "dimension"
T::SparseVector{Float64,Int}
get::Function
op::AbstractOrthoPoly

    # inner constructors
    function Tensor(dim::Int,mop::MultiOrthoPoly)
        tensorEntries = computeTensorizedSP(dim, mop)
        getfun(ind) = getentry(ind, tensorEntries, mop.ind, dim)
        new(dim, tensorEntries, getfun, mop)
    end
    function Tensor(dim::Int,opq::AbstractOrthoPoly)
        tensorEntries = computeTensorizedSP(dim, opq)
        getfun(ind) = getentry(ind, tensorEntries, calculateMultiIndices(1, opq.deg), dim)
        new(dim, tensorEntries, getfun, opq)
    end
end

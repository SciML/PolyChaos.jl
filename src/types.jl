export  AbstractMeasure,
        AbstractCanonicalMeasure,
        AbstractOrthoPoly,
        AbstractCanonicalOrthoPoly,
        AbstractQuad,
        OrthoPoly,
        MultiOrthoPoly,
        Quad,
        Measure,
        MultiMeasure,
        Tensor,
        LegendreMeasure,
        JacobiMeasure,
        LaguerreMeasure,
        genLaguerreMeasure,
        HermiteMeasure,
        genHermiteMeasure,
        GaussianMeasure,
        Beta01Measure,
        GammaMeasure,
        Uniform01Measure,
        LogisticMeasure,
        MeixnerPollaczekMeasure,
        LegendreOrthoPoly,
        JacobiOrthoPoly,
        LaguerreOrthoPoly,
        genLaguerreOrthoPoly,
        HermiteOrthoPoly,
        genHermiteOrthoPoly,
        GaussianOrthoPoly,
        Beta01OrthoPoly,
        GammaOrthoPoly,
        LogisticOrthoPoly,
        Uniform01OrthoPoly

abstract type AbstractMeasure end
abstract type AbstractCanonicalMeasure <: AbstractMeasure end
abstract type AbstractQuad end
abstract type AbstractOrthoPoly end
abstract type AbstractCanonicalOrthoPoly <: AbstractOrthoPoly end

struct Measure <: AbstractMeasure
    name::String
    w::Function
    dom::Tuple{Float64,Float64}
    symmetric::Bool
    pars::Dict
    function Measure(name::String,w::Function,dom::Tuple{Real,Real},symm::Bool,d::Dict=Dict())
        @assert dom[1] < dom[2] "Invalid support."
        new(lowercase(name),w,Float64.(dom),symm,d)
    end
end

struct MultiMeasure
    name::Vector{String}
    w::Function
    w_uni::Vector{Function}
    dom::Vector{Tuple{Float64,Float64}}
    symmetric::Vector{Bool}
    pars::Vector{Dict}
end

# constructor for classic distributions
struct LegendreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    function LegendreMeasure()
        new(w_legendre,(-1.,1.),true)
    end
end

struct JacobiMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    ashapeParameter::Real
    bshapeParameter::Real
    
    function JacobiMeasure(shape_a::Real, shape_b::Real)
        any
        shape_a <= -1 && throw(DomainError(shape_a, "shape parameter a must be > -1"))
        shape_b <= -1 && throw(DomainError(shape_b, "shape parameter b must be > -1"))
        new(build_w_jacobi(shape_a,shape_b), (-1,1), isapprox(shape_a, shape_b), shape_a, shape_b)
    end
end

struct LaguerreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function LaguerreMeasure()
        new(w_laguerre, (0,Inf), false)
    end
end

struct genLaguerreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    shapeParameter::Real

    function genLaguerreMeasure(shape::Real)
        shape <= -1 && throw(DomainError(shape, "invalid shape parameter"))
        new(build_w_genlaguerre(shape), (0,Inf), false, shape)
    end
end

struct HermiteMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function HermiteMeasure()
        new(w_hermite, (-Inf,Inf), true)
    end
end

struct genHermiteMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    muParameter::Real

    function genHermiteMeasure(mu::Real)
        mu <= -0.5 && throw(DomainError(mu, "invalid parameter value (must be > - 0.5)"))
        new(build_w_genhermite(mu), (-Inf,Inf), true, mu)
    end
end

struct MeixnerPollaczekMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    λParameter::Real
    ϕParameter::Real
    
    function MeixnerPollaczekMeasure(λ::Real, ϕ::Real)
        λ <= 0 && throw(DomainError(λ, "λ has to be positive"))
        !(0 < ϕ < pi) && throw(DomainError(ϕ, "ϕ has to be between 0 and pi"))
        new(build_w_meixner_pollaczek(λ, ϕ), (-Inf,Inf), false, λ, ϕ)
    end
end

struct GaussianMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function GaussianMeasure()
        new(w_gaussian, (-Inf,Inf), true)
    end
end

struct Uniform01Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function Uniform01Measure()
        new(w_uniform01, (0,1), true)
    end
end

struct Beta01Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    ashapeParameter::Real
    bshapeParameter::Real
    
    function Beta01Measure(a::Real, b::Real)
        a <= 0 && throw(DomainError(a, "shape parameter a must be positive"))
        b <= 0 && throw(DomainError(b, "shape parameter b must be positive"))
        new(build_w_beta(a,b), (0,1), isapprox(a,b), a, b)
    end
end

struct GammaMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    shapeParameter::Real
    rateParameter::Real

    function GammaMeasure(shape::Real, rate::Real)
        shape <= 0 && throw(DomainError(shape, "shape parameter needs to be positive"))
        rate != 1 && throw(DomainError(rate, "rate must be unity (currently!)"))
        new(build_w_gamma(shape), (0,Inf), false, shape, rate)
    end
end

struct LogisticMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function LogisticMeasure()
        new(w_logistic, (-Inf,Inf), true)
    end
end


function _checkConsistency(deg::Int, Nrec::Int)
    deg < 0 && throw(DomainError(deg, "degree has to be non-negative"))
    Nrec < deg + 1 && throw(DomainError(Nrec, "not enough recurrence coefficients specified (need >= $(deg + 1))"))
end

struct OrthoPoly <: AbstractOrthoPoly
    name::String
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad
    # inner constructor
    function OrthoPoly(name::String,deg::Int64,α::Vector{Float64},β::Vector{Float64},m::Measure)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert length(α) == length(β) "Different number of recursion coefficients α and β supplied."
        new(lowercase(name),deg,α,β,m,EmptyQuad())
    end
end

struct LegendreOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function LegendreOrthoPoly(deg::Int;Nrec::Int=deg+1)
        _checkConsistency(deg, Nrec)
        a, b = rm_legendre(Nrec)
        new(deg, a, b, EmptyQuad())
    end
end

struct JacobiOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    shapeParameters::Dict{Symbol, Real}
    quad::AbstractQuad

    # inner constructor
    function JacobiOrthoPoly(deg::Int, shapeParameters::Dict{Symbol, Real}; Nrec::Int=deg+1)
        _checkConsistency(deg, Nrec)
        a, b = rm_jacobi(Nrec, shapeParameters[:shape_a], shapeParameters[:shape_b])
        new(deg, a, b, shapeParameters, EmptyQuad())
    end
end

struct LaguerreOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function LaguerreOrthoPoly(deg::Int64,;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_laguerre(Nrec),EmptyQuad())
    end
end

struct genLaguerreOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function genLaguerreOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_laguerre(Nrec,Float64(d[:shape])),EmptyQuad())
    end
end

struct HermiteOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function HermiteOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_hermite(Nrec),EmptyQuad())
    end
end

struct genHermiteOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function genHermiteOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_hermite(Nrec,Float64(d[:mu])),EmptyQuad())
    end
end

struct MeixnerpollaczekOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function MeixnerpollaczekOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_meixner_pollaczek(Nrec,Float64(d[:lambda]),Float64(d[:phi])),EmptyQuad())
    end
end

struct GaussianOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function GaussianOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1/sqrt(2pi),rm_hermite_prob(Nrec)),EmptyQuad())
    end
end

struct Uniform01OrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function Uniform01OrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1.,rm_legendre01(Nrec)),EmptyQuad())
    end
end

struct Beta01OrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function Beta01OrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1/beta(d[:shape_a],d[:shape_b]),rm_jacobi01(Nrec,Float64(d[:shape_b])-1,Float64(d[:shape_a])-1)),EmptyQuad())
    end
end

struct GammaOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function GammaOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale((d[:rate]^d[:shape])/gamma(d[:shape]),rm_laguerre(Nrec,Float64(d[:shape]-1.))),EmptyQuad())
    end
end

struct LogisticOrthoPoly <: AbstractCanonicalOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function LogisticOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1.,rm_logistic(Nrec)),EmptyQuad())
    end
end

"""
# constructor for classic distributions
function OrthoPoly(name::String,deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
  @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
  name = lowercase(name)
  name == "legendre"         && return OrthoPoly(name,deg,rm_legendre(Nrec)...,Measure(name,Dict()))
  name == "jacobi"           && return OrthoPoly(name,deg,rm_jacobi(Nrec,Float64(d[:shape_a]), Float64(d[:shape_b]))...,Measure(name,d))
  name == "laguerre"         && return OrthoPoly(name,deg,rm_laguerre(Nrec)...,Measure(name,Dict()))
  name == "genlaguerre"      && return OrthoPoly(name,deg,rm_laguerre(Nrec,Float64(d[:shape]))...,Measure(name,d))
  name == "hermite"          && return OrthoPoly(name,deg,rm_hermite(Nrec)...,Measure(name,Dict()))
  name == "genhermite"       && return OrthoPoly(name,deg,rm_hermite(Nrec,Float64(d[:mu]))...,Measure(name,d))
  name == "meixnerpollaczek" && return OrthoPoly(name,deg,rm_meixner_pollaczek(Nrec,Float64(d[:lambda]),Float64(d[:phi]))...,Measure(name,d))
  # measures corresponding to probability density functions:
  name == "gaussian"  && return OrthoPoly(name,deg,r_scale(1/sqrt(2pi),rm_hermite_prob(Nrec)...)...,Measure(name,Dict()))
  name == "uniform01" && return OrthoPoly(name,deg,r_scale(1.,rm_legendre01(Nrec)...)...,Measure(name,Dict()))
  name == "beta01"    && return OrthoPoly(name,deg,r_scale(1/beta(d[:shape_a],d[:shape_b]),rm_jacobi01(Nrec,Float64(d[:shape_b])-1,Float64(d[:shape_a])-1)...)...,Measure(name,d))
  name == "gamma"     && return OrthoPoly(name,deg,r_scale((d[:rate]^d[:shape])/gamma(d[:shape]),rm_laguerre(Nrec,Float64(d[:shape]-1.))...)...,Measure(name,d))
  name == "logistic"  && return OrthoPoly(name,deg,r_scale(1.,rm_logistic(Nrec)...)...,Measure(name,Dict()))
end
"""
# constructor for known Measure
function OrthoPoly(name::String,deg::Int64,m::Measure;Nrec=deg+1,Nquad=10*Nrec,quadrature::Function=clenshaw_curtis,discretization::Function=stieltjes)
  @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
  name = lowercase(name)
  a,b = rm_compute(m,Nrec,Nquad,quadrature=quadrature,discretization=discretization)
  return OrthoPoly(name,deg,a,b,m)
end

# general constructor
function OrthoPoly(name::String,deg::Int64,w::Function,s::Tuple{Real,Real},symm::Bool,d::Dict=Dict();Nrec=deg+1,Nquad=10*Nrec,quadrature::Function=clenshaw_curtis,discretization::Function=stieltjes)
  name = lowercase(name)
  m = Measure(name,w,s,symm,d)
  OrthoPoly(name,deg,m;Nrec=Nrec,Nquad=Nquad,quadrature=quadrature,discretization=discretization)
end

struct Quad
    name::String              # name of quadrature
    Nquad::Int64              # number of qudrature points
    nodes::Vector{Float64}
    weights::Vector{Float64}
    meas::Measure
    function Quad(name::String,N::Int64,nodes::Vector{Float64},weights::Vector{Float64},m::Measure)
        @assert N >= 1 "Number of qudrature points has to be positive"
        @assert length(nodes) == length(weights) "Inconsistent number of nodes and weights inconsistent."
        new(lowercase(name),N,nodes,weights,m)
    end
end

struct EmptyQuad <: AbstractQuad
    function EmptyQuad()
        new()
    end
end

# general constructor
function Quad(N::Int,α::Vector{Float64},β::Vector{Float64},m::Measure)
    @assert length(α) == length(β) "Inconsistent length of recurrence coefficients."
    @assert N <= length(α) - 1 "Requested number of quadrature points $N cannot be provided with $(length(α)) recurrence coefficients"
    n,w = gauss(N,α,β)
    Quad("golubwelsch",N,n,w,m)
end
Quad(N::Int,op::OrthoPoly) = Quad(N,op.α,op.β,op.meas)

function Quad(N::Int64,weight::Function,α::Vector{Float64},β::Vector{Float64},supp::Tuple{Float64,Float64},symm::Bool,d::Dict=Dict())
    m = Measure("fun_"*String(nameof(weight)),weight,supp,symm,d)
    Quad(N,α,β,m)
end

# all-purpose constructor (last resort!)
function Quad(N::Int64,m::Measure;quadrature::Function=clenshaw_curtis)
  n, w = quadgp(m.w,m.dom[1],m.dom[2],N;quadrature=quadrature)
  Quad("quadgp",N,n,w,m)
end

function Quad(N::Int64,weight::Function,supp::Tuple{Real,Real},symm::Bool,d::Dict=Dict();quadrature::Function=clenshaw_curtis)
    @assert N >= 1 "Number of qudrature points has to be positive"
    m = Measure("fun_"*String(nameof(weight)),weight,supp,symm,d)
    Quad(N,m;quadrature=quadrature)
end

# Struct that contains pre-computed nodes and weights
struct OrthoPolyQ
    op::OrthoPoly
    quad::Quad
end

function OrthoPolyQ(op::OrthoPoly,N::Int64)
    q = Quad(N,op.α,op.β,op.meas)
    return OrthoPolyQ(op,q)
end
OrthoPolyQ(op::OrthoPoly) = OrthoPolyQ(op,length(op.α)-1)

function OrthoPolyQ(name::String,N::Int64,d::Dict=Dict();Nrec::Int64=N+1)
    op = OrthoPoly(name,N,d;Nrec=Nrec)
    OrthoPolyQ(op)
end


struct MultiOrthoPoly
name::Vector{String}
deg::Int64
dim::Int64
ind::Matrix{Int64} # multi-index
meas::MultiMeasure
uni::Union{Vector{OrthoPoly},Vector{OrthoPolyQ}}
    function MultiOrthoPoly(uni::Union{Vector{OrthoPoly},Vector{OrthoPolyQ}},deg::Int64)
      t = typeof(uni) == Vector{OrthoPoly}
      degs = [ t ? u.deg : u.op.deg for u in uni ]
      @assert deg <= minimum(degs) "Requested degree $deg is greater than smallest univariate degree $(minimum(degs))."
      Nuni::Int64 = length(uni)

      name_meas = [ t ? u.meas.name : u.op.meas.name for u in uni ]
      supp      = [ t ? u.meas.dom : u.op.meas.dom for u in uni       ]
      symm      = [ t ? u.meas.symmetric : u.op.meas.symmetric for u in uni ]
      pars      = [ t ? u.meas.pars : u.op.meas.pars for u in uni      ]
      w_uni     = [ t ? u.meas.w : u.op.meas.w for u in uni         ]
      w(t) = t ? prod([uni[i].meas.w(t[i]) for i=1:Nuni]) : prod([uni[i].op.meas.w(t[i]) for i=1:Nuni])
      m = MultiMeasure(name_meas,w,w_uni,supp,symm,pars)

      name = [ t ? u.name : u.op.name for u in uni ]
      ind = calculateMultiIndices(Nuni,deg)
      dim = size(ind,1)

      new(name,deg,dim,ind,m,uni)
    end
end

struct Tensor
dim::Int64          # "dimension"
T::SparseVector{Float64,Int64}
get::Function
op::Union{OrthoPolyQ,MultiOrthoPoly}
  function Tensor(dim::Int64,mop::MultiOrthoPoly)
    T = computeTensorizedSP(dim,mop)
    g(a) = getentry(a,T,mop.ind,dim)
    new(dim,T,g,mop)
  end
  function Tensor(dim::Int64,opq::OrthoPolyQ)
    T = computeTensorizedSP(dim,opq)
    g(a) = getentry(a,T,calculateMultiIndices(1,opq.op.deg),dim)
    new(dim,T,g,opq)
  end
end

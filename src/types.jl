export  OrthoPoly,
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
    pars::Dict
    function JacobiMeasure(d::Dict)
        @assert d[:shape_a] > -1 && d[:shape_b] > -1 "Invalid shape parameters."
        new(build_w_jacobi(d[:shape_a], d[:shape_b]),(-1,1),d[:shape_a] == d[:shape_b],d)
    end
end

struct LaguerreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function LaguerreMeasure()
        new(w_laguerre,(0,Inf),false)
    end
end

struct genLaguerreMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    pars::Dict
    function genLaguerreMeasure(d::Dict)
        new(build_w_genlaguerre(d[:shape]),(0,Inf),false,d)
    end
end

struct HermiteMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool

    function HermiteMeasure()
        new(w_hermite,(-Inf,Inf),true)
    end
end

struct genHermiteMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    pars::Dict
    function genHermiteMeasure(d::Dict)
        new(build_w_genhermite(Float64(d[:mu])),(-Inf,Inf),true,d)
    end
end

struct MeixnerpollaczekMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    pars::Dict
    function MeixnerpollaczekMeasure(d::Dict)
        new(build_w_meixner_pollaczek(d[:lambda], d[:phi]),(-Inf,Inf),false,d)
    end
end

struct GaussianMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    function GaussianMeasure()
        new(w_gaussian,(-Inf,Inf),true)
    end
end

struct Uniform01Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    function Uniform01Measure()
        new(w_uniform01,(0,1),true)
    end
end

struct Beta01Measure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    pars::Dict
    function Beta01Measure(d::Dict)
        @assert d[:shape_a] > 0 && d[:shape_b] > 0 "Invalid shape parameters."
        new(build_w_beta(d[:shape_a],d[:shape_b]),(0,1),d[:shape_a] == d[:shape_b],d)
    end
end

struct GammaMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    pars::Dict
    function GammaMeasure(d::Dict)
        @assert d[:rate] == 1.
        new(build_w_gamma(d[:shape]),(0, Inf),false,d)
    end
end

struct LogisticMeasure <: AbstractCanonicalMeasure
    w::Function
    dom::Tuple{Real,Real}
    symmetric::Bool
    function LogisticMeasure()
        new(w_logistic,(-Inf,Inf),true)
    end
end


"""    # error handling
    name == "jacobi" && !(d[:shape_a] > -1 && d[:shape_b] > -1) &&  throw(AssertionError("Invalid shape parameters."))
    name == "beta01" && !(d[:shape_a] >  0 && d[:shape_b] >  0) &&  throw(AssertionError("Invalid shape parameters."))
    name == "gamma"  && !(d[:rate] == 1.) &&  error("Rates different from one not supported currently.")
"""

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
        new(lowercase(name),deg,α,β,m)
    end
end

struct LegendreOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function LegendreOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_legendre(Nrec)...)
    end
end

struct JacobiOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function JacobiOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_jacobi(Nrec,Float64(d[:shape_a]), Float64(d[:shape_b]))...)
    end
end

struct LaguerreOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function LaguerreOrthoPoly(deg::Int64,;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_laguerre(Nrec)...)
    end
end

struct genLaguerreOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function genLaguerreOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_laguerre(Nrec,Float64(d[:shape]))...)
    end
end

struct HermiteOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function HermiteOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_hermite(Nrec)...)
    end
end

struct genHermiteOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function genHermiteOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_hermite(Nrec,Float64(d[:mu]))...)
    end
end

struct MeixnerpollaczekOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function MeixnerpollaczekOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,rm_meixner_pollaczek(Nrec,Float64(d[:lambda]),Float64(d[:phi]))...)
    end
end

struct GaussianOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function GaussianOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1/sqrt(2pi),rm_hermite_prob(Nrec)...)...)
    end
end

struct Uniform01OrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function Uniform01OrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1.,rm_legendre01(Nrec)...)...)
    end
end

struct Beta01OrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function Beta01OrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1/beta(d[:shape_a],d[:shape_b]),rm_jacobi01(Nrec,Float64(d[:shape_b])-1,Float64(d[:shape_a])-1)...)...)
    end
end

struct GammaOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function GammaOrthoPoly(deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale((d[:rate]^d[:shape])/gamma(d[:shape]),rm_laguerre(Nrec,Float64(d[:shape]-1.))...)...)
    end
end

struct LogisticOrthoPoly <: AbstractOrthoPoly
    deg::Int64          # maximum degree
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients
    quad::AbstractQuad

    # inner constructor
    function LogisticOrthoPoly(deg::Int64;Nrec::Int64=deg+1)
        @assert deg >= 0 "Degree has to be non-negative."
        @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
        new(deg,r_scale(1.,rm_logistic(Nrec)...)...)
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

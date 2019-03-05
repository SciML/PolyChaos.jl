export  OrthoPoly,
        MultiOrthoPoly,
        Quad,
        OrthoPolyQ,
        Measure,
        MultiMeasure,
        Tensor

struct Measure
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
function Measure(name::String,d::Dict=Dict())
    name = lowercase(name)
    if name == "legendre"
      s = (-1,1); symm = true
      return Measure(name,w_legendre,s,symm,d)
    elseif name == "jacobi"
        s = (-1,1)
        par1, par2 = d[:shape_a], d[:shape_b]
        par1==par2 ? symm=true : symm=false
        return Measure(name,build_w_jacobi(par1,par2),s,symm,d)
    elseif name == "laguerre"
        s = (0,Inf); symm = false
        return Measure(name,w_laguerre,s,symm,d)
    elseif name == "genlaguerre"
        s = (0,Inf); symm = false
        return Measure(name,build_w_genlaguerre(d[:shape]),s,symm,d)
    elseif name == "hermite"
        s = (-Inf,Inf); symm = true
        return Measure(name,w_hermite,s,symm,d)
    elseif name == "genhermite"
        s = (-Inf,Inf); symm = true
        return Measure(name,build_w_genhermite(Float64(d[:mu])),s,symm,d)
    elseif name == "meixnerpollaczek"
        s = (-Inf,Inf); symm = false
        par1, par2 = d[:lambda], d[:phi]
        return Measure(name,build_w_meixner_pollaczek(par1,par2),s,symm,d)
  ########################################################
  ########################################################
  ########################################################
  # measures corresponding to probability density functions:
    elseif name == "gaussian"
        s = (-Inf, Inf)
        return Measure(name,w_gaussian,s,true,d)
    elseif name == "uniform01"
        s = (0,1)
        return Measure(name,w_uniform01,s,true,d)
    elseif name == "beta01"  # parameters of beta distribution
        s = (0,1)
        par1, par2 = d[:shape_a], d[:shape_b]
        @assert par1 > 0 && par2 > 0 "Invalid shape parameters."
        symm = (par1 == par2)
    return Measure(name,build_w_beta(par1,par2),s,symm,d)
        elseif name == "gamma"
        shape, rate = d[:shape], d[:rate]
        @assert rate==1. "rates different from one not yet supported."
        s = (0,Inf)
        return Measure(name,build_w_gamma(shape),s,false,d)
    elseif name == "logistic"
        s = (-Inf,Inf)
        return Measure(name,w_logistic,s,true,d)
    else
        error("Measure named `$name` is not yet implemented.")
    end
end

struct OrthoPoly
  name::String
  deg::Int64          # maximum degree
  α::Vector{Float64}  # recurrence coefficients
  β::Vector{Float64}  # recurrence coefficients
  meas::Measure
  # inner constructor
  function OrthoPoly(name::String,deg::Int64,α::Vector{Float64},β::Vector{Float64},m::Measure)
    @assert deg >= 0 "Degree has to be non-negative."
    @assert length(α)==length(β) "Different number of recursion coefficients α and β supplied."
    new(lowercase(name),deg,α,β,m)
  end
end

# constructor for classic distributions
function OrthoPoly(name::String,deg::Int64,d::Dict=Dict();Nrec::Int64=deg+1)
  @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
  name = lowercase(name)
  if name == "legendre"
    d = Dict()
    a,b = rm_legendre(Nrec)
  elseif name == "jacobi"
    par1, par2 = d[:shape_a], d[:shape_b]
    symm = (par1 == par2)
    a,b = rm_jacobi(Nrec,Float64(par1),Float64(par2))
  elseif name == "laguerre"
    d = Dict()
    a,b = rm_laguerre(Nrec)
  elseif name == "genlaguerre"
    p = d[:shape]
    a,b = rm_laguerre(Nrec,Float64(p))
  elseif name == "hermite"
    d = Dict()
    a,b = rm_hermite(Nrec)
  elseif name == "genhermite"
    p = d[:mu]
    a,b = rm_hermite(Nrec,Float64(p))
  elseif name == "meixnerpollaczek"
    par1, par2 = d[:lambda], d[:phi]
    a,b = rm_meixner_pollaczek(Nrec,Float64(par1),Float64(par2))
  ########################################################
  ########################################################
  ########################################################
  # measures corresponding to probability density functions:
  elseif name == "gaussian"
    d = Dict()
    a,b = r_scale(1/sqrt(2pi),rm_hermite_prob(Nrec)...)
  elseif name == "uniform01"
    d = Dict()
    a,b = r_scale(1.,rm_legendre01(Nrec)...)
  elseif name == "beta01"  # parameters of beta distribution
    par1, par2 = d[:shape_a], d[:shape_b]
    a,b = r_scale(1/beta(par1,par2),rm_jacobi01(Nrec,Float64(par2)-1,Float64(par1)-1)...)
  elseif name == "gamma"
    shape, rate = d[:shape], d[:rate]
    @assert rate==1. "rates different from one not yet supported."
    a,b = r_scale((rate^shape)/gamma(shape),rm_laguerre(Nrec,Float64(shape-1.))...)
  elseif name == "logistic"
    d = Dict()
    a,b = r_scale(1.,rm_logistic(Nrec)...)
  else
    error("$name is not yet implemented.")
  end
  m = Measure(name,d)
  return OrthoPoly(name,deg,a,b,m)
end

# general constructor
function OrthoPoly(name::String,deg::Int64,m::Measure;Nrec=deg+1,Nquad=10*Nrec,quadrature::Function=clenshaw_curtis,discretization::Function=stieltjes)
  @assert Nrec >= deg + 1 "Not enough recurrence coefficients specified"
  name = lowercase(name)
  a,b = rm_compute(m,Nrec,Nquad,quadrature=quadrature,discretization=discretization)
  return OrthoPoly(name,deg,a,b,m)
end

function OrthoPoly(name::String,deg::Int64,w::Function,s::Tuple{Real,Real},symm::Bool,d::Dict=Dict();Nrec=deg+1,Nquad=10*Nrec,quadrature::Function=clenshaw_curtis,discretization::Function=stieltjes)
  @assert Nrec>=deg+1 "Not enough recurrence coefficients specified"
  name = lowercase(name)
  m = Measure(name,w,s,symm,d)
  a,b = rm_compute(m,Nrec,Nquad,quadrature=quadrature,discretization=discretization)
  return OrthoPoly(name,deg,a,b,m)
end

struct Quad
  name::String       # name of quadrature
  Nquad::Int64          #   number of qudrature points
  nodes::Vector{Float64}
  weights::Vector{Float64}
  meas::Measure
  function Quad(name::String,N::Int64,nodes::Vector{Float64},weights::Vector{Float64},m::Measure)
    @assert N >= 1 "Number of qudrature points has to be positive"
    @assert length(nodes) == length(weights) "Inconsistent number of nodes and weights inconsistent."
    new(lowercase(name),N,nodes,weights,m)
  end
end

function Quad(N::Int64,m::Measure;quadrature::Function=clenshaw_curtis)
  n,w = quadgp(m.w,m.dom[1],m.dom[2],N;quadrature=quadrature)
  Quad("quadgp",N,n,w,m)
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
function Quad(N::Int64,weight::Function,supp::Tuple{Real,Real},symm::Bool,d::Dict=Dict();quadrature::Function=clenshaw_curtis)
    @assert N >= 1 "Number of qudrature points has to be positive"
    m = Measure("fun_"*String(nameof(weight)),weight,supp,symm,d)
    n,w = quadgp(weight,Float64(supp[1]),Float64(supp[2]),N;quadrature=quadrature)
    Quad("quadgp",N,n,w,m)
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
    function MultiOrthoPoly(uni::Union{Vector{OrthoPoly}},deg::Int64)
      @assert deg>=0 "degree has to be non-negative"
      degs = [ op.deg for op in uni ]
      @assert deg<=minimum(degs) "Requested degree $deg is greater than smallest univariate degree $(minimum(degs))."
      Nuni::Int64 = length(uni)

      name_meas = [ op.meas.name for op in uni      ]
      supp      = [ op.meas.dom for op in uni       ]
      symm      = [ op.meas.symmetric for op in uni ]
      pars      = [ op.meas.pars for op in uni      ]
      w_uni     = [ op.meas.w for op in uni         ]
      w(t) = prod([uni[i].meas.w(t[i]) for i=1:Nuni])
      m = MultiMeasure(name_meas,w,w_uni,supp,symm,pars)

      name = [ op.name for op in uni ]
      ind = calculateMultiIndices(Nuni,deg)
      dim = size(ind,1)

      new(name,deg,dim,ind,m,uni)
    end
    function MultiOrthoPoly(uni::Union{Vector{OrthoPolyQ}},deg::Int64)
      @assert deg >= 0 "degree has to be non-negative"
      degs = [ s.op.deg for s in uni ]
      @assert deg<=minimum(degs) "Requested degree $deg is greater than smallest univariate degree $(minimum(degs))."
      Nuni::Int64 = length(uni)

      name_meas = [ s.op.meas.name for s in uni         ]
      supp      = [ s.op.meas.dom for s in uni          ]
      symm      = [ s.op.meas.symmetric for s in uni    ]
      pars      = [ s.op.meas.pars for s in uni         ]
      w_uni     = [ s.op.meas.w for s in uni            ]
      w(t) = prod([uni[i].op.meas.w(t[i]) for i=1:Nuni])
      m = MultiMeasure(name_meas,w,w_uni,supp,symm,pars)

      name = [ s.op.name for s in uni ]
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

# struct PolyChaosExpansion
#   name::String
#   coeffs::Vector{Float64}
#   dim::Int64
#   op::Union{OrthoPoly,OrthoPolyQ,MultiOrthoPoly}
#   function PolyChaosExpansion(n::String,c::Vector{Float64},d::Int64,op::Union{OrthoPoly,OrthoPolyQ,MultiOrthoPoly})
#     @assert length(c)==d "length of coefficients inconsistent with dimension ($(length(c)) vs. $d)"
#     new(n,c,d,op)
#   end
# end
#
# function PolyChaosExpansion(name)

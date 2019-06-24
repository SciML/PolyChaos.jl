export  convert2affinePCE,
        calculateAffinePCE,
        assign2multi,
        evaluatePCE,
        sampleMeasure,
        samplePCE,
        mean,
        var,
        std

"""
    calculateAffinePCE(α::Vector{Float64})::Vector{Float64}
Computes the affine PCE coefficients ``x_0`` and ``x_1`` from recurrence coefficients ``\alpha``.
"""
function calculateAffinePCE(α::Vector{Float64})
    @assert length(α)>=1 "not enough recursion coefficients"
    return [α[1], 1.]
end

calculateAffinePCE(op::OrthoPoly) = calculateAffinePCE(op.α)
calculateAffinePCE(opq::OrthoPolyQ) = calculateAffinePCE(opq.op)

function calculateAffinePCE(i::Int64,mop::MultiOrthoPoly)
    p = length(mop.name)
    @assert i<=p "mop is $p-variate, but the PC for the $i-variate entry was requested"
    calculateAffinePCE(mop.uni[i])
end

"""
```
convert2affinePCE(a::Vector{Float64},α0::Float64)
convert2affinePCE(name::String,p1::Float64,p2::Float64,d::Dict=Dict();kind::Symbol=:lbub)
```
Computes the affine PCE coefficients ``x_0`` and ``x_1`` from
```math
X = a_1 + a_2 \\Xi = x_0 + x_1 \\phi_1(\\Xi),
```
where ``\\phi_1(t) = t-\\alpha_0`` is the first-order monic basis polynomial.

For classical polynomials the `name` can be given directly. The keyword `kind in [:lbub, :μσ]`
specifies whether `p1` and `p2` have the meaning of lower/upper bounds or mean/standard deviation.
"""
function convert2affinePCE(a1::Float64,a2::Float64,α0::Float64)
    [
        a1+α0*a2;
        a2
    ]
end
function convert2affinePCE(name::String,p1::Float64,p2::Float64,d::Dict=Dict();kind::Symbol=:lbub)
    @assert kind in [:lbub, :μσ]
    name = lowercase(name)
    a1, a2, α = 0., 0., zeros(Float64,1)
    name=="beta01" && d[:shape_a]==d[:shape_b]==1 ? name="uniform01" : ()

    if name=="gaussian"
        # p1 --> mean
        # p2 --> standard deviation
        @assert p2>=0. "standard deviation has to be non-negative"
        abs(p2)<1e-4 ? (@warn "standard deviation close to zero: $p2") : ()
        a1, a2 = p1, p2
        α,~ = rm_hermite_prob(1)
    elseif name=="uniform01"
        α,~ = rm_legendre01(1)
        if kind==:lbub
            # p1 --> lower bound
            # p2 --> upper bound
            @assert p1<p2 "inconsistent bounds"
            a1, a2 = p1, p2-p1
        elseif kind==:μσ
            # p1 --> mean
            # p2 --> standard deviation
            @assert p2>=0. "standard deviation has to be nonnegative"
            abs(p2)<1e-4 ? (@warn "standard deviation close to zero: $p2") : ()
            a1, a2 = p1-sqrt(3)*p2 , 2*sqrt(3)*p2
        else
            error("not implemented")
        end
    elseif name=="beta01"
        s_α, s_β = d[:shape_a], d[:shape_b]
        α,~ = rm_jacobi01(1,Float64(s_β)-1,Float64(s_α)-1)
        if kind==:lbub
            # p1 --> lower bound
            # p2 --> upper bound
            @assert p1<p2 "inconsistent bounds"
            a1, a2 = p1, p2-p1
        elseif kind==:μσ
            # p1 --> mean
            # p2 --> standard deviation
            @assert p2>=0. "standard deviation has to be nonnegative"
            abs(p2)<1e-4 ? (@warn "standard deviation close to zero: $p2") : ()
            a1 = p1-sqrt(s_α/s_β)*sqrt(1+s_α+s_β)*p2
            a2 = (s_α+s_β)*sqrt((s_α+s_β+1)/(s_α*s_β))*p2
        else
            error("not implemented")
        end
    elseif name=="gamma"
        error("$name is not implemented")
    elseif name=="logistic"
        @assert p2>0
        α,~ = rm_logistic(1)
        a1, a2 = p1, p2
    else
        error("$name is not implemented")
    end
    convert2affinePCE(a1,a2,α[1])
end
# convert2affinePCE(name::String,p1::Float64,p2::Float64;kind::Symbol=:lbub) = convert2affinePCE(name,p1,p2,Dict();kind=kind)

function assign2multi(x::Vector{Float64},i::Int64,ind::Matrix{Int64})
    l,p = size(ind)
    nx, deg = length(x), ind[end,end]
    @assert nx<=deg+1 "inconsistent number of coefficients ($nx vs $(deg+1))"
    @assert i<=p "basis is $p-variate, you requested $i-variate"
    myind::Vector{Int64} = findUnivariateIndices(i,ind)[1:nx]
    y = spzeros(Float64,l)
    y[myind] = x
    return y
end

function sampleMeasure_byName(n::Int64,name::String,d::Dict=Dict())::Vector{Float64}
    @assert n >= 1 "invalid number $n of samples."
    name = lowercase(name)
    name == "gaussian" && return rand(Normal(),n)
    name == "uniform01" && return rand(Uniform(),n)
    name == "beta01" && return rand(Beta(d[:shape_a],d[:shape_b]),n)
    name == "gamma" && return rand(Gamma(d[:shape],1/d[:rate]),n)
    name == "logistic" && return rand(Logistic(),n)
    throw(error("$name is not implemented"))
end

function sampleMeasure_byFun(n::Int64,w::Function,dom::Tuple{Float64,Float64};method::String="adaptiverejection")
    @assert n>=1 "invalid number $n of samples."
    method = lowercase(method)
    @assert method in [ "adaptiverejection", "inversecdf"] "method $method not implemented"
    if method == "adaptiverejection"
        # works only if w is log-concave
        sampler = RejectionSampler(w,dom)
        return run_sampler!(sampler,n)
    elseif method == "rejection"
        # all purpose method but needs a solid envelope PDF
        error("method $method not yet implemented")
    elseif method == "inversecdf"
        error("method $method not yet implemented")
        # requires CDF and reliable root-finding method
    end
end

"""
__Univariate__
```
sampleMeasure(n::Int64,name::String,w::Function,dom::Tuple{Float64,Float64},symm::Bool,d::Dict;method::String="adaptiverejection")
sampleMeasure(n::Int64,m::Measure;method::String="adaptiverejection")
sampleMeasure(n::Int64,op::OrthoPoly;method::String="adaptiverejection")
sampleMeasure(n::Int64,opq::OrthoPolyQ;method::String="adaptiverejection")
```
Draw `n` samples from the measure `m` described by its
- `name`
- weight function `w`,
- domain `dom`,
- symmetry property `symm`,
- and, if applicable, parameters stored in the dictionary `d`.
By default an adaptive rejection sampling method is used (from [AdaptiveRejectionSampling.jl](https://github.com/mauriciogtec/AdaptiveRejectionSampling.jl)),
unless it is a common random variable for which [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) is used.

The function is multiply dispatched to accept `OrthoPoly` or `OrthoPolyQ`.

__Multivariate__
```
sampleMeasure(n::Int64,m::MultiMeasure;method::Vector{String}=["adaptiverejection" for i=1:length(m.name)])
sampleMeasure(n::Int64,mop::MultiOrthoPoly;method::Vector{String}=["adaptiverejection" for i=1:length(mop.meas.name)])
```
Multivariate extension which provides array of samples with `n` rows and
as many columns as the multimeasure has univariate measures.
"""
function sampleMeasure(n::Int64,name::String,w::Function,dom::Tuple{Float64,Float64},symm::Bool,d::Dict;method::String="adaptiverejection")
    # ToDo: symmetry is currently not exploited, but could be for rejection sampling
    name = lowercase(name)
    if name in ["gaussian","uniform01","beta01","gamma","logistic"]
        sampleMeasure_byName(n,name,d)
    else
        sampleMeasure_byFun(n,w,dom;method=method)
    end
end
sampleMeasure(n::Int64,m::Measure;method::String="adaptiverejection") = sampleMeasure(n,m.name,m.w,m.dom,m.symmetric,m.pars;method=method)
sampleMeasure(n::Int64,op::OrthoPoly;method::String="adaptiverejection") = sampleMeasure(n,op.meas;method=method)
sampleMeasure(n::Int64,opq::OrthoPolyQ;method::String="adaptiverejection") = sampleMeasure(n::Int64,opq.op;method=method)

function sampleMeasure(n::Int64, name::Vector{String}, w::Vector{Function},
                       dom::Vector{Tuple{Float64,Float64}}, symm::BitArray,d::Vector{Dict};
                       method::Vector{String}=["adaptiverejection" for i=1:length(name)])
    @assert length(name)==length(w)==length(dom)==length(symm)==length(d) "inconsistent number of parameters"
    m = length(name)
    ξ = zeros(n,m)
    for k=1:m
        ξ[:,k] = sampleMeasure(n,name[k],w[k],dom[k],symm[k],d[k];method=method[k])
    end
    return ξ
end
function sampleMeasure(n::Int64,m::MultiMeasure;method::Vector{String}=["adaptiverejection" for i=1:length(m.name)])
    s = falses(length(m.name))
    s[:] = m.symmetric[:]
    sampleMeasure(n,m.name,m.w_uni,m.dom,s,m.pars;method=method)
end
sampleMeasure(n::Int64,mop::MultiOrthoPoly;method::Vector{String}=["adaptiverejection" for i=1:length(mop.meas.name)]) = sampleMeasure(n,mop.meas;method=method)

"""
    evaluatePCE(x::Vector{Float64},ξ::Vector{Float64},α::Vector{Float64},β::Vector{Float64})
Evaluation of polynomial chaos expansion
```math
\\mathsf{x} = \\sum_{i=0}^{L} x_i \\phi_i{\\xi_j},
```
where `L+1 = length(x)` and ``x_j`` is the ``j``th sample where ``j=1,\\dots,m``
with `m = length(ξ)`.
"""
function evaluatePCE(x::Vector{Float64},ξ::Vector{Float64},α::Vector{Float64},β::Vector{Float64})
    @assert length(α) == length(β) "inconsistent number of recurrence coefficients"
    Nsmpl = length(ξ)
    @assert Nsmpl>=1 "inconsistent number $Nsmpl of samples"
    Nx, Nrec = length(x), length(α)
    @assert Nx<=Nrec "not enough recursion coefficients"
    if Nrec > Nx α,β = α[1:Nx], β[1:Nx] end
    ϕ = zeros(Float64,Nsmpl,Nx)
    for n in 1:Nx
        ϕ[:,n] = evaluate(n-1,ξ,α,β)
    end
    ϕ*x
end
evaluatePCE(x::Vector{Float64},ξ::Float64,α::Vector{Float64},β::Vector{Float64}) = evaluatePCE(x,[ξ],α,β)
evaluatePCE(x::Vector{Float64},ξ::Vector{Float64},op::OrthoPoly) = evaluatePCE(x,ξ,op.α,op.β)
evaluatePCE(x::Vector{Float64},ξ::Vector{Float64},opq::OrthoPolyQ) = evaluatePCE(x,ξ,opq.op)

function evaluatePCE(x::Vector{Float64},ξ::Matrix{Float64},α::Vector{Vector{Float64}},β::Vector{Vector{Float64}},ind::Matrix{Int64})
    Nsmpl = size(ξ,1)
    @assert Nsmpl>=1 "inconsistent number $Nsmpl of samples"
    @assert length(α)==length(β)==size(ξ,2)==size(ind,2) "inconsistent number of coefficients"
    Nx = length(x)
    @assert Nx <= size(ind,1) "too few pc coefficients (resp: too small basis)"
    ϕ = zeros(Float64,Nsmpl,Nx)
    for n in 1:Nx
        ϕ[:,n] = evaluate(ind[n,:],ξ,α,β)
    end
    ϕ*x
end
function evaluatePCE(x::Vector{Float64},ξ::Matrix{Float64},mOP::MultiOrthoPoly)
    a,b = coeffs(mOP)
    evaluatePCE(x,ξ,a,b,mOP.ind)
end

"""
__Univariate__
```
samplePCE(n::Int64,x::Vector{Float64},op::OrthoPoly;method::String="adaptiverejection")
samplePCE(n::Int64,x::Vector{Float64},opq::OrthoPolyQ;method::String="adaptiverejection")
```
Combines [`sampleMeasure`](@ref) and [`evaluatePCE`](@ref), i.e. it first draws `n` samples
from the measure, then evaluates the PCE for those samples.

__Multivariate__
```
samplePCE(n::Int64,x::Vector{Float64},mop::MultiOrthoPoly;method::Vector{String}=["adaptiverejection" for i=1:length(mop.meas.name)])
```
"""
function samplePCE(n::Int64,x::Vector{Float64},op::OrthoPoly;method::String="adaptiverejection")
    ξ = sampleMeasure(n,op;method=method)
    evaluatePCE(x,ξ,op)
end

samplePCE(n::Int64,x::Vector{Float64},opq::OrthoPolyQ;method::String="adaptiverejection") = samplePCE(n,x,opq.op;method=method)

function samplePCE(n::Int64,x::Vector{Float64},mop::MultiOrthoPoly;method::Vector{String}=["adaptiverejection" for i=1:length(mop.meas.name)])
    ξ = sampleMeasure(n,mop;method=method)
    evaluatePCE(x,ξ,mop)
end

"""
__Univariate__
```
mean(x::Vector{},op::OrthoPoly)
mean(x::Vector{},opq::OrthoPolyQ)
```
__Multivariate__
```
mean(x::Vector{},mop::MultiOrthoPoly)
```
compute mean of random variable with PCE `x`
"""
function mean(x::Vector{},op::OrthoPoly)
    x[1]*computeSP2(0,op.β)
end

mean(x::Vector{},opq::OrthoPolyQ) = mean(x,opq.op)

function mean(x::Vector{},mop::MultiOrthoPoly)
    nunc = length(mop.uni)
    x[1]*computeSP(zeros(Int64,nunc),mop)
end

"""
__Univariate__
```
var(x::Vector{},op::OrthoPoly)
var(x::Vector{},opq::OrthoPolyQ)
var(x::Vector{},t2::Tensor)
```
__Multivariate__
```
var(x::Vector{},mop::MultiOrthoPoly)
var(x::Vector{},t2::Tensor)
```
compute variance of random variable with PCE `x`
"""
function var(x::Vector{},op::OrthoPoly)
    t2 = computeSP2(op)
    @assert length(t2)>=length(x) "cannot compute variance; too many PCE coefficients"
    sum( x[i]^2*t2[i] for i=2:length(x) )
end
var(x::Vector{},opq::OrthoPolyQ) = var(x,opq.op)

function var(x::Vector{},t2::Tensor)
    sum( x[i]^2*t2.get([i-1,i-1]) for i in 2:length(x) )
end

function var(x::Vector{},mop::MultiOrthoPoly)
    @assert length(x) <= size(mop.ind,1) "cannot compute variance; too many PCE coefficients"
    var(x,Tensor(2,mop))
end

"""
__Univariate__
```
std(x::Vector{},op::OrthoPoly)
std(x::Vector{},opq::OrthoPolyQ)
```
__Multivariate__
```
std(x::Vector{},mop::MultiOrthoPoly)
```
compute standard deviation of random variable with PCE `x`
"""
function std(x::Vector{},op::OrthoPoly)
    sqrt(var(x,op))
end
std(x::Vector{},opq::OrthoPolyQ) = std(x,opq.op)
std(x::Vector{},mop::MultiOrthoPoly) = sqrt(var(x,mop))

function moment(n::Int64,x::Vector{},op::OrthoPoly)
    @assert 1<=n<=2 "for moments greater than 2 a quadrature rule is required"
    n==1 ? mean(x,op) : (mean(x,op)^2+var(x,op))
end

# function moment(n::Int64,x::Vector{Float64},opq::OrthoPolyQ)
#     @assert n>=1 "invalid moment"
#     if n==1
#         return mean(x,opq)
#     elseif n==2
#         return mean(x,opq)^2+var(x,opq)
#     else
#         t = Tensor(n,opq)
#         # needs IterTools for an elegant solution
#         error("This is not yet implemented")
#     end
# end

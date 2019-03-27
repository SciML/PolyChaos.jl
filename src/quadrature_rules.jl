
export  fejer,
        fejer2,
        clenshaw_curtis,
        quadgp,
        golubwelsch,
        gauss,
        radau,
        lobatto
"""
    fejer(N::Int64)
Fejer's first quadrature rule.
"""
function fejer(N::Int64)
    @assert N >= 1 "N has to be positive"
    N == 1 && return zeros(1), [2.]
    θ = map(x->(2x-1)*pi/(2N),1:N)
    M = N ÷ 2
    return cos.(θ), map(x->2/N*(1-2*sum(cos(2*n*x)/(4*n^2-1) for n=1:M)),θ)
end

"""
    fejer2(n::Int64)
Fejer's second quadrature rule according to [Waldvogel, J. Bit Numer Math (2006) 46: 195](https://doi.org/10.1007/s10543-006-0045-4).
"""
function fejer2(n::Int64)
    @assert n >= 2
    N = 1:2:n-1
    m = n - length(N)
    v0 = push!(push!(map(x->2/x/(x-2.),N),1/last(N)), zeros(m)...)
    @inbounds v2 = -v0[1:n] - v0[n+1:-1:2];
    wf2 = real(ifft(v2))
    return map(x->cos(x*pi/n),0:n), push!(wf2, first(wf2))
end

"""
    clenshaw_curtis(n::Int64)
Clenshaw-Curtis quadrature according to [Waldvogel, J. Bit Numer Math (2006) 46: 195](https://doi.org/10.1007/s10543-006-0045-4).
"""
function clenshaw_curtis(n::Int64)::Tuple{Vector{Float64},Vector{Float64}}
    @assert n >= 2
    N = 1:2:n-1
    l = length(N)
    m = n - l
    v0 = push!(push!(map(x->2/x/(x-2.),N),1/last(N)), zeros(m)...)
    @inbounds v2 = -v0[1:n] - v0[n+1:-1:2]
    g0= -ones(Float64,n)
    g0[1+l] += n
    g0[1+m] += n
    g = 1 / (n^2 - 1 + mod(n,2) ) * g0
    wcc = real(ifft(v2 + g));
    return map(i->cos(i/n*pi),0:n), push!(wcc,first(wcc))
end

"""
    quadgp(weight::Function,lb::Float64,ub::Float64,N::Int64=10;quadrature::Function=clenshaw_curtis,bnd::Float64=Inf)
general purpose quadrature based on Gautschi, "Orthogonal Polynomials: Computation and Approximation", Section 2.2.2, pp. 93-95

Compute the `N`-point quadrature rule for `weight` with support (`lb`, `ub`).
The quadrature rule can be specified by the keyword `quadrature`.
The keyword `bnd` sets the numerical value for infinity.
"""
function quadgp(weight::Function,lb::Float64,ub::Float64,N::Int64=10;quadrature::Function=clenshaw_curtis,bnd::Float64=Inf)
    @assert lb < ub "inconsistent interval bounds"
    t_fej,w_fej = quadrature(N)
    t, w, d = zeros(Float64,N),zeros(Float64,N),zeros(Float64,N)
    # transformation
    if (-bnd<lb<bnd) && (-bnd<ub<bnd)
        t = 0.5*((ub-lb)*t_fej .+ (ub+lb))
        d = 0.5*(ub-lb)*ones(Float64,size(t_fej))
    elseif -bnd>=lb && (-bnd<ub<bnd)
        t = ub .- (1 .- t_fej)./(1 .+ t_fej)
        d = 2 ./(1 .+ t_fej).^2
    elseif (-bnd<lb<bnd) && ub>=bnd
        t = lb .+ (1 .+ t_fej)./(1 .- t_fej)
        d = 2 ./(1 .- t_fej).^2
    elseif -bnd>=lb && ub>=bnd
        t = t_fej./(1 .- t_fej.^2)
        d = (1 .+ t_fej.^2)./( 1 .- t_fej.^2).^2
    end
    w = w_fej.*weight.(t).*d
    return t, w
end

####################################################
####################################################
####################################################
# Quadrature rules based on recurrence coefficients

function golubwelsch(α::Vector{Float64},β::Vector{Float64},maxiter::Int=30)
    N = length(α) - 1
    a, β0 = copy(α[1:N]), β[1]
    w = zero(a)
    special_eigenproblem!(a,copy(sqrt.(β)),w,maxiter)
    idx = sortperm(a)
    a[idx], (β0*w.^2)[idx]
end

golubwelsch(op::OrthoPoly) = golubwelsch(op.α,op.β)

"""
    gauss(N::Int64,α::Vector{Float64},β::Vector{Float64})
    gauss(α::Vector{Float64},β::Vector{Float64})
    gauss(N::Int64,op::OrthoPoly)
    gauss(op::OrthoPoly)
Gauss quadrature rule, also known as Golub-Welsch algorithm

`gauss()` generates the `N` Gauss quadrature nodes and weights for a given weight function.
The weight function is represented by the `N` recurrence coefficients for the monic polynomials orthogonal
with respect to the weight function.

!!! note
    The function `gauss` accepts at most `N = length(α) - 1` quadrature points,
        hence providing at most an `(length(α) - 1)`-point quadrature rule.

!!! note
    If no `N` is provided, then `N = length(α) - 1`.
"""
function gauss(N::Int64,α::Vector{Float64},β::Vector{Float64})
    N += 1
    @assert N > 0 "only positive N allowed"
    @assert length(α) == length(β) "inconsistent number of recurrence coefficients"
    N0 = length(α)
    @assert N0 >= N  "not enough recurrence coefficients"
    @inbounds golubwelsch(α[1:N],β[1:N])
end
gauss(α::Vector{Float64},β::Vector{Float64}) = gauss(length(α)-1,α,β)
gauss(N::Int64,op::OrthoPoly) = gauss(N::Int64,op.α,op.β)
gauss(op::OrthoPoly) = gauss(op.α,op.β)

"""
    radau(N::Int64,α::Vector{Float64},β::Vector{Float64},end0::Float64)
    radau(α::Vector{Float64},β::Vector{Float64},end0::Float64)
    radau(N::Int64,op::OrthoPoly,end0::Float64)
    radau(op::OrthoPoly,end0::Float64)
Gauss-Radau quadrature rule.
Given a weight function encoded by the recurrence coefficients `(α,β)`for the associated
orthogonal polynomials, the function generates the
nodes and weights `(N+1)`-point Gauss-Radau
quadrature rule for the weight function having a prescribed
node `end0` (typically at one of the end points of the support
interval of w, or outside thereof).

!!! note
    The function `radau` accepts at most `N = length(α) - 2` as an input,
        hence providing at most an `(length(α) - 1)`-point quadrature rule.

!!! note
    Reference: OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS AND RELATED QUADRATURE RULES by Walter Gautschi
"""
function radau(N::Int64,α_::Vector{Float64},β::Vector{Float64},end0::Float64)
    α = copy(α_)
    @assert N > 0 "only positive N allowed"
    @assert length(α) == length(β) > 0 "inconsistent number of recurrence coefficients"
    @assert length(α) >= N + 1 "not enough recurrence coefficients"
    p0 = 0.
    p1 = 1.
    for n in Base.OneTo(N)
      pm1 = p0
      p0 = p1;
      @inbounds p1 = (end0 - α[n])*p0 - β[n]*pm1;
    end
    @inbounds α[N+1] = end0 - β[N+1]*p0/p1
    gauss(N+1,α,β)
end
radau(α::Vector{Float64},β::Vector{Float64},end0::Float64) = radau(length(α)-2,α,β,end0)
radau(N::Int64,op::OrthoPoly,end0::Float64) = radau(N,op.α,op.β,end0)
radau(op::OrthoPoly,end0::Float64) = radau(op.α,op.β,end0::Float64)

"""
    lobatto(N::Int64,α::Vector{Float64},β::Vector{Float64},endl::Float64,endr::Float64)
    lobatto(α::Vector{Float64},β::Vector{Float64},endl::Float64,endr::Float64)
    lobatto(N::Int64,op::OrthoPoly,endl::Float64,endr::Float64)
    lobatto(op::OrthoPoly,endl::Float64,endr::Float64)
Gauss-Lobatto quadrature rule.
Given a weight function encoded by the recurrence coefficients for the associated
orthogonal polynomials, the function generates
the nodes and weights of the `(N+2)`-point Gauss-Lobatto
quadrature rule for the weight function, having two
prescribed nodes `endl`, `endr` (typically the left and right end
points of the support interval, or points to the left
resp. to the right thereof).

!!! note
    The function `radau` accepts at most `N = length(α) - 3` as an input,
    hence providing at most an `(length(α) - 1)`-point quadrature rule.

!!! note
    Reference: OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS AND RELATED QUADRATURE RULES by Walter Gautschi
"""
function lobatto(N::Int64,α_::Vector{Float64},β_::Vector{Float64},endl::Float64,endr::Float64)
    α, β = copy(α_), copy(β_)
    @assert N > 0 "only positive N allowed"
    @assert length(α) == length(β) > 0 "inconsistent number of recurrence coefficients"
    @assert endl < endr "inconsistent end points"
    @assert length(α) >= N+2 "not enough recurrence coefficients"
    p0l = 0.; p0r = 0.; p1l = 1.; p1r = 1.;
    for n=1:N+1
        pm1l=p0l
        p0l=p1l
        pm1r=p0r
        p0r=p1r
        @inbounds p1l=(endl-α[n])*p0l-β[n]*pm1l
        @inbounds p1r=(endr-α[n])*p0r-β[n]*pm1r
    end
    det = p1l * p0r - p1r * p0l
    @inbounds α[N+2] = (endl * p1l * p0r - endr * p1r * p0l) / det
    @inbounds β[N+2] = (endr - endl) * p1l * p1r / det;
    gauss(N+2,α,β)
end
lobatto(α::Vector{Float64},β::Vector{Float64},endl::Float64,endr::Float64) = lobatto(length(α)-3,α,β,endl,endr)
lobatto(N::Int64,op::OrthoPoly,endl::Float64,endr::Float64) = lobatto(N,op.α,op.β,endl,endr)
lobatto(op::OrthoPoly,endl::Float64,endr::Float64) = lobatto(op.α,op.β,endl,endr)

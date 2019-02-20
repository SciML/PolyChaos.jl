export  fejer,
        fejer2,
        clenshaw_curtis,
        quadgp,
        golubwelsch,
        gauss,
        radau,
        radau_jacobi,
        radau_laguerre,
        lobatto,
        lobatto_jacobi
"""
    fejer(N::Int64)
Fejer's first quadrature rule.
"""
function fejer(N::Int64)
    @assert N >= 1 "N has to be positive"
    θ = map(x->(2x-1)*pi/(2N),1:N)
    M = N ÷ 2
    return cos.(θ), map(x->2/N*(1-2*sum(cos(2*n*x)/(4*n^2-1) for n=1:M)),θ)
end

"""
    fejer2(n::Int64)
Fejer's second quadrature rule according to [Waldvogel, J. Bit Numer Math (2006) 46: 195](https://doi.org/10.1007/s10543-006-0045-4).
"""
function fejer2(n::Int64)
    N = 1:2:n-1
    m = n - length(N)
    v0 = push!(push!(map(x->2/x/(x-2.),N),1/(n-1)), zeros(m)...)
    @inbounds v2 = -v0[1:n] - v0[n+1:-1:2];
    wf2 = real(ifft(v2))
    return map(x->cos(x*pi/n),0:n), push!(wf2, first(wf2))
end

"""
    clenshaw_curtis(n::Int64)
Clenshaw-Curtis quadrature according to [Waldvogel, J. Bit Numer Math (2006) 46: 195](https://doi.org/10.1007/s10543-006-0045-4).
"""
function clenshaw_curtis(n::Int64)::Tuple{Vector{Float64},Vector{Float64}}
    @assert n>=1 "n must be > 1"
    x::Vector{Float64} = cos.(collect(0:n)/n*pi)
    N::Vector{Int64}=collect(1:2:n-1);
    l=length(N);
    m=n-l;
    K::Vector{Int64}=collect(0:m-1);
    v0::Vector{Float64}=[2 ./N./(N .- 2);  1/N[end]; zeros(m)];
    v2::Vector{Float64}=-v0[1:end-1]-v0[end:-1:2];
    g0::Vector{Float64}=-ones(n);
    g0[1+l]=g0[1+l]+n;
    g0[1+m]=g0[1+m]+n;
    g::Vector{Float64}=g0/(n^2-1+mod(n,2));
    wcc=real(ifft(v2+g));
    return x, [wcc; wcc[1]]
end

"""
    quadgp(weight::Function,a::Float64,b::Float64,N::Int64=10;quadrature::Function=clenshaw_curti
general purpose quadrature based on Gautschi, "Orthogonal Polynomials: Computation and Approximation", Section 2.2.2, pp. 93-95
"""
function quadgp(weight::Function,lb::Float64,ub::Float64,N::Int64=10;quadrature::Function=clenshaw_curtis,bnd::Float64=Inf)
    @assert lb < ub "inconsistent interval bounds"
    t_fej::Vector{Float64},w_fej::Vector{Float64} = quadrature(N)
    t::Vector{Float64},w::Vector{Float64},d::Vector{Float64} = zeros(Float64,N),zeros(Float64,N),zeros(Float64,N)
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
    size(w)
    return t, w
end

####################################################
####################################################
####################################################
# Quadrature rules based on recurrence coefficients

function golubwelsch(α::Vector{Float64},β::Vector{Float64})
    J = SymTridiagonal( α, sqrt.(β[2:end]) )
    nodes, V = eigen(J)
    # nodes, V = e.values, e.vectors
    weights = V[1,:].^2*β[1]
    return nodes, weights
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
    If no `N` is provided, then `N = length(α)`.
"""
function gauss(N::Int64,α::Vector{Float64},β::Vector{Float64})
    @assert N > 0 "only positive N allowed"
    @assert length(α) == length(β) "inconsistent number of recurrence coefficients"
    N0 = length(α)
    @assert N0 >= N "not enough recurrence coefficients"
    golubwelsch(α[1:N],β[1:N])
end
gauss(α::Vector{Float64},β::Vector{Float64}) = gauss(length(α),α,β)
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
nodes and weights `(n+1)`-point Gauss-Radau
quadrature rule for the weight function having a prescribed
node `end0` (typically at one of the end points of the support
interval of w, or outside thereof).

!!! note
    If no `N` is specified, then `N = length(α)-1`.

!!! note
    Reference: OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS AND RELATED QUADRATURE RULES by Walter Gautschi
"""
function radau(N::Int64,α::Vector{Float64},β::Vector{Float64},end0::Float64)
    @assert N > 0 "only positive N allowed"
    @assert length(α) == length(β) > 0. "inconsistent number of recurrence coefficients"
    @assert length(α) >= N + 1 "not enough recurrence coefficients"
    p0 = 0.
    p1 = 1.
    for n in Base.OneTo(N)
      pm1 = p0
      p0 = p1;
      p1 = (end0 - α[n])*p0 - β[n]*pm1;
    end
    α[N+1] = end0 - β[N+1]*p0/p1
    gauss(N+1,α,β)
end
radau(α::Vector{Float64},β::Vector{Float64},end0::Float64) = radau(length(α)-1,α,β,end0)
radau(N::Int64,op::OrthoPoly,end0::Float64) = radau(N,op.α,op.β,end0)
radau(op::OrthoPoly,end0::Float64) = radau(op.α,op.β,end0::Float64)

"""
    radau_jacobi(N::Int64,a::Float64,b::Float64;endpoint::String="left")
    endpoint in ["left", "right"]
    radau_jacobi(N::Int64,a::Float64;endpoint::String="left") = radau_jacobi(N,a,a;endpoint=endpoint)
    radau_jacobi(N::Int64;endpoint::String="left") = radau_jacobi(N,0.;endpoint=endpoint)
Gauss-Radau quadrature rule for Jacobi weight function, which
generates the `(n+1)`-point Gauss-Radau
rule for the Jacobi weight function on `[-1,1]` with parameters
`a` and `b`.

!!! note
    REFERENCE: W. Gautschi, ``Gauss-Radau formulae for Jacobi and
    Laguerre weight functions'', Math. Comput. Simulation 54
    (2000), 403-412.
"""
function radau_jacobi(N::Int64,a::Float64,b::Float64;endpoint::String="left")
    @assert N>0 "only positive N allowed"
    endpoint = lowercase(endpoint)
    @assert endpoint in ["left", "right"] "$endpoint is no valid specification"
    α, β = rm_jacobi(N+1,a,b);
    α[N+1] = endpoint == "left" ? -1+2*N*(N+a)/((2*N+a+b)*(2*N+a+b+1)) : 1-2*N*(N+b)/((2*N+a+b)*(2*N+a+b+1))
    gauss(α,β)
end
radau_jacobi(N::Int64,a::Float64;endpoint::String="left") = radau_jacobi(N,a,a;endpoint=endpoint)
radau_jacobi(N::Int64;endpoint::String="left") = radau_jacobi(N,0.;endpoint=endpoint)

"""
    radau_laguerre(N::Int64,a::Float64)
    radau_laguerre(N::Int64) = radau_laguerre(N,0.)
Gauss-Radau quadrature rule for Laguerre weight function, which
generates the `(n+1)`-point Gauss-Radau
rule for the Laguerre weight function on ``[0,\\infty]`` with parameter `a`.

!!! note
    REFERENCE: W. Gautschi, ``Gauss-Radau formulae for Jacobi and
    Laguerre weight functions'', Math. Comput. Simulation 54
    (2000), 403-412.
"""
function radau_laguerre(N::Int64,a::Float64)
    @assert N>=0 "only positive N allowed"
    α, β = rm_laguerre(N+1,a);
    α[N+1] = N
    gauss(α,β)
end
radau_laguerre(N::Int64) = radau_laguerre(N,0.)

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
resp. to the right therof).

!!! note
    If no `N` is specified, then `N = length(α)-2`.

!!! note
    Reference: OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS AND RELATED QUADRATURE RULES by Walter Gautschi
"""
function lobatto(N::Int64,α::Vector{Float64},β::Vector{Float64},endl::Float64,endr::Float64)
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
        p1l=(endl-α[n])*p0l-β[n]*pm1l
        p1r=(endr-α[n])*p0r-β[n]*pm1r
    end
    det = p1l * p0r - p1r * p0l
    α[N+2] = (endl * p1l * p0r - endr * p1r * p0l) / det
    β[N+2] = (endr - endl) * p1l * p1r / det;
    gauss(N+2,α,β)
end
lobatto(α::Vector{Float64},β::Vector{Float64},endl::Float64,endr::Float64) = lobatto(length(α)-2,α,β,endl,endr)
lobatto(N::Int64,op::OrthoPoly,endl::Float64,endr::Float64) = lobatto(N,op.α,op.β,endl,endr)
lobatto(op::OrthoPoly,endl::Float64,endr::Float64) = lobatto(op.α,op.β,endl,endr)

"""
    lobatto_jacobi(N::Int64,a::Float64,b::Float64)
    lobatto_jacobi(N::Int64,a::Float64) = lobatto_jacobi(N,a,a)
    lobatto_jacobi(N::Int64) = lobatto_jacobi(N,0.)
Gauss-Lobatto quadrature rule for Jacobi weight function, which
generates the `(n+2)`-point Gauss-
Lobatto rule for the Jacobi weight function on ``[-1,1]`` with
parameters `a` and `b`.

!!! note
    REFERENCE: W. Gautschi,``High-order Gauss-Lobatto formulae'',
    Numer. Algorithms 25 (2000), 213-222.
"""
function lobatto_jacobi(N::Int64,a::Float64,b::Float64)
    @assert N > 0 "only positive N allowed"
    α,β = rm_jacobi(N+2,a,b);
    α[N+2] = (a-b) / (2*N+a+b+2)
    β[N+2] = 4 * (N + a + 1) * (N + b + 1) * (N + a + b + 1) / ((2 * N + a + b + 1) * (2*N + a + b + 2)^2)
    gauss(α,β)
end
lobatto_jacobi(N::Int64,a::Float64) = lobatto_jacobi(N,a,a)
lobatto_jacobi(N::Int64) = lobatto_jacobi(N,0.)

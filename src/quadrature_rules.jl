export  fejer,
        fejer2,
        clenshaw_curtis,
        quadgp,
        golubwelsch,
        quadpts_beta01,
        quadpts_gamma,
        quadpts_gaussian,
        quadpts_logistic,
        quadpts_uniform01,
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
    @assert N>=1
    θ = (2*collect(1:N).-1)*pi/(2*N)
    t = cos.(θ)
    w = [ 2/N*( 1 - 2*sum( cos(2*n*θ[ν])/(4*n^2-1) for n=1:floor(Int,N/2)) ) for ν=1:N]
    return t, w
end

"""
    fejer2(n::Int64)
Fejer's second quadrature rule according to [Waldvogel, J. Bit Numer Math (2006) 46: 195](https://doi.org/10.1007/s10543-006-0045-4).
"""
function fejer2(n::Int64)
    x = cos.(collect(0:n)/n*pi)
    N=collect(1:2:n-1);
    l=length(N);
    m=n-l;
    K=collect(0:m-1);
    v0=[2 ./ N ./(N.-2);  1/N[end]; zeros(m)];
    v2=-v0[1:end-1]-v0[end:-1:2];
    wf2=real(ifft(v2));
    return x, [ wf2; wf2[1] ]
end

"""
    clenshaw_curtis(n::Int64)
Clenshaw-Curtis quadrature according to [Waldvogel, J. Bit Numer Math (2006) 46: 195](https://doi.org/10.1007/s10543-006-0045-4).
"""
function clenshaw_curtis(n::Int64)
    x = cos.(collect(0:n)/n*pi)
    N=collect(1:2:n-1);
    l=length(N);
    m=n-l;
    K=collect(0:m-1);
    v0=[2 ./N./(N .- 2);  1/N[end]; zeros(m)];
    v2=-v0[1:end-1]-v0[end:-1:2];
    g0=-ones(n);
    g0[1+l]=g0[1+l]+n;
    g0[1+m]=g0[1+m]+n;
    g=g0/(n^2-1+mod(n,2));
    wcc=real(ifft(v2+g));
    return x,[ wcc; wcc[1] ]
end

"""
    quadgp(weight::Function,a::Float64,b::Float64,N::Int64=10;quadrature::Function=clenshaw_curti
general purpose quadrature based on Gautschi, "Orthogonal Polynomials: Computation and Approximation", Section 2.2.2, pp. 93-95
"""
function quadgp(weight::Function,lb::Float64,ub::Float64,N::Int64=10;quadrature::Function=clenshaw_curtis,bnd::Float64=Inf)
    @assert lb<ub "inconsistent interval bounds"
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

function JacobiGW( n::Int64, a::Float64, b::Float64 )
    # # Golub-Welsch for Gauss--Jacobi quadrature. This is used when max(a,b)>5.
    # ab = a + b;
    # ii = 2:n-1;
    # abi = 2*ii + ab;
    # aa = Float64[(b - a)/(2 + ab);
    #       (b^2 - a^2)./((abi - 2).*abi);
    #       (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))] ::Vector{Float64}
    # bb = Float64[2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ;
    #       2 .*sqrt.(ii.*(ii .+ a).*(ii .+ b).*(ii .+ ab)./(abi.^2 .- 1))./abi] ::Vector{Float64}
    # TT = SymTridiagonal(aa, bb)  # Jacobi matrix.
    # x, V = eig( TT )                       # Eigenvalue decomposition.
    # # Quadrature weights:
    # #w = V[1,:].^2 .*( 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(2+ab) );
    # w = V[1,:].^2*2^(ab+1)*beta(a+1, b+1);
    # C = 2^(a+b+1)*beta(a+1, b+1)/sum(w)
    # #@show C
    # #w = C.*w;
    # x, vec(w)
    # Golub-Welsh for Gauss--Jacobi quadrature. This is used when max(a,b)>5.
   ab = a + b;
   ii = 2:n-1;
   abi = 2*ii .+ ab;
   aa = Float64[(b - a)/(2 + ab);
         (b^2 - a^2)./((abi .- 2).*abi);
         (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))] ::Vector{Float64}
   bb = Float64[2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ;
         2 .*sqrt.(ii.*(ii .+ a).*(ii .+ b).*(ii .+ ab)./(abi.^2 .- 1))./abi] ::Vector{Float64}
   TT = SymTridiagonal(aa, bb)  # Jacobi matrix.
   x, V = eigen( TT )                       # Eigenvalue decomposition.
   # Quadrature weights:
   w = V[1,:].^2 .*( 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(2+ab) );
   # w .= w./sum(w);
x, vec(w)
end

function golubwelsch(α::Vector{Float64},β::Vector{Float64})
    J = SymTridiagonal( α, sqrt.(β[2:end]) )
    nodes, V = eigen(J)
    # nodes, V = e.values, e.vectors
    weights = V[1,:].^2*β[1]
    return nodes, weights
end
golubwelsch(op::OrthoPoly) = golubwelsch(op.α,op.β)

function quadpts_hermite(N::Int64)
    if N==0
        return Float64[], Float64[]
    elseif N==1
        return [0.0], [sqrt(pi)]
    elseif 1<N<=20
        a,b= rm_hermite(N)
        return golubwelsch(a,b)
    else
        return gausshermite(N)
    end
end



"""
    quadpts_beta01(α::Float64,β::Float64,Nquad::Int64)
get quadrature points for beta distribution on ``(0,1)``` using Gauss-Jacobi quadrature
"""
function quadpts_beta01(Nquad::Int64,α::Float64,β::Float64)
    nodes, weights = gaussjacobi(Nquad,β-1,α-1)
    return 0.5*(1 .+ nodes), weights/(beta(α,β)*2^(β+α-1))
end
quadpts_beta01(op::OrthoPoly,Nquad::Int) = quadpts_beta01(Nquad,op.pars[:shape_a],op.pars[:shape_b])

"""
    quadpts_gaussian(Nquad::Int)
get quadrature points for normal distribution on ``(-\\infty,\\infty)`` using Gauss-Hermite quadrature
"""
function quadpts_gaussian(N::Int)
    if N==0
        return Float64[], Float64[]
    elseif N==1
        return [0.0], [1.0]
    elseif 1<N<=20
        a,b= r_scale(1/sqrt(2pi),rm_hermite_prob(N)...)
        return golubwelsch(a,b)
    else
        nodes, weights = gausshermite(N)
        return sqrt(2)*nodes, weights/sqrt(pi)
    end
end
quadpts_gaussian(op::OrthoPoly,Nquad::Int64) = quadpts_gaussian(Nquad)

"""
    quadpts_uniform01(Nquad::Int)
get quadrature points for uniform distribution on ``[0,1]`` using Gauss-Legendre quadrature
"""
function quadpts_uniform01(Nquad::Int)
    nodes, weights = gausslegendre(Nquad)
    return 0.5*(1 .+ nodes), 0.5*weights
end
quadpts_uniform01(op::OrthoPoly,Nquad::Int) = quadpts_uniform01(Nquad)

"""
    quadpts_gamma(α::Float64,Nquad::Int)
get quadrature points for gamma distribution on ``(0,\\infty)`` using Gauss-Laguerre quadrature
"""
function quadpts_gamma(Nquad::Int64,α::Float64)
    nodes, weights = gausslaguerre(Nquad,α-1)
    return nodes, weights/gamma(α)
end
quadpts_gamma(op::OrthoPoly,Nquad::Int64) = quadpts_gamma(Nquad,op.pars[:shape])

function quadpts_logistic(N::Int64)
    # get recursion coefficients
    α, β = rm_logistic(N)
    # get quadrature rule from recursion coefficients
    nodes, weights = golubwelsch(α,β)
    return nodes, weights
end
quadpts_logistic(op::OrthoPoly,N::Int64) = quadpts_logistic(N)

"""
    quadpts_logistic(N::Int,c1::Real,c2::Real,c3::Real=1.)
    quadpts_logistic(N::Int64)
`N`-point quadrature rule for weight function
```math
    w(t) = c_3 c_2\\frac{\\exp(-c_2(t-c_1))}{(1+\\exp(-c_2(t-c_1)))^2}
```

The default value for c3 is one.
In that case ``w(t)`` is the probability density function of
```math
    Y = \\frac{1}{c_2 X + c_1}, \\quad c_2>0
```
where X has the standard logistic density
```math
    ρ(t) = \\frac{\\exp(-t)}{(1+\\exp(-t))^2}
```
The `N`-point quadrature rule for ``ρ(t)`` is computed by calling `quadpts_logistic(N::Int64)`.

"""
function quadpts_logistic(N::Int,c1::Real,c2::Real,c3::Real=1.)
    @assert abs(1/c2)>=1000*eps() "Value for c2 almost zero."
    nodes, weights = quadpts_logistic(N)
    return c1+1/c2*nodes, weights*c3
end

####################################################
####################################################
####################################################
# Quadrature rules based on recurrence coefficients

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
    @assert N>0 "only positive N allowed"
    @assert length(α)==length(β) "inconsistent number of recurrence coefficients"
    N0 = length(α)
    @assert N0>=N "not enough recurrence coefficients"
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
    @assert N>0 "only positive N allowed"
    @assert length(α)==length(β)>0 "inconsistent number of recurrence coefficients"
    N0 = length(α)
    @assert N0>=N+1 "not enough recurrence coefficients"
    ab0::Matrix{Float64}=[α β];
    p0::Float64=0.; p1::Float64=1.;
    for n=1:N
      pm1=p0; p0=p1;
      p1=(end0-ab0[n,1])*p0-ab0[n,2]*pm1;
    end
    ab0[N+1,1]=end0-ab0[N+1,2]*p0/p1;
    gauss(N+1,ab0[:,1],ab0[:,2])
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
    α,β=rm_jacobi(N+1,a,b);
    ab::Matrix{Float64} = [α β]
    if endpoint=="left"
      ab[N+1,1]=-1+2*N*(N+a)/((2*N+a+b)*(2*N+a+b+1));
    else
      ab[N+1,1]=1-2*N*(N+b)/((2*N+a+b)*(2*N+a+b+1));
    end
    gauss(ab[:,1],ab[:,2]);
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
    α,β=rm_laguerre(N+1,a);
    ab::Matrix{Float64} = [α β]
    ab[N+1,1]=N;
    gauss(ab[:,1],ab[:,2]);
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
    @assert N>0 "only positive N allowed"
    @assert length(α)==length(β)>0 "inconsistent number of recurrence coefficients"
    @assert endl<endr "inconsistent end points"
    N0 = length(α)
    @assert N0>=N+2 "not enough recurrence coefficients"
    ab0::Matrix{Float64}=[α β];
    p0l::Float64=0.; p0r::Float64=0.; p1l::Float64=1.; p1r::Float64=1.;
    for n=1:N+1
        pm1l=p0l; p0l=p1l; pm1r=p0r; p0r=p1r;
        p1l=(endl-ab0[n,1])*p0l-ab0[n,2]*pm1l;
        p1r=(endr-ab0[n,1])*p0r-ab0[n,2]*pm1r;
    end
    det::Float64=p1l*p0r-p1r*p0l;
    ab0[N+2,1]=(endl*p1l*p0r-endr*p1r*p0l)/det;
    ab0[N+2,2]=(endr-endl)*p1l*p1r/det;
    gauss(N+2,ab0[:,1],ab0[:,2]);
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
    @assert N>0 "only positive N allowed"
    # if nargin<2, a=0; end; if nargin<3, b=a; end
    α,β=rm_jacobi(N+2,a,b);
    ab::Matrix{Float64} = [α β]
    ab[N+2,1]=(a-b)/(2*N+a+b+2);
    ab[N+2,2]=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)*(2*N+a+b+2)^2);
    gauss(ab[:,1],ab[:,2]);
end
lobatto_jacobi(N::Int64,a::Float64) = lobatto_jacobi(N,a,a)
lobatto_jacobi(N::Int64) = lobatto_jacobi(N,0.)

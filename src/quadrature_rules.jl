
"""
    fejer(N::Int64)
Fejer's first quadrature rule.
"""
function fejer(N::Int64)
    @assert N>=1
    θ = (2*collect(1:N)-1)*pi/(2*N)
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
    v0=[2 ./N./(N-2);  1/N[end]; zeros(m)];
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
function quadgp(weight::Function,lb::Float64,ub::Float64,N::Int64=10;quadrature::Function=clenshaw_curtis)
    @assert lb<ub "inconsistent interval bounds"
    bnd=1e4
    t_fej,w_fej = quadrature(N)
    t,w,d = zeros(N),zeros(N),zeros(N)
    # transformation
    if (-bnd<lb<bnd) && (-bnd<ub<bnd)
        t = 0.5*((ub-lb)*t_fej .+ (ub+lb))
        d = 0.5*(ub-lb)
    elseif -bnd>=lb && (-bnd<ub<bnd)
        t = ub-(1 .- t_fej)./(1 .+ t_fej)
        d = 2. /(1 .+ t_fej).^2
    elseif (-bnd<lb<bnd) && ub>=bnd
        t = lb+(1 .+ t_fej)./(1 .- t_fej)
        d = 2. /(1 .- t_fej).^2
    elseif -bnd>=lb && ub>=bnd
        t = t_fej./(1 .- t_fej.^2)
        d = (1 .+ t_fej.^2)./( 1 .- t_fej.^2).^2
    end
    w = w_fej.*weight.(t).*d
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

"""
    quadpts_logistic(N::Int64)
get quadrature points for logistic weight function on ``(-\\infty,\\infty)``
"""
function quadpts_logistic(N::Int64)
    # get recursion coefficients
    α, β = rm_logistic(N)
    # get quadrature rule from recursion coefficients
    nodes, weights = golubwelsch(α,β)
    return nodes, weights
end
quadpts_logistic(op::OrthoPoly,N::Int64) = quadpts_logistic(N)

"""
    Quadrature rule for weight function
        w(t) = c3*c2*exp(-c2*(t-c1))/(exp(-c2*(t-c1)))^2

    The default value for c3 is one, c3=1, and
        w(t) = c2*exp(-c2*(t-c1))/(1+exp(-c2*(t-c1)))^2.
    In that case w(t) is the probability density function of
        Y = 1/c2*X + c1, with c2>0
    where X has the standard logistic density
        ρ(t) = exp(-x)/(1+exp(-x))^2

    A value c3!=1 is needed, for example, when constructing sums of logistic densities.
"""
function quadpts_logistic(N::Int,c1::Real,c2::Real,c3::Real=1.)
    @assert abs(1/c2)>=1000*eps() "Value for c2 almost zero."
    nodes, weights = quadpts_logistic(N)
    return c1+1/c2*nodes, weights*c3
end

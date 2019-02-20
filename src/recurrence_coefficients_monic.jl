export  r_scale,
        rm_logistic,
        rm_laguerre,
        rm_logisticsum,
        rm_hermite,
        rm_hermite_prob,
        rm_laguerre,
        rm_jacobi,
        rm_jacobi01,
        rm_hahn,
        rm_meixner_pollaczek,
        rm_legendre,
        rm_legendre01,
        rm_chebyshev1,
        rm_compute
"""
    r_scale(c::Float64,β::Vector{Float64},α::Vector{Float64})
Given the recursion coefficients `(α,β)` for a system of orthogonal polynomials that are orthogonal with respect to some positive weight ``m(t)``,
this function returns the recursion coefficients `(α_,β_)` for the scaled measure ``c m(t)`` for some positive ``c``.
"""
function r_scale(c::Float64,a::Vector{Float64},b::Vector{Float64})
    @assert c>0 "Measure can only be scaled by positive number (provided c=$c)"
    return a, [c*b[1]; b[2:end]]
end

"""
    rm_compute(weight::Function,lb::Float64,ub::Float64;Npoly::Int64=4,Nquad::Int64=10,quadrature::Function=clenshaw_curtis)
Given a positive `weight` function with domain `(lb,ub)`, i.e. a function ``w: [lb, ub ] \\rightarrow \\mathbb{R}_{\\geq 0}``,
this function creates `Npoly` recursion coefficients `(α,β)`.

The keyword `quadrature` specifies what quadrature rule is being used.
"""
function rm_compute(weight::Function,lb::Float64,ub::Float64;Npoly::Int64=4,Nquad::Int64=10,quadrature::Function=clenshaw_curtis,discretization::Function=stieltjes)
    # @assert discretization in [ stieltjes, lanczos ]
    @assert Npoly<=Nquad
    Npoly==0 ? (return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)) : ()
    n,w=quadgp(weight,lb,ub,Nquad;quadrature=quadrature)
    a,b = discretization(Npoly,n,w)
    # display((b,maximum(b)))
    # the following assertions are taken from Theorem 1.28 from
    # Gautschi, W. "Orthogonal Polynomials: Computation and Approximation"
    @assert lb<minimum(a) && maximum(a)< ub "Not all recurrence coefficients α ($(minimum(a)),$(maximum(b))) are elements of the support s=($lb,$ub)."
    if length(b)>=2
        @assert 0<minimum(b) && maximum(b[2:end])<=max(lb^2,ub^2) "Not all recurrence coefficients β are elements of the interval (0, max(lb^2,ub^2)) = (0,$(max(lb^2,ub^2))), namely ($(minimum(b)), $(maximum(b)))."
    end
    return a,b
end

rm_compute(m::Measure;Npoly::Int64=4,Nquad::Int64=10,quadrature::Function=clenshaw_curtis,discretization::Function=stieltjes) = rm_compute(m.w,m.dom[1],m.dom[2];Npoly=Npoly,Nquad=Nquad,quadrature=quadrature,discretization=discretization)

##
function rm_logisticsum(n::Int,p1::Vector{Float64},p2::Vector{Float64};Mmax::Int=100,eps0::Real=1e-9)
    M0 = n
    Mcap = 0
    Mi = M0
    α, β = zeros(n), zeros(n)
    for i=1:Mmax
        nai, wai = quadrature_logistic(Mi,p1[1],p1[2],0.5) # 0.5 is needed to make the weight function a density, i.e. integral is 1.
        nbi, wbi = quadrature_logistic(Mi,p2[1],p2[2],0.5) # 0.5 is needed to make the weight function a density, i.e. integral is 1.
        ni, wi   = [nai; nbi], [wai; wbi]
        α_, β_ = stieltjes(n,ni,wi)
        err = maximum( abs(β_[k]-β[k])/β_[k] for k=1:n )
        # display("err = $err")
        if err<=eps0
            Mcap = i-1
            return α_, β_
        else
            α, β, = α_, β_
            i==1 ? Mi = M0+1 : Mi = Int(2^(floor((i-1)/5))*n)
        end
    end
    warn("Algorithm did not terminate after $Mmax iterations.")
end

"""
    rm_logistic(N::Int)
Creates `N` recurrence coefficients for monic polynomials that are orthogonal
on ``(-\\infty,\\infty)`` relative to ``w(t) = \\frac{\\mathrm{e}^{-t}}{(1 - \\mathrm{e}^{-t})^2}``
"""
function rm_logistic(N::Int)
    @assert N >= 0 "parameter(s) out of range."
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    zeros(N), map(k->k^4*pi^2/(4*k^2-1),Base.OneTo(N-1))
end

"""
    rm_hermite(N::Int,mu::Float64)
    rm_hermite(N::Int)
Creates `N` recurrence coefficients for monic generalized Hermite polynomials
that are orthogonal on ``(-\\infty,\\infty)`` relative to ``w(t) = |t|^{2 \\mu} \\mathrm{e}^{-t^2}``

The call `rm_hermite(N)` is the same as `rm_hermite(N,0)`.
"""
function rm_hermite(N::Int,mu::Float64)
    @assert N >= 0 && mu > -0.5 "parameter(s) out of range."
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    m0 = mu != 0. ? gamma(mu + 0.5) : sqrt(π)
    N == 1 && return [0.], [m0]
    return zeros(N), pushfirst!(map(x->isodd(x) ? 0.5*x+mu : 0.5*x,1:N-1), m0)
end
rm_hermite(N::Int) = rm_hermite(N,0.)

"""
    rm_hermite_prob(N::Int)
Creates `N` recurrence coefficients for monic probabilists' Hermite polynomials
that are orthogonal on ``(-\\infty,\\infty)`` relative to ``w(t) = \\mathrm{e}^{-0.5t^2}``
"""
function rm_hermite_prob(N::Int)
    @assert N >= 0 "parameter(s) out of range."
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    # return zeros(N), [sqrt(2*pi); collect(1.:N-1) ]
    return zeros(N), pushfirst!( collect(1.:N-1), sqrt(2*pi))
end

"""
    rm_laguerre(N::Int,a::Float64)
    rm_laguerre(N::Int)
Creates `N` recurrence coefficients for monic generalized Laguerre polynomials
that are orthogonal on ``(0,\\infty)`` relative to ``w(t) = t^a \\mathrm{e}^{-t}``.

The call `rm_laguerre(N)` is the same as `rm_laguerre(N,0)`.
"""
function rm_laguerre(N::Int,a::Float64)
    @assert N >= 0 && a > -1. "parameter(s) out of range"
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    N==1 && return [a+1.], [gamma(a+1)]
    n = 1.:N-1
    return pushfirst!(map(x->2x+a+1.,n)::Vector{Float64},a+1.), pushfirst!(map(x->x^2+a*x,n)::Vector{Float64},gamma(a+1))
end
function rm_laguerre(N::Int)
    rm_laguerre(N,0.)
end

"""
    rm_jacobi(N::Int,a::Float64,b::Float64)
    rm_jacobi(N::Int,a::Float64)
    rm_jacobi(N::Int)
Creates `N` recurrence coefficients for monic Jacobi polynomials
that are orthogonal on ``(-1,1)`` relative to ``w(t) = (1-t)^a (1+t)^b``.

The call `rm_jacobi(N,a)` is the same as `rm_jacobi(N,a,a)` and `rm_jacobi(N)` the same as
`rm_jacobi(N,0,0)`.
"""
function rm_jacobi(N::Int,a::Float64,b::Float64)
    @assert N >= 0 && a >- 1. && b >- 1. "parameter(s) out of range"
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    nu = (b - a) / ( a+ b + 2.);
    mu::Float64 = a + b + 2. <= 128. ? 2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2)) : exp((a+b+1)*log(2)+((log(gamma(a+1))+log(gamma(b+1)))-log(gamma(a+b+2))))
    N == 1 && return [nu], [mu]
    n = 1:N-1
    nab = map(x->2x+a+b,n)
    A = pushfirst!(map(x->(b^2-a^2)/(x*(x+2)),nab),nu)
    B1 = map(x->4*(x+a)*(x+b)*x*(x+a+b),n[2:N-1])
    B2 = map(x->(x^2.)*(x+1.)*(x-1.),nab[2:N-1])
    B3 = 4. * (a + 1) * (b + 1) / ((a + b + 2)^2 * (a + b + 3));
    A, pushfirst!(pushfirst!(B1./B2,B3),mu)::Vector{Float64}
end

rm_jacobi(N::Int,a::Float64) = rm_jacobi(N,a,a)
rm_jacobi(N::Int) = rm_jacobi(N,0.,0.)
"""
    rm_jacobi01(N::Int,a::Float64,b::Float64)
    rm_jacobi01(N::Int,a::Float64)
    rm_jacobi01(N::Int)

Creates `N` recurrence coefficients for monic Jacobi polynomials
that are orthogonal on ``(0,1)`` relative to ``w(t) = (1-t)^a t^b``.

The call `rm_jacobi01(N,a)` is the same as `rm_jacobi01(N,a,a)` and `rm_jacobi01(N)` the same as
`rm_jacobi01(N,0,0)`.
"""
function rm_jacobi01(N::Int,a::Float64,b::Float64)
    @assert N>=0 && a>-1. && b>-1. "parameter(s) out of range"
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    c, d = rm_jacobi(N,a,b)
    map(x->(1+x)/2,c), pushfirst!(0.25*d[2:N], d[1]/2^(a+b+1.))
end
rm_jacobi01(N::Int,a::Float64) = rm_jacobi01(N,a,a)
rm_jacobi01(N::Int) = rm_jacobi01(N,0.,0.)

"""
    rm_legendre(N::Int)

Creates `N` recurrence coefficients for monic Legendre polynomials
that are orthogonal on ``(-1,1)`` relative to ``w(t) = 1``.
"""
rm_legendre(N::Int) = rm_jacobi(N)

"""
    rm_legendre01(N::Int)

Creates `N` recurrence coefficients for monic Legendre polynomials
that are orthogonal on ``(0,1)`` relative to ``w(t) = 1``.
"""
rm_legendre01(N::Int) = rm_jacobi01(N)


function rm_chebyshev1(N::Int64)
    @assert N>=0 "N has to be non-negative"
    α = zeros(Float64,N)
    if N == 1
        return α, [pi]
    elseif N == 2
        return α, [pi; 0.5]
    else
        return α, pushfirst!(pushfirst!(0.25*ones(N-2),0.5),pi)
    end
end
#     rm_hahn(N::Int,a::Float64,b::Float64)
#     rm_hahn(N::Int,a::Float64)
#     rm_hahn(N::Int)
#
#  Creates `N` recurrence coefficients for monic Hahn polynomials.
#
#     ab=rm_hahn(N,a,b) generates the N+1 recurrence coefficients
#     for the monic Hahn polynomials with parameters a and b. These
#     are orthogonal on the discrete set of N+1 points 0,1,2,...,N
#     with weights (a+k choose k)(b+N-k choose N-k), k=0,1,2,...,N.
#     The N+1 alpha-coefficients are stored in the first column,
#     the N+1 beta-coefficients in the second column, of the
#     (N+1)x2 array ab. The call ab=rm_hahn(N,a) is the same as ab=
#     rm_hahn(N,a,a), and ab=HAHN(N) the same as ab=rm_hahn(N,0,0),
#     which produces the recurrence coefficients of the discrete
#     Chebyshev polynomials for the points 0,1,2,...,N.
# """
# function rm_hahn(N::Int,a::Float64,b::Float64)
#     N-=1
#     @assert N>=0 && a>=-1. && b>=-1. "parameter(s) out of range"
#     N==0 ? (return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)) : ()
#     ab = zeros(N+1,2)
#     n=1:N; ab[1,2]=prod(1 .+ (a+b+1.)./n);
#     if a+b == 0.
#       n=(1:N+1)';
#       ab[n,1]=((2*n .+ (a+b-1.))*N+(b-a)*n .+ a)./(2*(2*n .- 1));
#       n=(1:N)';
#       ab[n.+1,2]=.25*((N+1)^2)*(1 .+ a./n).*(1 .+ b./n).*(1 .-(n./(N+1)).^2)./(4 .- (1 ./n).^2);
#     elseif a+b+1. == 0.
#       n=(1:N+1)';
#       ab[n,1]=((2*(n .- 1).^2 .+ b)*N .+ (2*b+1)*(n .- 1).^2)./(4*(n .- 1).^2 .- 1);
#       n=(1:N)';
#       ab[n.+1,2]=(.25*((N+1)^2)).*(1 .+ a*(1 ./n)).*(1 .+ b./n).*(1 .- n./(N+1)).*(1 .+ (n .- 1)./(N+1))./(4 .- (1 ./n).^2);
#     else
#       n=(1:N+1)';
#       ab[n,1]=((n .+ (a+b)).*(n .+ a).*((N+1) .- n)./(2*n .+ (a+b))+(n .- 1).*(n .+ (b-1)).*(N .+ n .+ (a+b))./(2*n .+ (a+b-2)))./(2*n .+ (a+b-1));
#       n=(1:N)';
#       ab[n.+1,2]=((N+1)^2)*(1 .+ a./n).*(1 .+ b./n).*(1 .+ (a+b)./n).*(1 .- n./(N+1)).*(1 .+ (n .+ (a+b))./(N .+ 1))./(((2 .+ (a+b)./n).^2).*((2 .+ (a+b)./n).^2 .- (1 ./n).^2));
#     end
#     return ab[:,1], ab[:,2]
# end
# rm_hahn(N::Int,a::Float64) = rm_hahn(N,a,a)
# rm_hahn(N::Int) = rm_hahn(N,0.,0.)

"""
    rm_meixner_pollaczek(N::Int,lambda::Float64,phi::Float64)
    rm_meixner_pollaczek(N::Int,lambda::Float64)

 Creates `N` recurrence coefficients for monic
 Meixner-Pollaczek polynomials with parameters λ and ϕ. These are orthogonal on
 ``[-\\infty,\\infty]`` relative to the weight function ``w(t)=(2 \\pi)^{-1} \\exp{(2 \\phi-\\pi)t} |\\Gamma(\\lambda+ i t)|^2``.

 The call `rm_meixner_pollaczek(n,lambda)` is the same as `rm_meixner_pollaczek(n,lambda,pi/2)`.
"""
function rm_meixner_pollaczek(N::Int,lambda::Float64,phi::Float64)
    @assert N>=0 && lambda>0. && phi>0. "parameter(s) out of range"
    N == 0 && return Array{Float64,1}(undef,0), Array{Float64,1}(undef,0)
    n=1:N;
    sinphi=sin(phi); lam2=2*lambda;
    ab=zeros(Float64,N,2)
    if sinphi==1
      ab[:,1]=zeros(N);
    else
        ab[:,1]= -(n .+ (lambda-1))/tan(phi);
    end
    ab[:,2]= (n .- 1).*(n .+ (lam2-2))/(4*sinphi^2);
    ab[1,2]=gamma(lam2)/(2*sinphi)^lam2;
    return ab[:,1], ab[:,2]
end
rm_meixner_pollaczek(N::Int,lambda::Float64) = rm_meixner_pollaczek(N,lambda,pi/2)

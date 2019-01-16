export stieltjes, lanczos

function removeZeroWeights(n::Vector{Float64},w::Vector{Float64})
    nw = [n w]
    inds = findall(w->abs(w)>=1e3*eps(), nw[:,2])
    nw = nw[inds,:]
    # nw = sortrows(nw, by=x->x[1])
    nw = sortslices(nw, dims=1, by=x->x[1])
    nw[:,1], nw[:,2]
end

"""
    stieltjes(N::Int64,nodes::Vector{Float64},weights::Vector{Float64})
Description based on W. Gautschi *OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS
AND RELATED QUADRATURE RULES*
Given the discrete inner product (with nodes and weights) the function
generates the first`N` recurrence coefficients of the
corresponding discrete orthogonal polynomials.
"""
function stieltjes(N::Int64,nodes_::Vector{Float64},weights_::Vector{Float64})
    tiny = 10*floatmin()
    huge = 0.1*floatmax()
    α, β = zeros(Float64,N), zeros(Float64,N)
    nodes, weights = removeZeroWeights(nodes_,weights_)
    Ncap = length(nodes)
    @assert N>0 && N<=Ncap "N is out of range."
    s0 = sum(weights)
    α[1] = dot(nodes,weights)/s0
    β[1] = s0

    if N==1  return (α, β) end
    p1, p2 = zeros(Ncap), ones(Ncap)
    for k=1:N-1
        p0 = p1
        p1 = p2
        p2 = (nodes .- α[k]).*p1 - β[k]*p0
        s1 = dot(weights,p2.^2)
        s2 = dot(nodes, weights.*p2.^2)
        # s1 = sum( weights[ν]*p2[ν]^2 for ν=1:Ncap )
        # s2 = sum( weights[ν]*nodes[ν]*p2[ν]^2 for ν=1:Ncap )
        @assert abs(s1)>=tiny "Underflow in stieltjes() for k=$k"
        @assert maximum(abs.(p2))<=huge && abs(s2)<=huge "Overflow in stieltjes for k=$k"
        α[k+1] = s2/s1
        β[k+1] = s1/s0
        s0 = s1
    end
    return α, β
end

"""
  lanczos(N::Int64,nodes_::Vector{Float64},weights_::Vector{Float64})
Description based on W. Gautschi *OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS
AND RELATED QUADRATURE RULES*

Given the discrete inner product (with nodes and weights) the function
generates the first `N` recurrence coefficients of the
corresponding discrete orthogonal polynomials.

The script is adapted from the routine RKPW in
W.B. Gragg and W.J. Harrod, ``The numerically stable
reconstruction of Jacobi matrices from spectral data'',
Numer. Math. 44 (1984), 317-335.
"""
function lanczos(N::Int64,nodes_::Vector{Float64},weights_::Vector{Float64})
    @assert length(nodes_)==length(weights_)>0 "inconsistent number of nodes and weights"
    nodes, weights = removeZeroWeights(nodes_,weights_)
    Ncap = length(nodes)
    N<=0 || N>Ncap ? error("$N out of range") : ()
    p0 = copy(nodes)
    p1 = zeros(Float64,Ncap);
    p1[1]= weights[1];
    for n=1:Ncap-1
      pn = weights[n+1]
      gam=1.; sig=0.; t=0.;
      xlam = nodes[n+1]
      for k=1:n+1
        rho=p1[k]+pn; tmp=gam*rho; tsig=sig;
        if rho<=0.
          gam=1.; sig=0.;
        else
          gam=p1[k]/rho;
          sig=pn/rho;
        end
        tk=sig*(p0[k]-xlam)-gam*t;
        p0[k]=p0[k]-(tk-t); t=tk;
        if sig<=0
          pn=tsig*p1[k];
        else
          pn=(t^2)/sig;
        end
        tsig=sig; p1[k]=tmp;
      end
    end
    p0[1:N], p1[1:N]
end


function mcdis_twologistics(n::Int,p1::Vector{Float64},p2::Vector{Float64};Mmax::Int=100,eps0::Real=1e-9)
    M0 = n
    Mcap = 0
    Mi = M0
    α, β = zeros(n), zeros(n)
    for i=1:Mmax
        nai, wai = quadrature_logistic(Mi,p1[1],p1[2])
        nbi, wbi = quadrature_logistic(Mi,p2[1],p2[2])
        wai, wbi = 0.5*wai, 0.5*wbi
        ni, wi   = [nai; nbi], [wai; wbi]
        α_, β_ = stieltjes(n,ni,wi)
        err = maximum( abs(β_[k]-β[k])/β_[k] for k=1:n )
        display("err = $err")
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

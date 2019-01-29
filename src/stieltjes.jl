export stieltjes, lanczos, mcdiscretization

function removeZeroWeights(n::Vector{Float64},w::Vector{Float64})
    nw = [n w]
    inds = findall(w->abs(w)>=eps(), nw[:,2])
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
    # nodes, weights = removeZeroWeights(nodes_,weights_)
    nodes, weights = nodes_, weights_
    Ncap = length(nodes)
    @assert N>0 && N<=Ncap "N is out of range."
    s0::Float64 = sum(weights)
    α[1] = dot(nodes,weights)/s0
    β[1] = s0

    if N==1  return (α, β) end
    p0, p1, p2 = zeros(Float64,Ncap), zeros(Float64,Ncap), ones(Float64,Ncap)
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
    # nodes, weights = removeZeroWeights(nodes_,weights_)
    nodes, weights = nodes_, weights_
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

"""
    mcdiscretization()
This routine performs a sequence of discretizations of the
given weight function (or measure), each discretization being
followed by an application of the Stieltjes, or Lanczos,
procedure to produce approximations to the desired recurrence
coefficients. The fineness of the discretization is
characterized by a discretization parameter `N`. The support of
the continuous part of the weight function is decomposed into
a given number `mc` of subintervals (some or all of which may
be identical). The routine then applies to each subinterval
an `N`-point quadrature rule to discretize the weight function
on that subinterval. The discrete part of the weight function
(if there is any) is added on to the discretized continuous
weight function. The sequence of discretizations, if chosen
judiciously, leads to convergence of the recurrence
coefficients for the discretized measures to those of the
given measure. If convergence to within a prescribed accuracy
`eps0` occurs before `N` reaches its maximum allowed value `Nmax`,
then the value of `N` that yields convergence is output as
`Ncap`, and so is the number of iterations, `kount`. If there is
no convergence, the routine displays the message "Ncap
exceeds Nmax in mcdis" prior to exiting.

The choice between the Stieltjes and the Lanczos procedure is
made by setting the parameter `discretization`.

The details of the discretization are to be specified prior
to calling the procedure. They are embodied in the following
global parameters:

mc     = the number of component intervals
mp     = the number of points in the discrete part of the
         measure (mp=0 if there is none)
iq     = a parameter to be set equal to 1, if the user
         provides his or her own quadrature routine, and
         different from 1 otherwise
δ = a parameter whose default value is 1, but is
         preferably set equal to 2, if iq=1 and the user
         provides Gauss-type quadrature routines

The component intervals have to be specified (in the order
left to right) by a global mcx2 array AB=[[a1 b1];[a2 b2];
...;[amc bmc]],  where for infinite extreme intervals a1=-Inf
resp. bmc=Inf. The discrete spectrum (if mp>0) is similarly
specified by a global mpx2 array DM=[[x1 y1];[x2 y2];...;
;[xmp ymp]] containing the abscissae and jumps.

If the user provides his or her own quadrature routine
"quadown", the routine mcdis must be called with the input
parameter "quad" replaced by "@quadown", otherwise with
"quad" replaced by "@quadgp", a general-purpose routine
provided in the package. The quadrature routine must have
the form

                         function xw=quad(N,i)

where N is the number of nodes and i identifies the interval
to which the routine is to be applied.

The routine mcdis also applies to measures given originally
in multi-component form.

"""
function mcdiscretization(N::Int64,
                        quads::Vector{},
                        discretemeasure::Matrix{Float64}=zeros(0,2);
                        discretization::Function=stieltjes,
                        Nmax::Integer=300,
                        ε::Float64=1e-8,
                        gaussquad::Bool=false)
    @assert Nmax>0 && Nmax>N "invalid choice of Nmax=$Nmax."
    @assert ε>0 "invalid choice of ε=$ε"
    @assert discretization in [stieltjes, lanczos] "unknown discretization $discretization"
    @assert length(quads)>0 "no quadrature rule specified"
    # @assert quadrature in [clenshaw_curtis, fejer, fejer2] "unknown quadrature $quadrature"
    # @assert size(AB)==(mc,2) "dimensions of continuous intervals are off"
    # @assert size(DM)==(mp,2) "dimensions of discrete intervals are off"
    # f='Ncap exceeds Nmax in mcdis with irout=%2.0f\n';
    δ::Int64=1
    gaussquad ? δ=2 : ()
    mc, mp = size(quads,1), size(discretemeasure,1)
    Δ::Int64=1; kount::Int64=-1;
    α, β, b = zeros(N), zeros(N), ones(N)
    Mi=floor(Int64,(2*N-1)/δ);
    while any(abs.(β-b).>ε*abs.(β))
        b=copy(β);
        @show kount=kount+1;
        kount>1 ? Δ=2^(floor(Int64,kount/5))*N : ()
        Mi+=Δ;
        Mi>Nmax ? error("Mi=$Mi exceeds Nmax=$Nmax") : ()
        Ntot=mc*Mi;
        xx, ww = zeros(Float64,Ntot), zeros(Float64,Ntot)
        for i=1:mc
            nn=(i-1)*Mi;
            x, w = quads[i](Mi)
            size(w)
            xx[nn+1:nn+Mi], ww[nn+1:nn+Mi] = x, w
        end
        if mp>0
            xx[Ntot+1:Ntot+mp], ww[Ntot+1:Ntot+mp] = discretemeasure[:,1], discretemeasure[:,2]
        end
        α,β = discretization(N,xx,ww)
        @show sum(ww)
        end
        return α, β
end

# create recursion coefficient of logistic density numerically
# and compare against analytic solution

using PolyChaos, Test
import LinearAlgebra: norm

function myquad(N::Int64,pos::Bool)
    a,b = rm_laguerre(N+1)
    n,w = golubwelsch(a,b)
    w = w ./ ((1 .+ exp.(-n)).^2)
    pos ? (n,w) : (-n,w)
end

ε, tol = 1e-10, 1e-6
N = 1:3:15
quads = [ n->myquad(n,true); n->myquad(n,false) ]

@testset "Logistic density" begin
    for n in N, disc in [stieltjes, lanczos]
        a,b = mcdiscretization(n,quads;discretization=disc,Nmax=300,gaussquad=true,ε=ε)
        aa,bb = rm_logistic(n)
        @test isapprox(norm(b-bb,Inf),0.;atol=tol)
    end
end


function hrhermite(t)
    exp(-t^2)
end

AB=[[0. 3.];[3. 6.];[6. 9.];[9. Inf]];
quads = [ n->quadgp(hrhermite,AB[i,1],AB[i,2],n;quadrature=fejer) for i=1:size(AB,1)]
α, β = mcdiscretization(40,quads,Nmax=400,gaussquad=false,discretization=lanczos,ε=1e-7)

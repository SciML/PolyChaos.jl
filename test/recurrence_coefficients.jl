# Test whether the numerically computed recursion coefficients match the analytic expressions

using PolyChaos, Test, LinearAlgebra
tol=1e-4
N = collect(0:10)
# hermite
mu = 0.:0.2:1

@testset "generalized Hermite" begin
    for N_ in N, mu_ in mu
        a_ana,b_ana = rm_hermite(N_,mu_)
        w(t) = build_w_hermite(t,mu_)
        a_num,b_num = rm_compute(w,-Inf,Inf;Npoly=N_,Nquad=max(3000,10*N_))
        @test isapprox(norm(a_ana-a_num,Inf),0.; atol=tol)
        @test isapprox(norm(b_ana-b_num,Inf),0.; atol=tol)
    end
end

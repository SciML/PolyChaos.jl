# Test whether the numerically computed recursion coefficients match the analytic expressions
using PolyChaos, Test
import LinearAlgebra: norm
tol=1e-4
N = collect(1:6)
# hermite
mu = 0.:0.2:1

test = []

newtest = @testset "generalized Hermite" begin
    for N_ in N, mu_ in mu, disc in [stieltjes]
        # display("N = $N_, μ = $mu_, $disc")
        a_ana,b_ana = rm_hermite(N_,mu_)
        w = build_w_genhermite(mu_)
        a_num,b_num = rm_compute(w,-Inf,Inf;Npoly=N_,Nquad=max(10000,10*N_),discretization=disc)
        @test isapprox(norm(a_ana-a_num,Inf),0.; atol=tol)
        @test isapprox(norm(b_ana-b_num,Inf),0.; atol=tol)
    end
end
push!(test,newtest)

a = 0.:1.:4
b = 0.:1.:4
# if shape parameters are negative, the weight till go to infinity at its bounds!
newtest = @testset "Jacobi" begin
    for N_ in N, a_ in a, b_ in b, disc in [stieltjes]
        # display("N = $N_, a = $a_, b = $b_, disc=$disc")
        a_ana,b_ana = rm_jacobi(N_,a_,b_)
        w = build_w_jacobi(a_,b_)
        a_num,b_num = rm_compute(w,-1.,1.;Npoly=N_,Nquad=max(1000,10*N_),discretization=disc)
        @test isapprox(norm(a_ana-a_num,Inf),0.; atol=tol)
        @test isapprox(norm(b_ana-b_num,Inf),0.; atol=tol)
    end
end

push!(test,newtest)

newtest = @testset "Jacobi01" begin
    for N_ in N, a_ in a, b_ in b, disc in [stieltjes]
        # display("N = $N_, a = $a_, b = $b_, disc=$disc")
        a_ana,b_ana = rm_jacobi01(N_,a_,b_)
        w = build_w_jacobi01(a_,b_)
        a_num,b_num = rm_compute(w,0.,1.;Npoly=N_,Nquad=max(1000,10*N_),discretization=disc)
        # @show a_ana-a_num
        @test isapprox(norm(a_ana-a_num,Inf),0.; atol=tol)
        @test isapprox(norm(b_ana-b_num,Inf),0.; atol=tol)
    end
end

push!(test,newtest)

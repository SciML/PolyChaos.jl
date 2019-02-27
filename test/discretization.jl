# testing discretization procedures
# example is based on Walter Gautschi "Orthogonal Polynomials: Computation and Approximation"
# examples 2.36

using PolyChaos, Test
import LinearAlgebra: norm

nodes = [40, 80, 160, 320]
tol = 1e-14

@time @testset "Stieltjes procedure" begin
    for n in nodes
        myfile = open("fej$n.txt")
        βref = parse.(Float64,readlines(myfile))
        βcom = stieltjes(n,fejer(n)...)[2]
        @test isapprox(norm(βref-βcom,Inf),0.;atol=tol)
        close(myfile)
    end
end

@time @testset "Lanczos procedure" begin
    for n in nodes
        myfile = open("fej$n.txt")
        βref = parse.(Float64,readlines(myfile))
        βcom = lanczos(n,fejer(n)...)[2]
        @test isapprox(norm(βref-βcom,Inf),0.;atol=tol)
        close(myfile)
    end
end

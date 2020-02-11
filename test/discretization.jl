# testing discretization procedures
# example is based on Walter Gautschi "Orthogonal Polynomials: Computation and Approximation"
# examples 2.36

using PolyChaos, Test, StaticArrays
import LinearAlgebra: norm

nodes = [40, 80, 160, 320]
tol = 1e-14

@testset "Stieltjes procedure" begin
    for n in nodes
        myfile = open("fej$n.txt")
        βref = parse.(Float64,readlines(myfile))
        βcom = stieltjes(n,fejer(n)...)[2]
        nodes, weights = fejer(n)
        @test stieltjes(n,SVector(nodes...),SVector(weights...)) == stieltjes(n,nodes,weights) == stieltjes(n,nodes,SVector(weights...)) == stieltjes(n,SVector(nodes...),weights)
        @test isapprox(norm(βref-βcom,Inf),0.;atol=tol)
        close(myfile)
    end
end

@testset "Lanczos procedure" begin
    for n in nodes
        myfile = open("fej$n.txt")
        βref = parse.(Float64,readlines(myfile))
        nodes, weights = fejer(n)
        βcom = lanczos(n,nodes,weights)[2]
        @test lanczos(n,SVector(nodes...),SVector(weights...)) == lanczos(n,nodes,weights) == lanczos(n,nodes,SVector(weights...)) == lanczos(n,SVector(nodes...),weights)
        @test isapprox(norm(βref-βcom,Inf),0.;atol=tol)
        close(myfile)
    end
end

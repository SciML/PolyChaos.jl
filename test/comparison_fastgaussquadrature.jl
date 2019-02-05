using PolyChaos, FastGaussQuadrature, Test, LinearAlgebra

mytol = 1e-10
N = 5:10:100

names = ["hermite", "legendre", "laguerre"]
funs = [gausshermite, gausslegendre, gausslaguerre]

@testset "Compare $(names...) against FastGaussQuadrature" begin
    for i=1:length(names), n in N
        opq = OrthoPolyQ(names[i],n-1)
        mynw = nw(opq)
        n, w = funs[i](n)
        @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
        @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
    end
end

a = -0.5:0.1:2.0
b = copy(a)

@testset "Compare Jacobi against FastGaussQuadrature" begin
    for aa in a, bb in b, n in N
        opq = OrthoPolyQ("jacobi",n-1,Dict(:shape_a=>aa,:shape_b=>bb))
        mynw = nw(opq)
        n, w = gaussjacobi(n,aa,bb)
        @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
        @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
    end
end

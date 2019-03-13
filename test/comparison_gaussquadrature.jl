using PolyChaos, GaussQuadrature, Test, LinearAlgebra

mytol = 1e-10
N = 5:10:100

names = ["hermite", "legendre"]
funs = [hermite, legendre, laguerre]

@testset "Compare $(names...) against GaussQuadrature" begin
    for i=1:length(names), n in N
        opq = OrthoPolyQ(names[i],n)
        mynw = nw(opq)
        n, w = funs[i](n)
        @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
        @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
    end
end

a = -0.5:0.1:2.0
b = copy(a)

@testset "Compare Jacobi against GaussQuadrature" begin
    for aa in a, bb in b, n in N
        opq = OrthoPolyQ("jacobi",n,Dict(:shape_a=>aa,:shape_b=>bb))
        mynw = nw(opq)
        n, w = jacobi(n,aa,bb)
        @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
        @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
    end
end

@testset "Compare Laguerre against GaussQuadrature" begin
    for n in N
        opq = OrthoPolyQ("laguerre",n,Dict(:shape=>0.,:rate=>1.))
        mynw = nw(opq)
        n, w = laguerre(n,0.)
        @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
        @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
    end
end

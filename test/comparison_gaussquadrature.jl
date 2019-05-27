using PolyChaos, GaussQuadrature, Test, LinearAlgebra

@testset "Test Fejer rule" begin
    for m = 1:20
        n, w = fejer(m)
        @test length(n) == length(w) == m
    end
    for m = 1:20
        n, w = golubwelsch(rm_hermite(m+1)...)
        @test length(n) == length(w) == m
    end
end

@testset "N+1 quadrature rules" begin
    @testset "Test Fejer2 rule" begin
        for m = 2:20
            n, w = fejer2(m)
            @test length(n) == length(w) == m + 1
        end
    end
    @testset "Test Clenshaw Curtis rule" begin
        for m = 2:20
            n, w = clenshaw_curtis(m)
            @test length(n) == length(w) == m + 1
        end
    end
    @testset "Test Radau rule" begin
        for m = 1:20
            n, w = radau(rm_hermite(m+2)...,4.)
            @test length(n) == length(w) == m + 1
        end
    end
end

@testset "N+2 quadrature rule (lobatto)" begin
    for m = 1:20
        n, w = lobatto(rm_legendre(m+3)...,-1.,1.)
        @test length(n) == length(w) == m + 2
    end
end




# mytol = 1e-10
# N = 5:10:100
#
# names = ["hermite", "legendre"]
# funs = [hermite, legendre, laguerre]
#
# @testset "Compare $(names...) against GaussQuadrature" begin
#     for i=1:length(names), n in N
#         opq = OrthoPolyQ(names[i],n)
#         mynw = nw(opq)
#         n, w = funs[i](n)
#         @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
#         @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
#     end
# end
#
# a = -0.5:0.1:2.0
# b = copy(a)
#
# @testset "Compare Jacobi against GaussQuadrature" begin
#     for aa in a, bb in b, n in N
#         opq = OrthoPolyQ("jacobi",n,Dict(:shape_a=>aa,:shape_b=>bb))
#         mynw = nw(opq)
#         n, w = jacobi(n,aa,bb)
#         @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
#         @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
#     end
# end
#
# @testset "Compare Laguerre against GaussQuadrature" begin
#     for n in N
#         opq = OrthoPolyQ("laguerre",n,Dict(:shape=>0.,:rate=>1.))
#         mynw = nw(opq)
#         n, w = laguerre(n,0.)
#         @test isapprox(norm(n-mynw[:,1],Inf),0;atol=mytol)
#         @test isapprox(norm(w-mynw[:,2],Inf),0;atol=mytol)
#     end
# end

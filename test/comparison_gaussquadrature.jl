using PolyChaos, GaussQuadrature, Test, LinearAlgebra

@testset "Test Fejer rule" begin
    for m in 1:20
        n, w = fejer(m)
        @test length(n) == length(w) == m
    end
    for m in 1:20
        n, w = golubwelsch(rm_hermite(m + 1)...)
        @test length(n) == length(w) == m
    end
end

@testset "N+1 quadrature rules" begin
    @testset "Test Fejer2 rule" begin
        for m in 2:20
            n, w = fejer2(m)
            @test length(n) == length(w) == m + 1
        end
    end
    @testset "Test Clenshaw Curtis rule" begin
        for m in 2:20
            n, w = clenshaw_curtis(m)
            @test length(n) == length(w) == m + 1
        end
    end
    @testset "Test Radau rule" begin
        for m in 1:20
            n, w = radau(rm_hermite(m + 2)..., 4.0)
            @test length(n) == length(w) == m + 1
        end
    end
end

@testset "N+2 quadrature rule (lobatto)" begin
    for m in 1:20
        n, w = lobatto(rm_legendre(m + 3)..., -1.0, 1.0)
        @test length(n) == length(w) == m + 2
    end
end

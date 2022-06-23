using PolyChaos, Test, QuadGK, StaticArrays
import LinearAlgebra: norm

names = [GaussOrthoPoly, Uniform01OrthoPoly, Uniform_11OrthoPoly, LogisticOrthoPoly,
         LegendreOrthoPoly, HermiteOrthoPoly, LaguerreOrthoPoly]
deg = 10
tol = 1e-7

@testset "Norms of orthgonal polynomials ($names)" begin
    for op in names
        opq = op(deg)
        s = computeSP2(opq)
        t = [quadgk(x -> evaluate(d, x, opq)^2 * opq.measure.w(x), opq.measure.dom...;
                    order = max(10, 10d))[1] for d in 0:deg]
        u = computeSP2(deg, opq)
        @test computeSP2(deg, opq) == computeSP2(deg, opq.β) == computeSP2(deg, SVector(opq.β...))

        @test isapprox(norm(s - t, Inf) / maximum(s), 0; atol = tol)
        @test isapprox(norm(s - u, Inf) / maximum(s), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0], opq), t[1]; atol = tol)
        @test isapprox(computeSP([0, 0, 0, 1], opq), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0, deg - 1, deg - 1], opq), s[deg]; atol = tol)
        @test isapprox(computeSP(SVector([0, 0, 0, deg - 1, deg - 1]...), opq.α, SVector(opq.β...),
                                 opq.quad.nodes, SVector(opq.quad.weights...)), s[deg]; atol = tol)
    end
end

a, b = 1.1:0.5:4, 1.1:0.5:4
names = [Beta01OrthoPoly, JacobiOrthoPoly, genHermiteOrthoPoly]

@testset "Norms of orthgonal polynomials ($names)" begin
    for α in a, β in b, name in names
        opq = if name in [Beta01OrthoPoly, JacobiOrthoPoly]
            name(deg, α, β)
        elseif name in [genHermiteOrthoPoly]
            name(deg, α)
        end

        s = computeSP2(opq)
        t = [quadgk(x -> evaluate(d, x, opq)^2 * opq.measure.w(x), opq.measure.dom...;
                    order = max(10, 10d))[1] for d in 0:deg]
        u = computeSP2(deg, opq)
        @test computeSP2(deg, opq) == computeSP2(deg, opq.β) == computeSP2(deg, SVector(opq.β...))
        @test isapprox(norm(s - t, Inf) / maximum(s), 0; atol = tol)
        @test isapprox(norm(s - u, Inf) / maximum(s), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0], opq), t[1]; atol = tol)
        @test isapprox(computeSP([0, 0, 0, 1], opq), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0, deg - 1, deg - 1], opq), s[deg]; atol = tol)
        @test isapprox(computeSP(SVector([0, 0, 0, deg - 1, deg - 1]...), opq.α, SVector(opq.β...),
                                 opq.quad.nodes, SVector(opq.quad.weights...)), s[deg]; atol = tol)
    end
end

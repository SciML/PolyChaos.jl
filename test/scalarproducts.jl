using PolyChaos, Test, QuadGK
import LinearAlgebra: norm

names = ["gaussian", "uniform01", "logistic", "legendre", "hermite", "laguerre"]
deg = 10
tol = 1e-7

@testset "Norms of orthgonal polynomials ($names)" begin
    for name in names
        # display(name)
        opq = OrthoPolyQ(name,deg)
        s = computeSP2(opq)
        t = [ quadgk(x->evaluate(d,x,opq)^2*opq.op.meas.w(x), opq.op.meas.dom...; order = max(10,10d))[1] for d=0:deg ]
        u = computeSP2(deg,opq.op)
        # display([s t])
        @test isapprox(norm(s-t,Inf)/maximum(s), 0; atol = tol)
        @test isapprox(norm(s-u,Inf)/maximum(s), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0], opq), t[1]; atol = tol)
        @test isapprox(computeSP([0, 0, 0, 1], opq), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0, deg-1, deg-1], opq), s[deg]; atol = tol)
    end
end

a, b = 1.1:0.5:4, 1.1:0.5:4
names = ["beta01", "jacobi", "genhermite", "genlaguerre"]

@testset "Norms of orthgonal polynomials ($names)" begin
    for α in a, β in b, name in names
        if name == "beta01" || name == "jacobi"
            pars = Dict(:shape_a=>α,:shape_b=>β)
        elseif name == "genhermite"
            pars = Dict(:mu=>α)
        elseif name == "genlaguerre"
            pars = Dict(:shape=>α)
        end
        opq = OrthoPolyQ(name,deg,pars)
        s = computeSP2(opq)
        t = [ quadgk(x->evaluate(d,x,opq)^2*opq.op.meas.w(x), opq.op.meas.dom...; order = max(10,10d))[1] for d=0:deg ]
        u = computeSP2(deg,opq.op)
        @test isapprox(norm(s-t,Inf)/maximum(s), 0; atol = tol)
        @test isapprox(norm(s-u,Inf)/maximum(s), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0], opq), t[1]; atol = tol)
        @test isapprox(computeSP([0, 0, 0, 1], opq), 0; atol = tol)
        @test isapprox(computeSP([0, 0, 0, deg-1, deg-1], opq), s[deg]; atol = tol)
    end
end

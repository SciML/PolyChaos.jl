using PolyChaos, Test, QuadGK
import LinearAlgebra: norm

names = ["gaussian", "uniform01", "logistic"]
deg = 10
tol = 1e-7

@testset "Norms of orthgonal polynomials ($names)" begin
    for name in names
        # display(name)
        opq = OrthoPolyQ(name,deg)
        s = computeSP2(opq)
        t = [ quadgk(x->evaluate(d,x,opq)^2*opq.op.meas.w(x), opq.op.meas.dom...; order = max(10,10d))[1] for d=0:deg ]
        # display([s t])
        @test isapprox(norm(s-t,Inf)/maximum(s), 0; atol = tol)
    end
end

a, b = 1.1:0.5:4, 1.1:0.5:4
names = ["beta01", "jacobi"]

@testset "Norms of orthgonal polynomials ($names)" begin
    for α in a, β in b, name in names
        opq = OrthoPolyQ(name,deg,Dict(:shape_a=>α,:shape_b=>β))
        s = computeSP2(opq)
        t = [ quadgk(x->evaluate(d,x,opq)^2*opq.op.meas.w(x), opq.op.meas.dom...; order = max(10,10d))[1] for d=0:deg ]
        @test isapprox(norm(s-t,Inf)/maximum(s), 0; atol = tol)
    end
end

using Test, PolyChaos

d, nunc = 5, 10
op = OrthoPoly("genlaguerre",d,Dict(:shape=>1.23))
opq = OrthoPolyQ(op)
mop = MultiOrthoPoly([op for i in 1:nunc],d)
mopq = MultiOrthoPoly([opq for i in 1:nunc],d)
@testset "dimenions" begin
    @test dim(op) == d + 1
    @test dim(opq) == d + 1
    @test dim(mop) == factorial(d + nunc) / (factorial(d) * factorial(nunc)) == dim(mopq)
end

# @test deg(op) == deg(opq) == deg(mop)

NW = [ opq.quad.nodes opq.quad.weights ]
@testset "nodes and weights" begin
    @test nw(opq) == NW
    @test nw([opq for i in 1:nunc])[1] == [NW[:,1] for i in 1:nunc]
    @test nw([opq for i in 1:nunc])[2] == [NW[:,2] for i in 1:nunc]
    @test nw(mopq)[1] == [NW[:,1] for i in 1:nunc]
    @test nw(mopq)[2] == [NW[:,2] for i in 1:nunc]
end

c = [ op.α op.β ]
@testset "coefficients" begin
    @test coeffs(op) == c
    @test coeffs(opq) == c
    @test coeffs([op for i in 1:nunc])[1] == [ c[:,1] for i in 1:nunc ]
    @test coeffs([op for i in 1:nunc])[2] == [ c[:,2] for i in 1:nunc ]
    @test coeffs([opq for i in 1:nunc])[1] == [ c[:,1] for i in 1:nunc ]
    @test coeffs([opq for i in 1:nunc])[2] == [ c[:,2] for i in 1:nunc ]
    @test coeffs(mop)[1] == [ c[:,1] for i in 1:nunc ]
    @test coeffs(mop)[2] == [ c[:,2] for i in 1:nunc ]
end

@test issymmetric(op) == issymmetric(opq) == issymmetric(opq.quad)

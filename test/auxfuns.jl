using Test, PolyChaos

d, nunc = 5, 10
op = genLaguerreOrthoPoly(d,1.23,addQuadrature = false)
opq = genLaguerreOrthoPoly(d,1.23)
mop = MultiOrthoPoly([op for i in 1:nunc],d)
mopq = MultiOrthoPoly([opq for i in 1:nunc],d)
@testset "dimenions" begin
    @test isequal(PolyChaos.dim(op), d + 1)
    @test isequal(PolyChaos.dim(opq), d + 1)
    @test isequal(PolyChaos.dim(mop), factorial(d + nunc) / (factorial(d) * factorial(nunc)))
end

@testset "degrees" begin
    @test isequal(PolyChaos.deg(op),d)
    @test isequal(PolyChaos.deg(opq),d)
    @test isequal(PolyChaos.deg(mop),d)
end

NW = [ opq.quad.nodes opq.quad.weights ]
@testset "nodes and weights" begin
    @test PolyChaos.nw(opq) == NW
    @test PolyChaos.nw([opq for i in 1:nunc])[1] == [NW[:,1] for i in 1:nunc]
    @test PolyChaos.nw([opq for i in 1:nunc])[2] == [NW[:,2] for i in 1:nunc]
    @test PolyChaos.nw(mopq)[1] == [NW[:,1] for i in 1:nunc]
    @test PolyChaos.nw(mopq)[2] == [NW[:,2] for i in 1:nunc]
end

c = [ op.α op.β ]
@testset "coefficients" begin
    @test PolyChaos.coeffs(op) == c
    @test PolyChaos.coeffs(opq) == c
    @test PolyChaos.coeffs([op for i in 1:nunc])[1] == [ c[:,1] for i in 1:nunc ]
    @test PolyChaos.coeffs([op for i in 1:nunc])[2] == [ c[:,2] for i in 1:nunc ]
    @test PolyChaos.coeffs([opq for i in 1:nunc])[1] == [ c[:,1] for i in 1:nunc ]
    @test PolyChaos.coeffs([opq for i in 1:nunc])[2] == [ c[:,2] for i in 1:nunc ]
    @test PolyChaos.coeffs(mop)[1] == [ c[:,1] for i in 1:nunc ]
    @test PolyChaos.coeffs(mop)[2] == [ c[:,2] for i in 1:nunc ]
end

@test issymmetric(op) == issymmetric(opq) 

using Test, PolyChaos

@test nw(EmptyQuad()) |> typeof == Array{Float64,2}
@test nw(EmptyQuad()) |> length == 0

nodes, weights = gauss(rm_legendre01(5)...)
myQuad = Quad("myQuad", length(nodes), nodes, weights)
foo(x) = 6 * x^5
res = 1
@test isapprox(integrate(foo, nodes, weights), res)
@test isapprox(integrate(foo, myQuad), res)
@test_throws DomainError integrate(foo, EmptyQuad())
@test isapprox(integrate(foo, Uniform01OrthoPoly(4)), res)

foo(x) = 5 * x^4
res = 1
@test isapprox(integrate(foo, Uniform_11OrthoPoly(4)), res)

d, nunc = 5, 10
op = genLaguerreOrthoPoly(d, 1.23; addQuadrature = false)
opq = genLaguerreOrthoPoly(d, 1.23)
mop = MultiOrthoPoly([op for i in 1:nunc], d)
mopq = MultiOrthoPoly([opq for i in 1:nunc], d)
@testset "dimenions" begin
    @test isequal(PolyChaos.dim(op), d + 1)
    @test isequal(PolyChaos.dim(opq), d + 1)
    @test isequal(PolyChaos.dim(mop), factorial(d + nunc) / (factorial(d) * factorial(nunc)))
end

@testset "degrees" begin
    @test isequal(PolyChaos.deg(op), d)
    @test isequal(PolyChaos.deg(opq), d)
    @test isequal(PolyChaos.deg(mop), d)
end

NW = [opq.quad.nodes opq.quad.weights]
@testset "nodes and weights" begin
    @test PolyChaos.nw(opq) == NW
    @test PolyChaos.nw([opq for i in 1:nunc])[1] == [NW[:, 1] for i in 1:nunc]
    @test PolyChaos.nw([opq for i in 1:nunc])[2] == [NW[:, 2] for i in 1:nunc]
    @test PolyChaos.nw(mopq)[1] == [NW[:, 1] for i in 1:nunc]
    @test PolyChaos.nw(mopq)[2] == [NW[:, 2] for i in 1:nunc]
end

c = [op.α op.β]
@testset "coefficients" begin
    @test PolyChaos.coeffs(op) == c
    @test PolyChaos.coeffs(opq) == c
    @test PolyChaos.coeffs([op for i in 1:nunc])[1] == [c[:, 1] for i in 1:nunc]
    @test PolyChaos.coeffs([op for i in 1:nunc])[2] == [c[:, 2] for i in 1:nunc]
    @test PolyChaos.coeffs([opq for i in 1:nunc])[1] == [c[:, 1] for i in 1:nunc]
    @test PolyChaos.coeffs([opq for i in 1:nunc])[2] == [c[:, 2] for i in 1:nunc]
    @test PolyChaos.coeffs(mop)[1] == [c[:, 1] for i in 1:nunc]
    @test PolyChaos.coeffs(mop)[2] == [c[:, 2] for i in 1:nunc]
end

@test issymmetric(op) == issymmetric(opq)

## Accounting for StaticArrays
using StaticArrays

max_degree = 4

op = Uniform01OrthoPoly(max_degree; Nrec = max_degree + 3)
nodes, weights = SVector(op.quad.nodes...), SVector(op.quad.weights...)

myfun(t) = t^2

@test integrate(myfun, op) == integrate(myfun, nodes, weights)

n_multi = 4
mop = MultiOrthoPoly([op for i in 1:n_multi], max_degree)
L = size(mop.ind)[1]

test_vector = [rand(0:(L - 1)) for i in 0:L]
reference = multi2uni(test_vector, mop.ind)

@testset "multi2uni for StaticArrays" begin
    @test multi2uni(SVector(test_vector...), mop.ind) == reference
    @test multi2uni(test_vector, SMatrix{size(mop.ind)...}(mop.ind)) == reference
    @test multi2uni(SVector(test_vector...), SMatrix{size(mop.ind)...}(mop.ind)) == reference
end

tensor_dim = 2
tensor_entries = computeTensorizedSP(tensor_dim, mop)
my_getentry(row, entries) = getentry(row, tensor_entries, entries, tensor_dim)

if VERSION >= VersionNumber("1.1.1")
    @testset "getentry for StaticArrays" begin
        for row in eachrow(rand(0:(L - 1), L, tensor_dim))
            reference = my_getentry(Vector(row), mop.ind)
            @test my_getentry(row, mop.ind) == reference
            @test my_getentry(row, SMatrix{size(mop.ind)...}(mop.ind)) == reference
            @test my_getentry(Vector(row), SMatrix{size(mop.ind)...}(mop.ind)) == reference
        end
    end
end

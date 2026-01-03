using Test, PolyChaos

@test nw(EmptyQuad()) |> typeof == Array{Float64, 2}
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
op = genLaguerreOrthoPoly(d, 1.23, addQuadrature = false)
opq = genLaguerreOrthoPoly(d, 1.23)
mop = MultiOrthoPoly([op for i in 1:nunc], d)
mopq = MultiOrthoPoly([opq for i in 1:nunc], d)
@testset "dimensions" begin
    @test isequal(PolyChaos.dim(op), d + 1)
    @test isequal(PolyChaos.dim(opq), d + 1)
    @test isequal(PolyChaos.dim(mop),
        factorial(d + nunc) / (factorial(d) * factorial(nunc)))
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

@test integrate(myfun, op) ≈ integrate(myfun, nodes, weights)

n_multi = 4
mop = MultiOrthoPoly([op for i in 1:n_multi], max_degree)
L = size(mop.ind)[1]

test_vector = [rand(0:(L - 1)) for i in 0:L]
reference = multi2uni(test_vector, mop.ind)

@testset "multi2uni for StaticArrays" begin
    @test multi2uni(SVector(test_vector...), mop.ind) == reference
    @test multi2uni(test_vector, SMatrix{size(mop.ind)...}(mop.ind)) == reference
    @test multi2uni(SVector(test_vector...), SMatrix{size(mop.ind)...}(mop.ind)) ==
          reference
end

tensor_dim = 2
tensor_entries = computeTensorizedSP(tensor_dim, mop)
my_getentry(row, entries) = getentry(row, tensor_entries, entries, tensor_dim)

@testset "getentry for StaticArrays" begin
    for row in eachrow(rand(0:(L - 1), L,
        tensor_dim))
        reference = my_getentry(Vector(row), mop.ind)
        @test my_getentry(row, mop.ind) == reference
        @test my_getentry(row, SMatrix{size(mop.ind)...}(mop.ind)) == reference
        @test my_getentry(Vector(row), SMatrix{size(mop.ind)...}(mop.ind)) == reference
    end
end

@testset "multivariate integration with MultiOrthoPoly (issue #48)" begin
    # 2D case: Gaussian x Uniform(0,1)
    op_gauss = GaussOrthoPoly(4)
    op_unif = Uniform01OrthoPoly(4)
    mop2 = MultiOrthoPoly([op_gauss, op_unif], 4)

    # Constant function should integrate to 1.0 (normalized measures)
    @test isapprox(integrate((x, y) -> 1.0, mop2), 1.0)

    # E[X^2] * E[Y] where X ~ N(0,1), Y ~ U(0,1)
    # E[X^2] = 1, E[Y] = 0.5
    @test isapprox(integrate((x, y) -> x^2 * y, mop2), 0.5; atol = 1e-10)

    # E[X] * E[Y^2] where E[X] = 0, E[Y^2] = 1/3
    # Result should be 0
    @test isapprox(integrate((x, y) -> x * y^2, mop2), 0.0; atol = 1e-14)

    # 3D case from original issue
    op1 = GaussOrthoPoly(3, addQuadrature = true)
    op2 = Uniform01OrthoPoly(5, addQuadrature = true)
    op3 = Beta01OrthoPoly(6, 2, 1.2, addQuadrature = true)
    mop3 = MultiOrthoPoly([op1, op2, op3], 3)

    # Constant function should integrate to 1.0
    @test isapprox(integrate((a, b, c) -> 1.0, mop3), 1.0)

    # Compare with nested univariate integration
    nested_result = integrate(c -> integrate(b -> integrate(a -> a * b * c, op1), op2), op3)
    mop_result = integrate((a, b, c) -> a * b * c, mop3)
    @test isapprox(mop_result, nested_result; atol = 1e-14)
end

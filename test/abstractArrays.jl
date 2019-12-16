using PolyChaos, StaticArrays, Test

max_degree = 4

op = Uniform01OrthoPoly(max_degree; Nrec = max_degree + 3)
nodes, weights = SVector(op.quad.nodes...), SVector(op.quad.weights...)

myfun(t) = t^2

@test integrate(myfun, op) == integrate(myfun, nodes, weights)

n_multi = 4
mop = MultiOrthoPoly([op for i in 1:n_multi], max_degree)
L = size(mop.ind)[1]

test_vector = [ rand(0:L-1) for i in 0:L ]
reference = multi2uni(test_vector, mop.ind)

@testset "multi2uni for StaticArrays" begin
    @test multi2uni(SVector(test_vector...), mop.ind) == reference
    @test multi2uni(test_vector, SMatrix{size(mop.ind)...}(mop.ind)) == reference
    @test multi2uni(SVector(test_vector...), SMatrix{size(mop.ind)...}(mop.ind)) == reference
end

tensor_dim = 2
tensor_entries = computeTensorizedSP(tensor_dim, mop)
my_getentry(row, entries) = getentry(row, tensor_entries, entries, tensor_dim)

@testset "getentry for StaticArrays" begin
    for row in eachrow(rand(0:L-1, L, tensor_dim))
        reference = my_getentry(Vector(row), mop.ind)
        @test my_getentry(row, mop.ind) == reference
        @test my_getentry(row, SMatrix{size(mop.ind)...}(mop.ind)) == reference
        @test my_getentry(Vector(row), SMatrix{size(mop.ind)...}(mop.ind)) == reference
    end
end

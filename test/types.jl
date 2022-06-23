using Test, PolyChaos
import LinearAlgebra: norm, dot

##########################################################
# MEASURE
##########################################################
@testset "Error handling for Measure" begin
    @test_throws DomainError JacobiMeasure(-1, 2.0)
    @test_throws DomainError JacobiMeasure(-1.0, -2.0)
    @test_throws DomainError JacobiMeasure(1, -2)

    @test_throws DomainError Beta01Measure(0.0, 2.0)
    @test_throws DomainError Beta01Measure(0.0, -2.0)
    @test_throws DomainError Beta01Measure(1, -2)

    @test_throws MethodError Measure("something")

    @test_throws DomainError Measure("test", t -> t, (3, -1), true)
end

##########################################################
# ORTHOPOLY
##########################################################
# @test_throws ErrorException OrthoPoly("idonotexist",3)

quads = [fejer, fejer2, clenshaw_curtis]
discs = [lanczos, stieltjes]
degs = 1:5:20
m = Uniform01Measure()
tol = 1e-6
@testset "General constructor for Gaussian OrthoPoly" begin for quad in quads,
                                                                disc in discs, deg in degs

    op = Uniform01OrthoPoly(deg)
    op1 = OrthoPoly("myuniform", deg, m; quadrature = quad, discretization = disc,
                    Nquad = 100 * deg)
    @test isapprox(norm(coeffs(op) - coeffs(op1), Inf), 0; atol = tol)
    op2 = OrthoPoly("myuniform", deg, t -> 1, (0.0, 1.0), true; quadrature = quad,
                    discretization = disc, Nquad = 100 * deg)
    @test isapprox(norm(coeffs(op) - coeffs(op2), Inf), 0; atol = tol)
end end

##########################################################
# QUAD
##########################################################
Ns = 10:5:30
m = Uniform01Measure()
myf(t) = 6t^5
result = 1.0
tol = 1e-8
@testset "Last-resort-constructor for Quad using quadgp" begin for quad in quads, N in Ns
    op = OrthoPoly("myOP", N, m)
    NW = nw(Quad(N, m; quadrature = quad))
    @test isapprox(dot(myf.(NW[:, 1]), NW[:, 2]), result; atol = tol)
    NW = nw(Quad(N, t -> 1.0, (0, 1); quadrature = quad))
    @test isapprox(dot(myf.(NW[:, 1]), NW[:, 2]), result; atol = tol)

    NW = nw(Quad(N, op))
    @test isapprox(dot(myf.(NW[:, 1]), NW[:, 2]), result; atol = tol)
    # NW = nw(Quad(N, m.w, op.α, op.β, m.dom, m.symmetric))
    # @test isapprox(dot(myf.(NW[:,1]),NW[:,2]), result; atol=tol)
end end

d = 5

ops = [HermiteOrthoPoly(d),
    Uniform01OrthoPoly(d),
    Beta01OrthoPoly(d, 3.2, 5.34)]
mop = MultiOrthoPoly(ops, d)

is = 0:d
n = length(ops)

x = [1.0, 2.0, 3.0]

@testset "Evaluation of multivariate basis" begin for ind in Iterators.product([collect(0:d)
                                                                                for i in 1:n]...)
    @test isapprox(prod(map(i -> evaluate(ind[i], x[i], ops[i]), 1:n)) -
                   evaluate([ind...], x, mop)[1], 0; atol = tol)
end end

##########################################################
# TENSOR
##########################################################
deg, M = 5, 1:4
opq = LegendreOrthoPoly(deg; Nrec = 3 * deg)
mop = MultiOrthoPoly([
                         GaussOrthoPoly(deg; Nrec = 3 * deg),
                         LogisticOrthoPoly(deg; Nrec = 3 * deg),
                         Uniform01OrthoPoly(deg; Nrec = 3 * deg),
                     ], deg)

# @test_throws AssertionError Tensor(0,mop)
@testset "Check tensor" begin for m in M, op in [opq, mop]
    tensor = Tensor(m, op)
    for ind in Iterators.product([collect(0:deg) for i in 1:m]...)
        @test isapprox(tensor.get(collect(ind)), computeSP(collect(ind), op); atol = tol)
    end
end end

# test that error is thrown in case a quadrature rule is missing
opq = LegendreOrthoPoly(deg; Nrec = 3 * deg, addQuadrature = false)
mop = MultiOrthoPoly([
                         GaussOrthoPoly(deg; Nrec = 3 * deg, addQuadrature = false),
                         LogisticOrthoPoly(deg; Nrec = 3 * deg),
                         Uniform01OrthoPoly(deg; Nrec = 3 * deg),
                     ], deg)

@test_throws InconsistencyError Tensor(3, opq)
@test_throws InconsistencyError Tensor(3, mop)

using Test, PolyChaos
import LinearAlgebra: norm


##########################################################
# MEASURE
##########################################################
@testset "Error handling for Measure" begin
    @test_throws AssertionError Measure("jacobi",Dict(:shape_a=>-1.,:shape_b=>2.))
    @test_throws AssertionError Measure("jacobi",Dict(:shape_a=>-1.,:shape_b=>-2.))
    @test_throws AssertionError Measure("jacobi",Dict(:shape_a=>1,:shape_b=>-2))

    @test_throws AssertionError Measure("beta01",Dict(:shape_a=>0.,:shape_b=>2.))
    @test_throws AssertionError Measure("beta01",Dict(:shape_a=>0.,:shape_b=>-2.))
    @test_throws AssertionError Measure("beta01",Dict(:shape_a=>1,:shape_b=>-2))

    @test_throws ErrorException Measure("something")

    @test_throws AssertionError Measure("test",t->t,(3,-1),true)
end

##########################################################
# ORTHOPOLY
##########################################################
@test_throws ErrorException OrthoPoly("idonotexist",3)

quads = [ fejer, fejer2, clenshaw_curtis ]
discs = [ lanczos, stieltjes ]
degs  = 1:5:20
m = Measure("uniform01")
tol = 1e-6
@testset "General constructor for Gaussian OrthoPoly" begin
    for quad in quads, disc in discs, deg in degs
        op = OrthoPoly("uniform01",deg)
        op1 = OrthoPoly("myuniform",deg,m; quadrature=quad, discretization = disc, Nquad = 100*deg)
        @test isapprox(norm(coeffs(op)-coeffs(op1),Inf),0; atol = tol)
        op2 = OrthoPoly("myuniform",deg,t->1,(0,1),true; quadrature=quad, discretization = disc, Nquad = 100*deg)
        @test isapprox(norm(coeffs(op)-coeffs(op2),Inf),0; atol = tol)
    end
end

##########################################################
# QUAD
##########################################################
Ns = 10:5:30
m = Measure("uniform01")
myf(t) = 6t^5
result = 1.
tol = 1e-8
@testset "Last-resort-constructor for Quad using quadgp" begin
    for quad in quads, N in Ns
        NW = nw(Quad(N,m;quadrature=quad))
        @test isapprox(dot(myf.(NW[:,1]),NW[:,2]), result; atol=tol)
        NW = nw(Quad(N,t->1.,(0,1),true;quadrature=quad))
        @test isapprox(dot(myf.(NW[:,1]),NW[:,2]), result; atol=tol)
    end
end

d = 5
ops = [ OrthoPoly("hermite",d),
        OrthoPoly("uniform01",d),
        OrthoPoly("beta01",d,Dict(:shape_a=>3.2,:shape_b=>5.34)) ]
mop = MultiOrthoPoly(ops,d)

is = 0:d
n = length(ops)

x = [1., 2., 3.]

@testset "Evaluation of multivariate basis" begin
    for ind in Iterators.product([collect(0:d) for i = 1:n]...)
        @test isapprox(prod( map(i->evaluate(ind[i],x[i],ops[i]),1:n) ) -  evaluate([ind...],x,mop)[1],0;atol=tol)

    end
end

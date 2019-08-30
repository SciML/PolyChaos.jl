using PolyChaos, Test

deg = 5
op = genLaguerreOrthoPoly(deg,Dict(:shape=>1.23))
opq = OrthoPolyQ(op)

n = 10
X = [ exp(3 + rand()), exp.(3 .+ rand(n)) ]
@testset "evaluate univariate basis at point(s)" begin
    for x in X
        @test evaluate(collect(0:deg),x,op.α,op.β) == evaluate(collect(0:deg),x,op) == evaluate(collect(0:deg),x,opq)
        @test evaluate(x,op) == evaluate(x,opq)
        for d in 0:deg
            @test evaluate(d,x,op.α,op.β) == evaluate(d,x,op) == evaluate(d,x,opq)
        end
    end
end

ops = [  genLaguerreOrthoPoly(deg,Dict(:shape=>1.23)),
        genHermiteOrthoPoly(,deg,Dict(:mu=>3.4)),
        LegendreOrthoPoly(deg),
        JacobiOrthoPoly(,deg,Dict(:shape_a=>4.3,:shape_b=>10.))
     ]
opqs = map(OrthoPolyQ, ops)
mop = MultiOrthoPoly(ops,deg)
mopq = MultiOrthoPoly(opqs,deg)

nop = length(ops)
X = [ exp.(3 .+ rand(nop)) , exp.(3 .+ rand(n,nop)) ]

@testset "evaluate multivariate basis at point(s)" begin
    for x in X, op in [mop, mopq]
        @test evaluate(op.ind,x,coeffs(op)...) == evaluate(op.ind,x,op) == evaluate(x,op)
        for d in eachrow(op.ind[1:3,:])
            @test evaluate(collect(d),x,coeffs(op)...) == evaluate(collect(d),x,op)
        end
    end
end

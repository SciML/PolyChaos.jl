using PolyChaos, Test

deg = 5
op = genLaguerreOrthoPoly(deg,1.23)

n = 10
X = [ exp(3 + rand()), exp.(3 .+ rand(n)) ]
@testset "evaluate univariate basis at point(s)" begin
    for x in X
        @test evaluate(collect(0:deg),x,op.α,op.β) == evaluate(collect(0:deg),x,op)
        for d in 0:deg
            @test evaluate(d,x,op.α,op.β) == evaluate(d,x,op)
        end
    end
end

ops = [  genLaguerreOrthoPoly(deg,1.23),
        genHermiteOrthoPoly(deg,3.),
        LegendreOrthoPoly(deg),
        JacobiOrthoPoly(deg,4.3,10.)
     ]

mop = MultiOrthoPoly(ops,deg)


nop = length(ops)
X = [ exp.(3 .+ rand(nop)) , exp.(3 .+ rand(n,nop)) ]

@testset "evaluate multivariate basis at point(s)" begin
    for x in X
        @test evaluate(mop.ind,x,coeffs(mop)...) == evaluate(mop.ind,x,mop) == evaluate(x,mop)
        for d in eachrow(mop.ind[1:3,:])
            @test evaluate(collect(d),x,coeffs(mop)...) == evaluate(collect(d),x,mop)
        end
    end
end

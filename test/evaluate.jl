using PolyChaos, Test

deg = 5
op = genLaguerreOrthoPoly(deg,1.23)

n = 10
X = [ exp(3 + rand()), exp.(3 .+ rand(n)) ]

@testset "evaluate univariate basis at point(s)" begin
    for x in X
        @test evaluate(collect(0:deg),x,op.α,op.β) == evaluate(collect(0:deg),x,op)
        for d in 0:deg
            if typeof(x) <: Real
                @test evaluate(d,x,op.α,op.β) == evaluate(d,x,op) == evaluate(d,[x],op.α,op.β)
            else
                @test evaluate(d,x,op.α,op.β) == evaluate(d,x,op)
            end
        end

        res = hcat(map(i->evaluate(i,x,op.α,op.β),collect(0:op.deg))...)
        if typeof(x) <: Real
            @test evaluate(x,op) == evaluate([x], op) == evaluate(collect(0:op.deg),[x],op) == evaluate(collect(0:op.deg),[x],op.α,op.β) == res
        else
            @test evaluate(x, op) == evaluate(collect(0:op.deg),x,op) == evaluate(collect(0:op.deg),x,op.α,op.β) == res
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
        if typeof(x) <: Vector{<:Real}
            @test evaluate(mop.ind,x,coeffs(mop)...) == evaluate(mop.ind,x,mop) == evaluate(x,mop) == evaluate(mop.ind,reshape(x,1,length(x)),coeffs(mop)...)
            @test evaluate(x,mop) == evaluate(mop.ind,reshape(x,1,length(x)),mop) == evaluate(mop.ind,reshape(x,1,length(x)),coeffs(mop)...)
        else
            @test evaluate(x,mop) == evaluate(mop.ind,x,mop) == evaluate(mop.ind,x,coeffs(mop)...)
        end
        for d in eachrow(mop.ind[1:3,:])
            @test evaluate(collect(d),x,coeffs(mop)...) == evaluate(collect(d),x,mop)
        end
    end
end

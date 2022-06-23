using PolyChaos, Test, StaticArrays

deg = 5
op = genLaguerreOrthoPoly(deg, 1.23)

n = 10
X = [exp(3 + rand()), exp.(3 .+ rand(n))]

α, β, αs, βs = op.α, op.β, SVector(op.α...), SVector(op.β...)
combinations = Iterators.product([α, αs], [β, βs])

@testset "evaluate univariate basis at point(s)" begin
    for x in X
        [@test evaluate(collect(0:deg), x, op) == evaluate(collect(0:deg), x, a, b)
         for (a, b) in combinations]
        [@test evaluate(SVector(0:deg...), x, op) == evaluate(SVector(0:deg...), x, a, b)
         for (a, b) in combinations]
        [@test evaluate(0:deg, x, op) == evaluate(0:deg, x, a, b) for (a, b) in combinations]
        [@test evaluate(0:2:deg, x, op) == evaluate(0:2:deg, x, a, b) for (a, b) in combinations]
        for d in 0:deg
            if typeof(x) <: Real
                [@test evaluate(d, x, a, b) == evaluate(d, x, op) == evaluate(d, [x], a, b)
                 for (a, b) in combinations]
            else
                [@test evaluate(d, x, a, b) == evaluate(d, x, op) for (a, b) in combinations]
            end
        end

        res = hcat(map(i -> evaluate(i, x, op.α, op.β), collect(0:(op.deg)))...)
        if typeof(x) <: Real
            [@test evaluate(x, op) == evaluate([x], op) == evaluate(collect(0:(op.deg)), [x], op) ==
                   evaluate(collect(0:(op.deg)), [x], a, b) == res for (a, b) in combinations]
            [@test evaluate(x, op) == evaluate([x], op) == evaluate(0:(op.deg), [x], op) ==
                   evaluate(0:(op.deg), [x], a, b) == res for (a, b) in combinations]
        else
            [@test evaluate(x, op) == evaluate(collect(0:(op.deg)), x, op) ==
                   evaluate(collect(0:(op.deg)), x, a, b) == res for (a, b) in combinations]
            [@test evaluate(x, op) == evaluate(0:(op.deg), x, op) ==
                   evaluate(0:(op.deg), x, a, b) == res for (a, b) in combinations]
        end
    end
end

ops = [genLaguerreOrthoPoly(deg, 1.23),
       genHermiteOrthoPoly(deg, 3.0),
       LegendreOrthoPoly(deg),
       JacobiOrthoPoly(deg, 4.3, 10.0)]

mop = MultiOrthoPoly(ops, deg)

nop = length(ops)
X = [exp.(3 .+ rand(nop)), exp.(3 .+ rand(n, nop))]

α, β = coeffs(mop)
αs, βs = SVector(α...), SVector(β...)
combinations = Iterators.product([α, αs], [β, βs])

@testset "evaluate multivariate basis at point(s)" begin
    for x in X
        [@test evaluate(mop.ind, x, a, b) == evaluate(mop.ind, x, mop) == evaluate(x, mop)
         for (a, b) in combinations]
        if typeof(x) <: AbstractVector{<:Real}
            [@test evaluate(mop.ind, x, a, b) == evaluate(mop.ind, x, mop) == evaluate(x, mop) ==
                   evaluate(mop.ind, reshape(x, 1, length(x)), coeffs(mop)...)
             for (a, b) in combinations]
            [@test evaluate(x, mop) == evaluate(mop.ind, reshape(x, 1, length(x)), mop) ==
                   evaluate(mop.ind, reshape(x, 1, length(x)), a, b) for (a, b) in combinations]
        else
            [@test evaluate(x, mop) == evaluate(mop.ind, x, mop) == evaluate(mop.ind, x, a, b)
             for (a, b) in combinations]
        end
        if VERSION >= VersionNumber("1.1.1")
            for d in eachrow(mop.ind[1:3, :])
                [@test evaluate(collect(d), x, a, b) == evaluate(collect(d), x, mop)
                 for (a, b) in combinations]
                [@test evaluate(d, x, a, b) == evaluate(collect(d), x, mop)
                 for (a, b) in combinations]
            end
        end
    end
end

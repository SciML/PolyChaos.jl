using PolyChaos, Test, StaticArrays
import LinearAlgebra: norm
myfile = open("dataQuadratureRules/config.txt")
ns = parse.(Int, readlines(myfile))
close(myfile)
myfile = open("dataQuadratureRules/endPts.txt")
endPts = parse.(Float64, readlines(myfile))
close(myfile)
nodes = ns[1]:ns[2]:ns[3]
high = ns[3]
tol = 1e-7

@time @testset "Fejer" begin
    for n in nodes
        myfile = open("dataQuadratureRules/fejer$n.txt")
        αβref = parse.(Float64, readlines(myfile))
        αβcom = fejer(n)
        @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
        close(myfile)
    end
end

@time @testset "Fejer2" begin
    for n in nodes
        myfile = open("dataQuadratureRules/fejer2_$n.txt")
        αβref = parse.(Float64, readlines(myfile))
        αβcom = fejer2(n)
        @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
        close(myfile)
    end
end

@time @testset "Clenshaw Curtis" begin
    for n in nodes
        myfile = open("dataQuadratureRules/cc$n.txt")
        αβref = parse.(Float64, readlines(myfile))
        αβcom = clenshaw_curtis(n)
        @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
        close(myfile)
    end
end

@time @testset "Gauss" begin
    data1 = open("dataQuadratureRules/logHigh.txt")
    d1 = reshape(parse.(Float64, readlines(data1)), :, 2)
    data2 = open("dataQuadratureRules/hermiteHigh.txt")
    d2 = reshape(parse.(Float64, readlines(data2)), :, 2)
    close(data1)
    close(data2)
    for n in ns
        myfile = open("dataQuadratureRules/gaussLog$n.txt")
        αβref = parse.(Float64, readlines(myfile))
        αβcom = gauss(n, d1[:, 1], d1[:, 2])
        @test gauss(n, d1[:, 1], d1[:, 2]) == gauss(n, SVector(d1[:, 1]...), SVector(d1[:, 2]...))
        @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
        q = Quad(n, d1[:, 1], d1[:, 2])
        @test isapprox(norm(αβref - vcat(nw(q)...), Inf), 0.0; atol = tol)
        close(myfile)
        myfile = open("dataQuadratureRules/gaussHerm$n.txt")
        αβref = parse.(Float64, readlines(myfile))
        αβcom = gauss(n, d2[:, 1], d2[:, 2])
        @test gauss(n, d2[:, 1], d2[:, 2]) == gauss(n, SVector(d2[:, 1]...), SVector(d2[:, 2]...))
        @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
        q = Quad(n, d2[:, 1], d2[:, 2])
        @test isapprox(norm(αβref - vcat(nw(q)...), Inf), 0.0; atol = tol)
        close(myfile)
    end
end

@time @testset "Radau with Logistic" begin
    data1 = open("dataQuadratureRules/logHigh.txt")
    d1 = reshape(parse.(Float64, readlines(data1)), :, 2)
    close(data1)
    for n in ns
        for i in [endPts; -endPts]
            myfile = open("dataQuadratureRules/radauLog$(n)pt$i.txt")
            αβref = parse.(Float64, readlines(myfile))
            αβcom = radau(n, d1[:, 1], d1[:, 2], i)
            @test radau(n, d1[:, 1], d1[:, 2], i) ==
                  radau(n, MVector(d1[:, 1]...), MVector(d1[:, 2]...), i) ==
                  radau(n, MVector(d1[:, 1]...), d1[:, 2], i) ==
                  radau(n, d1[:, 1], MVector(d1[:, 2]...), i)
            @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
            close(myfile)
        end
    end
end

@time @testset "Radau with Hermite" begin
    data2 = open("dataQuadratureRules/hermiteHigh.txt")
    d2 = reshape(parse.(Float64, readlines(data2)), :, 2)
    close(data2)
    for n in ns
        for i in [endPts; -endPts]
            myfile = open("dataQuadratureRules/radauHerm$(n)pt$i.txt")
            αβref = parse.(Float64, readlines(myfile))
            αβcom = radau(n, d2[:, 1], d2[:, 2], i)
            @test radau(n, d2[:, 1], d2[:, 2], i) ==
                  radau(n, MVector(d2[:, 1]...), MVector(d2[:, 2]...), i) ==
                  radau(n, MVector(d2[:, 1]...), d2[:, 2], i) ==
                  radau(n, d2[:, 1], MVector(d2[:, 2]...), i)
            @test isapprox(norm(αβref - [αβcom[1]; αβcom[2]], Inf), 0.0; atol = tol)
            close(myfile)
        end
    end
end

function testInt(lb, ub)
    # analytic integral of f(x) = exp(-|x|)
    @assert lb < ub "lower bound has to be smaller than upper bound"
    lb > 0 && return -exp(-ub) + exp(-lb)
    ub < 0 && return exp(ub) - exp(lb)
    return 2 - exp(lb) - exp(-ub)
end

ns[3] = 700
intervals = [(-5.0, 5.0), (-Inf, Inf), (-Inf, -1.0), (1.0, Inf), (-3.0, 5.0), (-10.0, 1.0)]
@time @testset "gpquad" begin
    for n in intervals
        nod1, wei1 = quadgp(x -> exp.(-abs.(x)), n[1], n[2], ns[3]; quadrature = fejer)
        #print("Nodes: ",nod1,"\n")
        #print("Weights: ",wei1,"\n")
        analytic = testInt(n[1], n[2])
        @test isapprox(sum(wei1) - analytic, 0; atol = 1e-4)
    end
    for n in intervals
        nod2, wei2 = quadgp(x -> 1.0, n[1], n[2], ns[3]; quadrature = fejer)
        analytic = testInt(n[1], n[2])
        @test isapprox(exp.(-abs.(nod2))' * wei2 - analytic, 0; atol = 1e-4)
    end
end

degree = 9
N = degree + 1
α, β = rm_legendre01(degree + 1)
op = Uniform01OrthoPoly(degree)
endl, endr = 0, 1

@testset "OrthoPoly overloading" begin
    @test gauss(N - 1, α, β) == gauss(α, β) == gauss(N - 1, op) == gauss(op) ==
          gauss(N - 1, SVector(α...), β) == gauss(N - 1, α, SVector(β...))
    @test radau(N - 2, α, β, endr) == radau(α, β, endr) == radau(N - 2, op, endr) == radau(op, endr)
    @test lobatto(N - 3, α, β, endl, endr) == lobatto(α, β, endl, endr) ==
          lobatto(N - 3, op, endl, endr) == lobatto(op, endl, endr) ==
          lobatto(N - 3, MVector(α...), MVector(β...), endl, endr)
end

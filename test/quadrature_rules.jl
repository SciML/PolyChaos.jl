using PolyChaos, Test
import LinearAlgebra: norm
myfile = open("dataQuadratureRules/config.txt")
ns = parse.(Int,readlines(myfile))
close(myfile)
myfile = open("dataQuadratureRules/endPts.txt")
endPts = parse.(Float64,readlines(myfile))
close(myfile)
nodes = ns[1]:ns[2]:ns[3]
high = ns[3]
tol = 1e-7


@time @testset "Fejer" begin
    for n in nodes
            myfile = open("dataQuadratureRules/fejer$n.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = fejer(n)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
    end
end

@time @testset "Fejer2" begin
    for n in nodes
            myfile = open("dataQuadratureRules/fejer2_$n.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = fejer2(n)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
    end
end

@time @testset "Clenshaw Curtis" begin
    for n in nodes
            myfile = open("dataQuadratureRules/cc$n.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = clenshaw_curtis(n)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
    end
end

@time @testset "Gauss" begin
        data1 = open("dataQuadratureRules/logHigh.txt")
        d1 = reshape(parse.(Float64,readlines(data1)),:,2)
        data2 = open("dataQuadratureRules/hermiteHigh.txt")
        d2 = reshape(parse.(Float64,readlines(data2)),:,2)
        close(data1)
        close(data2)
    for n in ns
            myfile = open("dataQuadratureRules/gaussLog$n.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = gauss(n,d1[:,1],d1[:,2])
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
            myfile = open("dataQuadratureRules/gaussHerm$n.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = gauss(n,d2[:,1],d2[:,2])
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
    end
end

@time @testset "Radau with Logistic" begin
        data1 = open("dataQuadratureRules/logHigh.txt")
        d1 = reshape(parse.(Float64,readlines(data1)),:,2)
        close(data1)
    for n in ns
        for i in [endPts;-endPts]
            myfile = open("dataQuadratureRules/radauLog$(n)pt$i.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = radau(n,d1[:,1],d1[:,2],i)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)

        end
    end
end

@time @testset "Radau with Hermite" begin
        data2 = open("dataQuadratureRules/hermiteHigh.txt")
        d2 = reshape(parse.(Float64,readlines(data2)),:,2)
        close(data2)
    for n in ns
        for i in [endPts;-endPts]
            myfile = open("dataQuadratureRules/radauHerm$(n)pt$i.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = radau(n,d2[:,1],d2[:,2],i)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
        end
    end
end

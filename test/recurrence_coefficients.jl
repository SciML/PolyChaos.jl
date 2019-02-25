# Test whether the numerically computed recursion coefficients are correct
using PolyChaos, Test
import LinearAlgebra: norm
# load config
myfile = open("dataRecCoeffs/config.txt")
ns = parse.(Int,readlines(myfile))

nodes = ns[1]:ns[2]:ns[3]
mus = 0.:0.1:1
albe = 0.:0.2:2
tol = 1e-9


@time for n in nodes
        @testset "Hermite $n nodes" begin
            for m in mus
            myfile = open("dataRecCoeffs/hermite$(n)mu$m.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = rm_hermite(n,m)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
        end
    end
end

@time for n in nodes
        @testset "Logistic $n nodes" begin
            myfile = open("dataRecCoeffs/log$(n).txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = rm_logistic(n)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
    end
end

@time for n in nodes
        @testset "Jacobi $n nodes" begin
            for al in albe
                for be in albe
                    myfile = open("dataRecCoeffs/jac$(n)al$(al)be$(be).txt")
                    αβref = parse.(Float64,readlines(myfile))
                    αβcom = rm_jacobi(n,al,be)
                    @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
                end
            end
    end
end

@time for n in nodes
        @testset "Jacobi01 $n nodes" begin
            for al in albe
                for be in albe
                    myfile = open("dataRecCoeffs/jac01$(n)al$(al)be$(be).txt")
                    αβref = parse.(Float64,readlines(myfile))
                    αβcom = rm_jacobi01(n,al,be)
                    @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
                    close(myfile)
                end
            end
    end
end

@time for n in nodes
        @testset "Laguerre $n nodes" begin
            for m in mus
            myfile = open("dataRecCoeffs/laguerre$(n)a$m.txt")
            αβref = parse.(Float64,readlines(myfile))
            αβcom = rm_laguerre(n,m)
            @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
            close(myfile)
        end
    end
end

meipol = 0.1:0.2:2
@time for n in nodes
        @testset "meixner_pollaczek $n nodes" begin
            for lambda in meipol
                for phi in meipol
                    myfile = open("dataRecCoeffs/meixpol$(n)la$(lambda)phi$(phi).txt")
                    αβref = parse.(Float64,readlines(myfile))
                    αβcom = rm_meixner_pollaczek(n,lambda,phi)
                    @test isapprox(norm(αβref-[αβcom[1];αβcom[2]],Inf),0.;atol=tol)
                    close(myfile)
                end
            end
    end
end

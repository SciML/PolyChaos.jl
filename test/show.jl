using PolyChaos, Test, LinearAlgebra

# values taken from https://en.wikipedia.org/wiki/Hermite_polynomials#Definition
c = [
    [-0.0],
    [-1.0, -0.0],
    [0.0, -3.0, -0.0],
    [3.0, 0.0, -6.0, -0.0],
    [-0.0, 15.0, 0.0, -10.0, -0.0],
    [-15.0, -0.0, 45.0, 0.0, -15.0, -0.0],
    [0.0, -105.0, -0.0, 105.0, 0.0, -21.0, -0.0],
    [105.0, 0.0, -420.0, -0.0, 210.0, 0.0, -28.0, -0.0],
    [-0.0, 945.0, 0.0, -1260.0, -0.0, 378.0, 0.0, -36.0, -0.0],
    [-945.0, -0.0, 4725.0, 0.0, -3150.0, -0.0, 630.0, 0.0, -45.0, -0.0],
    [0.0, -10395.0, -0.0, 17325.0, 0.0, -6930.0, -0.0, 990.0, 0.0, -55.0, -0.0],
]
α, β = rm_hermite_prob(length(c))
cc = rec2coeff(α, β)

@testset "showpoly" begin
    for (i, c_) in enumerate(c)
        @test isapprox(norm(c_ - cc[i]), 0)
    end
end

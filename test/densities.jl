using PolyChaos, Test

outliers = [-1.3, 1.3]

@test_throws DomainError PolyChaos._throwError(3)

densities = [w_legendre,
    w_uniform01,
    w_uniform_11,
    build_w_jacobi(1.2, 3.4),
    build_w_beta(1.2, 3.4),
]

@testset "Supports of densities" begin for ρ in densities, x in outliers
    @test_throws DomainError ρ(x)
end end

@test_throws DomainError w_laguerre(-1)
@test_throws DomainError build_w_gamma(1)(-1)

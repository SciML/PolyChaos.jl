using PolyChaos, Test, StaticArrays, Random, QuadGK
import Statistics

function betaMoments(α, β)
    # moments of beta distribution, analytic solution
    return α / (α + β), sqrt(α * β / ((α + β)^2 * (1 + α + β)))
end

@testset "adaptive rejection sampling through sampleMeasure" begin
    weight(x) = exp(-((x - 0.3)^2) / 0.5)
    domain = (-1.0, 1.0)
    normalization, _ = quadgk(weight, domain...)
    target_mean, _ = quadgk(x -> x * weight(x), domain...)
    target_mean /= normalization
    target_variance, _ = quadgk(x -> (x - target_mean)^2 * weight(x), domain...)
    target_variance /= normalization
    target_fourth_central_moment, _ = quadgk(
        x -> (x - target_mean)^4 * weight(x), domain...
    )
    target_fourth_central_moment /= normalization

    Random.seed!(1234)
    samples = sampleMeasure(20_000, weight, domain)

    @test length(samples) == 20_000
    @test all(x -> domain[1] <= x <= domain[2], samples)
    @test abs(Statistics.mean(samples) - target_mean) <
        4 * sqrt(target_variance / length(samples))
    variance_se = sqrt(
        (target_fourth_central_moment - target_variance^2) / length(samples)
    )
    @test abs(Statistics.var(samples) - target_variance) < 4 * variance_se

    Random.seed!(1234)
    infinite_samples = sampleMeasure(1_000, x -> exp(-x^2), (-Inf, Inf))
    @test length(infinite_samples) == 1_000
    @test all(isfinite, infinite_samples)
end

degs, Nsamples = 1:5, 10000

Random.seed!(1234)
α, β = rand():2:10, rand():0.3:7

@testset "Mean and variance of beta distribution" begin
    for a in α, b in β, deg in degs
        op = Beta01OrthoPoly(deg, a, b)

        @test calculateAffinePCE(op) == calculateAffinePCE(op.α) ==
            calculateAffinePCE(SVector(op.α...))

        coeffs = calculateAffinePCE(op)
        coeffs_static = SVector(coeffs...)

        @test (mean(coeffs, op), std(coeffs, op)) ==
            (mean(coeffs_static, op), std(coeffs_static, op))

        mu, sigma = mean(coeffs, op), std(coeffs, op)
        mu_ana, sigma_ana = betaMoments(a, b)
        @test isapprox(mu, mu_ana; atol = 1.0e-5)
        @test isapprox(sigma, sigma_ana; atol = 1.0e-5)

        samples = sampleMeasure(Nsamples, op)
        @test isapprox(mu, mean(samples); atol = 1.0e-2)
        @test isapprox(sigma, std(samples); atol = 1.0e-2)

        evals = evaluatePCE(coeffs, samples, op)

        @test evals == evaluatePCE(SVector(coeffs...), SVector(samples...), op) ==
            evaluatePCE(SVector(coeffs...), samples, op) ==
            evaluatePCE(SVector(coeffs...), samples, op.α, op.β)

        @test isapprox(mu, mean(evals); atol = 1.0e-2)
        @test isapprox(sigma, std(evals); atol = 1.0e-2)
    end
end

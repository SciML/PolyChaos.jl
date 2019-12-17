using PolyChaos, Test, StaticArrays

function betaMoments(α, β)
	# moments of beta distribution, analytic solution
	α / (α + β), sqrt(α*β / ((α+β)^2 * (1 + α + β)))
end

degs, Nsamples = 1:5, 10000

α, β = rand():2:10, rand():0.3:7


@testset "Mean and variance of beta distribution" begin
	for a in α, b in β, deg in degs
		op = Beta01OrthoPoly(deg, a, b)

		@test calculateAffinePCE(op) == calculateAffinePCE(op.α) == calculateAffinePCE(SVector(op.α...))

		coeffs = calculateAffinePCE(op)
		coeffs_static = SVector(coeffs...)

		@test (mean(coeffs, op), std(coeffs, op)) == (mean(coeffs_static, op), std(coeffs_static, op))

		mu, sigma = mean(coeffs, op), std(coeffs, op)
		mu_ana, sigma_ana = betaMoments(a, b)
		@test isapprox(mu, mu_ana; atol=1e-5)
		@test isapprox(sigma, sigma_ana; atol=1e-5)

		samples = sampleMeasure(Nsamples, op)
		@test isapprox(mu, mean(samples); atol=1e-2)
		@test isapprox(sigma, std(samples); atol=1e-2)

		evals = evaluatePCE(coeffs, samples, op)

		@test evals == evaluatePCE(SVector(coeffs...), SVector(samples...), op) == evaluatePCE(SVector(coeffs...), samples, op) == evaluatePCE(SVector(coeffs...), samples, op.α, op.β)

		@test isapprox(mu, mean(evals); atol=1e-2)
		@test isapprox(sigma, std(evals); atol=1e-2)

	end
end

using PolyChaos, Test

@testset "Inverse Transform Sampling" begin
    Nsamples = 5000
    atol_mean = 0.05
    atol_std = 0.1

    @testset "sampleInverseCDF basic functionality" begin
        # Uniform distribution
        pdf_uniform = x -> 1.0
        samples = sampleInverseCDF(1000, pdf_uniform, (0.0, 1.0))
        @test length(samples) == 1000
        @test all(0 .<= samples .<= 1)
        @test isapprox(mean(samples), 0.5; atol = atol_mean)

        # Beta-like distribution
        pdf_beta = x -> x * (1 - x)
        samples = sampleInverseCDF(Nsamples, pdf_beta, (0.0, 1.0))
        @test isapprox(mean(samples), 0.5; atol = atol_mean)

        # Truncated Gaussian
        pdf_gauss = x -> exp(-x^2 / 2)
        samples = sampleInverseCDF(Nsamples, pdf_gauss, (-3.0, 3.0))
        @test isapprox(mean(samples), 0.0; atol = atol_mean)
        @test isapprox(std(samples), 1.0; atol = atol_std)
    end

    @testset "sampleMeasure with inversecdf method" begin
        # Custom Measure
        w = x -> exp(-x^2)
        m = Measure("custom_gaussian", w, (-5.0, 5.0), true)
        samples = sampleMeasure(Nsamples, m; method = "inversecdf")
        @test isapprox(mean(samples), 0.0; atol = atol_mean)

        # Beta01Measure
        m = Beta01Measure(2.0, 5.0)
        samples = sampleMeasure(Nsamples, m; method = "inversecdf")
        expected_mean = 2 / (2 + 5)
        @test isapprox(mean(samples), expected_mean; atol = atol_mean)
    end

    @testset "Canonical measures with inversecdf" begin
        # GaussMeasure (infinite domain)
        m = GaussMeasure()
        samples = sampleMeasure(Nsamples, m; method = "inversecdf")
        @test isapprox(mean(samples), 0.0; atol = atol_mean)
        @test isapprox(std(samples), 1.0; atol = atol_std)

        # Uniform01Measure
        m = Uniform01Measure()
        samples = sampleMeasure(Nsamples, m; method = "inversecdf")
        @test isapprox(mean(samples), 0.5; atol = atol_mean)

        # GammaMeasure (semi-infinite domain)
        m = GammaMeasure(2.0, 1.0)
        samples = sampleMeasure(Nsamples, m; method = "inversecdf")
        @test isapprox(mean(samples), 2.0; atol = atol_mean)
    end

    @testset "OrthoPoly interface with inversecdf" begin
        # Test via AbstractOrthoPoly
        op = Beta01OrthoPoly(3, 3.0, 2.0)
        samples = sampleMeasure(Nsamples, op; method = "inversecdf")
        expected_mean = 3 / (3 + 2)
        @test isapprox(mean(samples), expected_mean; atol = atol_mean)
    end

    @testset "Error handling" begin
        pdf = x -> 1.0
        @test_throws DomainError sampleInverseCDF(0, pdf, (0.0, 1.0))
        @test_throws DomainError sampleInverseCDF(-1, pdf, (0.0, 1.0))
        @test_throws DomainError sampleInverseCDF(10, pdf, (1.0, 0.0))
    end
end

using Test
using PolyChaos

@testset "Interface Compatibility" begin
    @testset "BigFloat support" begin
        # Test that evaluate returns BigFloat when given BigFloat inputs
        op = LegendreOrthoPoly(5)
        x_bf = BigFloat[BigFloat("0.1"), BigFloat("0.5"), BigFloat("0.9")]
        result = evaluate(2, x_bf, op)
        @test eltype(result) == BigFloat
        @test length(result) == 3

        # Test single point evaluation
        x_single = BigFloat("0.5")
        result_single = evaluate(2, x_single, op)
        @test result_single isa BigFloat

        # Test evaluate with multiple degrees
        result_multi = evaluate([0, 1, 2], x_bf, op)
        @test eltype(result_multi) == BigFloat
        @test size(result_multi) == (3, 3)

        # Test evaluatePCE with BigFloat
        x_pce = BigFloat[BigFloat("1.0"), BigFloat("0.5")]
        xi_bf = BigFloat[BigFloat("0.2"), BigFloat("0.4")]
        result_pce = evaluatePCE(x_pce, xi_bf, op)
        @test eltype(result_pce) == BigFloat

        # Test computeSP2 with BigFloat
        beta_bf = BigFloat[BigFloat("2.0"), BigFloat("0.5"), BigFloat("0.25")]
        result_sp2 = computeSP2(2, beta_bf)
        @test eltype(result_sp2) == BigFloat
        @test length(result_sp2) == 3

        # Test with BigFloat recurrence coefficients
        alpha_bf = BigFloat[BigFloat("0.0") for _ in 1:6]
        beta_full_bf = BigFloat[BigFloat("2.0"); [BigFloat(n^2) / (4 * BigFloat(n)^2 - 1) for n in 1:5]]
        op_bf = OrthoPoly("test_bigfloat", 5, alpha_bf, beta_full_bf, LegendreMeasure())
        @test eltype(op_bf.α) == BigFloat
        @test eltype(op_bf.β) == BigFloat

        # Test gauss with BigFloat coefficients
        nodes, weights = gauss(5, alpha_bf, beta_full_bf)
        @test eltype(nodes) == BigFloat
        @test eltype(weights) == BigFloat
    end

    @testset "MultiOrthoPoly BigFloat support" begin
        # Test multivariate evaluation with BigFloat
        ops = [LegendreOrthoPoly(3), HermiteOrthoPoly(3)]
        mop = MultiOrthoPoly(ops, 2)
        x_bf = BigFloat[BigFloat("0.1") BigFloat("0.2")]
        result = evaluate([1, 0], x_bf, mop)
        @test eltype(result) == BigFloat
    end

    @testset "Type preservation in mean/var/std" begin
        op = LegendreOrthoPoly(5)
        # mean with BigFloat should return BigFloat
        x_bf = BigFloat[BigFloat("1.0"), BigFloat("0.5"), BigFloat("0.1")]
        m = mean(x_bf, op)
        @test m isa BigFloat

        # var with BigFloat should return BigFloat
        v = var(x_bf, op)
        @test v isa BigFloat

        # std with BigFloat should return BigFloat
        s = std(x_bf, op)
        @test s isa BigFloat
    end

    @testset "assign2multi type preservation" begin
        x_bf = BigFloat[BigFloat("1.0"), BigFloat("0.5")]
        ind = calculateMultiIndices(2, 3)
        result = assign2multi(x_bf, 1, ind)
        @test eltype(result) == BigFloat
    end
end

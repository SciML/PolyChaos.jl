using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    points = [-0.5, 0.0, 0.5]
    multivariate_points = [-0.5 0.25; 0.0 -0.25; 0.5 0.75]
    pce_coefficients = collect(1.0:6.0)
    square = x -> x^2
    product = (x, y) -> x * y

    @compile_workload begin
        univariate_basis = LegendreOrthoPoly(4; Nrec = 8)
        multivariate_basis = MultiOrthoPoly(
            [LegendreOrthoPoly(2; Nrec = 4), HermiteOrthoPoly(2; Nrec = 4)], 2
        )

        evaluate(points, univariate_basis)
        integrate(square, univariate_basis)
        computeSP([1, 1], univariate_basis)
        computeTensorizedSP(3, univariate_basis)

        evaluate(multivariate_points, multivariate_basis)
        integrate(product, multivariate_basis)
        evaluatePCE(pce_coefficients, multivariate_points, multivariate_basis)
        mean(pce_coefficients, multivariate_basis)
        var(pce_coefficients, multivariate_basis)
        calculateMultiIndices(3, 2)
    end
end

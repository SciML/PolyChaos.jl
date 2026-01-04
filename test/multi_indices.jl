using Test, PolyChaos

A = [
    0 0 0 0;
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    2 0 0 0;
    1 1 0 0;
    1 0 1 0;
    1 0 0 1;
    0 2 0 0;
    0 1 1 0;
    0 1 0 1;
    0 0 2 0;
    0 0 1 1;
    0 0 0 2
]

@test calculateMultiIndices(4, 2) == A
@test calculateMultiIndices(3, 0) == zeros(Int64, 1, 3)
@test_throws DomainError calculateMultiIndices(3, -1)
@test_throws DomainError calculateMultiIndices(0, 1)
@test_throws DomainError calculateMultiIndices(-1, 1)

@test PolyChaos.numberPolynomials(3, 4) == factorial(3 + 4) ÷ (factorial(3) * factorial(4))
@test findUnivariateIndices(2, A) == [1, 3, 10]

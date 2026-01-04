export calculateMultiIndices,
    findUnivariateIndices

function calculateMultiIndices(d::Int, n::Int)
    # d denotes dimension of random variables/number of sources of uncertainty,
    # n the maximum degree of multivariate basis
    # function to calculate indices of multivariate basis following the algorithm
    # from "Spectral Methods for Uncertainty Quantiﬁcation; Le Maitre, Knio;2014" p.516-517
    # No::BigInt = factorial(BigInt(d+n))/(factorial(Bigint(d))*factorial(BigInt(n)));  #number of polynomials of multivariate basis
    n < 0 && throw(DomainError(n, "maximum degree must be non-negative"))
    d <= 0 && throw(DomainError(d, "number of uncertainties must be positive"))
    # catch case n == 0 --> No-d==0
    n == 0 && return zeros(Int64, 1, d)
    # non-pathological cases begin here
    No = numberPolynomials(d, n)
    inds = vcat(zeros(Int64, 1, d), Matrix(1I, d, d), zeros(Int64, No - d - 1, d))  #initiate index matrix for basis
    pi = ones(Int64, No, d)

    for k in 2:No
        g = 0
        for l in 1:d
            pi[k, l] = sum(pi[k - 1, :]) - g
            g = g + pi[k - 1, l]
        end
    end

    P = d + 1
    for k in 2:n
        L = P
        for j in 1:d, m in (L - pi[k, j] + 1):L

            P += 1
            inds[P, :] = inds[m, :]
            inds[P, j] = inds[P, j] + 1
        end
    end

    return inds
end

"""
    computes the number of polynomials with a multivariate basis
    `(d+n)!/(d!+n!)`
"""

function numberPolynomials(d::Int64, n::Int64)
    x, y = max(d, n), min(d, n)
    return UInt128(prod(UInt128(x + 1):UInt128(d + n)) ÷ factorial(UInt128(y)))
end

"""
    findUnivariateIndices(i::Int,ind::AbstractMatrix{Int64,2})

Given the multi-index `ind` this function returns all entries of the multivariate basis
that correspond to the `i`th univariate basis.
"""
function findUnivariateIndices(i::Int, ind::AbstractMatrix{Int})
    l, p = size(ind)
    i > p && throw(DomainError((i, p), "basis is $p-variate, you requested $i-variate"))
    deg = ind[end, end]
    deg < 0 && throw(DomainError(deg, "invalid degree"))
    col = ind[:, i]
    myind = zeros(Int64, deg)
    for deg_ in 1:deg
        myind[deg_] = findfirst(x -> x == deg_, col)
    end
    return pushfirst!(myind, 1)
end

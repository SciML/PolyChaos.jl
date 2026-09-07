export w_legendre,
    build_w_jacobi,
    w_jacobi,
    w_laguerre,
    w_hermite,
    build_w_genhermite,
    build_w_genlaguerre,
    w_meixner_pollaczek,
    build_w_meixner_pollaczek, w_gaussian,
    w_uniform01,
    w_uniform_11,
    w_logistic,
    w_genhermite,
    build_w_beta,
    build_w_gamma

_throwError(t) = throw(DomainError(t, "not in support"))

"""
    w_legendre(t)

Evaluate the unit Legendre weight on `[-1, 1]`.

# Arguments

- `t`: evaluation point in `[-1, 1]`.

Throws `DomainError` outside the support.
"""
function w_legendre(t)
    return -1.0 <= t <= 1.0 ? 1.0 : _throwError(t)
end

"""
    build_w_jacobi(a, b)

Return the Jacobi weight closure `t -> w_jacobi(t, a, b)`.

# Arguments

- `a`, `b`: exponents of `(1 - t)` and `(1 + t)`.
"""
function build_w_jacobi(a, b)
    return t -> w_jacobi(t, a, b)
end

"""
    w_jacobi(t, a, b)

Evaluate the unnormalized Jacobi weight `(1 - t)^a * (1 + t)^b` on
`[-1, 1]`.

# Arguments

- `t`: evaluation point.
- `a`, `b`: Jacobi shape parameters.
"""
function w_jacobi(t, a, b)
    return -1.0 <= t <= 1.0 ? (1 - t)^a * (1 + t)^b : _throwError(t)
end

"""
    w_hermite(t)

Evaluate the Hermite weight `exp(-t^2)`.

# Arguments

- `t`: real evaluation point.
"""
function w_hermite(t)
    return exp(-t^2)
end

"""
    build_w_genhermite(mu)

Return the generalized-Hermite weight closure.

# Arguments

- `mu`: exponent parameter used by `w_genhermite`.
"""
function build_w_genhermite(mu)
    return t -> w_genhermite(t, mu)
end

"""
    w_genhermite(t, mu)

Evaluate `abs(t)^(2mu) * exp(-t^2)`.

# Arguments

- `t`: real evaluation point.
- `mu`: generalized-Hermite parameter.
"""
function w_genhermite(t, μ)
    return abs(t)^(2 * μ) * exp(-t^2)
end

"""
    build_w_genlaguerre(a)

Return the generalized-Laguerre weight closure.

# Arguments

- `a`: exponent parameter used by the returned function.
"""
function build_w_genlaguerre(a)
    return t -> w_genlaguerre(t, a)
end

"""
    w_laguerre(t)

Evaluate the Laguerre weight `exp(-t)` on `[0, Inf)`.

# Arguments

- `t`: nonnegative evaluation point.
"""
function w_laguerre(t)
    return t >= 0.0 ? exp(-t) : _throwError(t)
end

"""
    w_meixner_pollaczek(t, lambda, phi)

Evaluate the Meixner-Pollaczek weight.

# Arguments

- `t`: real evaluation point.
- `lambda`: positive family parameter.
- `phi`: angle parameter.
"""
function w_meixner_pollaczek(t, lambda, phi)
    return 1 / (2pi) * exp((2 * phi - pi) * t) * abs(gamma(lambda + im * t))^2
end

"""
    build_w_meixner_pollaczek(lambda, phi)

Return a Meixner-Pollaczek weight closure.

# Arguments

- `lambda`: positive family parameter.
- `phi`: angle parameter.
"""
function build_w_meixner_pollaczek(lambda, phi)
    return t -> w_meixner_pollaczek(t, lambda, phi)
end

##################################################
# probability density functions
"""
    w_gaussian(t)

Evaluate the standard normal probability density at `t`.

# Arguments

- `t`: real evaluation point.

# Returns

The standard normal density evaluated at `t`.
"""
function w_gaussian(t)
    return 1 / (sqrt(2 * pi)) * exp(-0.5 * t^2)
end

"""
    build_w_beta(alpha, beta)

Return a beta-density closure on `(0, 1)`.

# Arguments

- `alpha`, `beta`: beta shape parameters.
"""
function build_w_beta(α, β)
    return t -> w_beta(t, α, β)
end

function w_beta(t, α, β)
    return -1 <= t <= 1 ? t^(α - 1) * (1 - t)^(β - 1) / beta(α, β) : _throwError(t)
end

"""
    build_w_gamma(alpha)

Return a unit-rate gamma-density closure.

# Arguments

- `alpha`: positive gamma shape parameter.
"""
function build_w_gamma(α)
    return t -> w_gamma(t, α)
end

function w_gamma(t, α)
    return t >= 0.0 ? (1 / gamma(α) * Float64(t)^(α - 1) * exp(-t)) : _throwError(t)
end

"""
    w_uniform01(t)

Evaluate the unit uniform density on `(0, 1)`.

# Arguments

- `t`: evaluation point in `(0, 1)`.

# Returns

The unit uniform density, or a `DomainError` outside the support.
"""
function w_uniform01(t)
    return 0.0 <= t <= 1.0 ? 1.0 : _throwError(t)
end

"""
    w_uniform_11(t)

Evaluate the uniform density `1 / 2` on `(-1, 1)`.

# Arguments

- `t`: evaluation point in `(-1, 1)`.

# Returns

The uniform density, or a `DomainError` outside the support.
"""
function w_uniform_11(t)
    return -1.0 <= t <= 1.0 ? 0.5 : _throwError(t)
end

"""
    w_logistic(t)

Evaluate the standard logistic probability density at `t`.

# Arguments

- `t`: real evaluation point.

# Returns

The standard logistic density evaluated at `t`.
"""
function w_logistic(t)
    return 0.25 * sech(0.5t)^2
end

export  w_legendre,
        build_w_jacobi,
        w_jacobi,
        w_laguerre,
        w_hermite,
        build_w_genhermite,
        build_w_genlaguerre,
        w_meixner_pollaczek,
        build_w_meixner_pollaczek,

        w_gaussian,
        w_uniform01,
        w_uniform_11,
        w_logistic,
        w_genhermite,
        build_w_beta,
        build_w_gamma

_throwError(t) = throw(DomainError(t, "not in support"))

function w_legendre(t)
    -1. <= t <= 1. ? 1. : _throwError(t)
end

function build_w_jacobi(a,b)
    return t->w_jacobi(t,a,b)
end

function w_jacobi(t,a,b)
    -1. <= t <= 1. ? (1-t)^a*(1+t)^b : _throwError(t)
end

function w_hermite(t)
    exp(-t^2)
end

function build_w_genhermite(mu)
    return t->w_genhermite(t,mu)
end

function w_genhermite(t,μ)
    abs(t)^(2*μ)*exp(-t^2)
end

function build_w_genlaguerre(a)
    return t -> w_genlaguerre(t,a)
end

function w_laguerre(t)
    t >= 0. ? exp(-t) : _throwError(t)
end

function w_meixner_pollaczek(t,lambda,phi)
    1 / (2pi) * exp((2*phi - pi)*t) * abs(gamma(lambda+im*t))^2
end

function build_w_meixner_pollaczek(lambda,phi)
    t->w_meixner_pollaczek(t,lambda,phi)
end

##################################################
# probability density functions
function w_gaussian(t)
    1 / (sqrt(2*pi))*exp(-0.5*t^2)
end

function build_w_beta(α,β)
    return t->w_beta(t,α,β)
end

function w_beta(t,α,β)
    -1 <= t <= 1 ? t^(α-1)*(1-t)^(β-1)/beta(α,β) : _throwError(t)
end

function build_w_gamma(α)
    return t->w_gamma(t,α)
end

function w_gamma(t,α)
    t >= 0. ? (1/gamma(α)*Float64(t)^(α-1)*exp(-t)) : _throwError(t)
end

function w_uniform01(t)
    0. <= t <= 1. ? 1. : _throwError(t)
end

function w_uniform_11(t)
    -1. <= t <= 1. ? 0.5 : _throwError(t)
end

function w_logistic(t)
    0.25*sech(0.5t)^2
end

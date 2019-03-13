export  w_legendre,
        build_w_jacobi,
        build_w_jacobi01,
        w_jacobi,
        w_jacobi01,
        w_laguerre,
        w_hermite,
        build_w_genhermite,
        build_w_genlaguerre,
        w_meixner_pollaczek,
        build_w_meixner_pollaczek,

        w_gaussian,
        w_uniform01,
        w_logistic,
        w_genhermite,
        build_w_beta,
        build_w_gamma

function w_legendre(t)
    -1. <= t <= 1. ? 1. : throw(error("$t not in support"))
end

function build_w_jacobi(a,b)
    @assert a > -1. && b > -1. "Invalid shape parameters"
    # return w(t) = (1-Float64(t))^a*(1+Float64(t))
    return t->w_jacobi(t,a,b)
end
function w_jacobi(t,a,b)
    -1. <=t<= 1. ? ((1-t)^a*(1+t)^b) : throw(error("$t not in support"))
end

function build_w_jacobi01(a,b)
    @assert a > -1. && b > -1. "Invalid shape parameters"
    # return w(t) = (1-Float64(t))^a*(1+Float64(t))
    return t->w_jacobi01(t,a,b)
end
function w_jacobi01(t,a,b)
    0. <=t<= 1. ? ((1-t)^a*t^b) : throw(error("$t not in support"))
end

function w_hermite(t)
    exp(-t^2)
end

function build_w_genhermite(mu)
    @assert mu>-0.5 "invalid parameter ($mu !> -0.5)"
    return t->w_genhermite(t,mu)
end
function w_genhermite(t,μ::Float64)
    abs(t)^(2*μ)*exp(-t^2)
end
function build_w_genlaguerre(a)
    @assert a>-1 "Invalid shape parameter"
    # return w(t) = Float64(t)^a*exp(-Float64(t))
    return t->w_genlaguerre(t,a)
end
function w_genlaguerre(t,a)
    t>=0. ? Float64(t)^a*exp(-Float64(t)) : throw(error("$t not in support"))
end

function w_laguerre(t)
    t>=0. ? exp(-t) : throw(error("$t not in support"))
end

function w_meixner_pollaczek(t,lambda,phi)
    1/(2*pi)*exp((2*phi-pi)*t)*abs(gamma(lambda+im*t))^2
end

function build_w_meixner_pollaczek(lambda,phi)
    @assert lambda>0 "lambda has to be positive"
    @assert 0<phi<pi "phi has to be between 0 and pi"
    t->w_meixner_pollaczek(t,lambda,phi)
end


##################################################
# probability density functions
function w_gaussian(t)
    1/(sqrt(2*pi))*exp(-0.5*t^2)
end

function build_w_beta(α,β)
    @assert α>-0. && β>-0. "Invalid shape parameters"
    return t->w_beta(t,α,β)
end
function w_beta(t,α,β)
    # 0. <=t<= 1. ? (t^(α-1)*(1-t)^(β-1)/beta(α,β)) : (error("$t not in support"))
    t^(α-1)*(1-t)^(β-1)/beta(α,β)
end



function build_w_gamma(α)
    @assert α>0 "Invalid shape parameter"
    return t->w_gamma(t,α)
end
function w_gamma(t,α)
    t>=0. ? (1/gamma(α)*Float64(t)^(α-1)*exp(-t)) : (error("$t not in support"))
end


function w_uniform01(t)
    0. <=t<= 1. ? 1. : (error("$t not in support"))
end


function w_logistic(t::Float64)
    0.25*sech(0.5t)^2
end

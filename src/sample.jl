export sampleInverseCDF

import LinearAlgebra: norm
import FFTW: fft

"""
    sampleInverseCDF(n::Int, pdf::Function, dom::Tuple{<:Real,<:Real})

Inverse transform sampling using Chebyshev technology.

Draw `n` samples from a probability density function `pdf` with support `dom`.
This method uses Chebyshev polynomial approximation of the PDF, constructs the
CDF via integration, and then uses bisection to sample from the inverse CDF.

Based on:
- Olver, Sheehan, and Alex Townsend. "Fast inverse transform sampling in one and two dimensions." arXiv preprint arXiv:1307.1223 (2013).
- https://github.com/dlfivefifty/InverseTransformSampling

This method works for any bounded PDF and does not require log-concavity
like adaptive rejection sampling.

# Arguments
- `n::Int`: Number of samples to draw
- `pdf::Function`: Probability density function (need not be normalized)
- `dom::Tuple{<:Real,<:Real}`: Domain (support) of the PDF as (lower, upper)

# Returns
- `Vector{Float64}`: Vector of `n` samples from the distribution

# Example
```julia
# Sample from a truncated normal
pdf = x -> exp(-x^2/2)
samples = sampleInverseCDF(1000, pdf, (-3.0, 3.0))
```
"""
function sampleInverseCDF(n::Int, pdf::Function, dom::Tuple{<:Real, <:Real})
    n < 1 && throw(DomainError(n, "Number of samples must be positive"))
    !(dom[1] < dom[2]) && throw(DomainError(dom, "Invalid domain bounds"))

    # Handle infinite domains by truncating
    a, b = _handle_domain(dom)

    # Scale PDF onto [-1, 1]
    map_fn(t) = (b - a) * (t + 1) / 2 + a
    f(x) = pdf(map_fn(x))

    # Construct Chebyshev approximation
    Y = _simple_constructor(f)
    integral = _sum_unit_interval(Y)
    if integral <= 0
        throw(ArgumentError("PDF integral is non-positive, check the PDF function"))
    end
    Y = Y ./ integral

    # Compute CDF coefficients
    cdf_coeffs = _simple_cumsum(Y)

    # Evaluate CDF at Chebyshev nodes to find effective support
    v = _simple_chebpolyval(cdf_coeffs)
    tol = 100 * eps()

    idx1 = something(findfirst(v_ -> v_ > tol, v), 1)
    idx2 = something(findlast(v_ -> v_ < 1 - tol, v), length(v))

    k = length(v) - 1
    chebnodes = sin.(π * (-k:2:k) / (2k))

    # Refined domain
    dnew_a = (b - a) * (chebnodes[idx1] + 1) / 2 + a
    dnew_b = (b - a) * (chebnodes[idx2] + 1) / 2 + a

    # Reconstruct on refined domain
    map_fn_new(t) = (dnew_b - dnew_a) * (t + 1) / 2 + dnew_a
    f_new(x) = pdf(map_fn_new(x))

    Y_new = _simple_constructor(f_new)
    integral_new = _sum_unit_interval(Y_new)
    if integral_new <= 0
        throw(ArgumentError("PDF integral is non-positive on refined domain"))
    end
    Y_new = Y_new ./ integral_new

    # Generate samples and map back
    samples_canonical = _generate_random_samples(Y_new, n)
    samples = map_fn_new.(samples_canonical)

    return samples
end

"""
Handle potentially infinite domains by truncating to a reasonable finite interval.

For unbounded distributions, we use a heuristic to find a good truncation point
that captures most of the probability mass while maintaining numerical stability.
"""
function _handle_domain(dom::Tuple{<:Real, <:Real})
    a, b = Float64(dom[1]), Float64(dom[2])
    # Use smaller truncation for better numerical stability
    # Most PDFs will have negligible mass outside [-20, 20]
    if isinf(a)
        a = -20.0
    end
    if isinf(b)
        b = 20.0
    end
    return (a, b)
end

"""
Generate random samples using inverse CDF via bisection.
"""
function _generate_random_samples(Y::Vector{Float64}, n::Int)
    # Compute CDF coefficients
    cdf = _simple_cumsum(Y)

    # Generate uniform random samples
    r = rand(n)

    # Bisection method
    a = -ones(n)
    b = ones(n)

    while norm(b - a, Inf) > eps()
        mid = (a + b) / 2
        vals = _clenshaw_evaluate(cdf, mid)
        I1 = (vals .- r) .<= -eps()
        I2 = (vals .- r) .>= eps()
        I3 = .~I1 .& .~I2
        a = I1 .* mid .+ I2 .* a .+ I3 .* mid
        b = I1 .* b .+ I2 .* mid .+ I3 .* mid
    end

    return (a + b) / 2
end

"""
Construct Chebyshev coefficients for a function on [-1, 1].
"""
function _simple_constructor(f::Function)
    Y = Float64[]

    for k_exp in 3:18
        k = 2^k_exp
        x = sin.(π * (-k:2:k) / (2k))
        vals = f.(x)

        # Laurent fold via FFT (DCT-I)
        Y = real.(fft([reverse(vals); vals[2:(end - 1)]]) / (2 * length(x) - 2))
        idx = findlast(y -> abs(y) > 10 * log2(k) * eps(), Y[1:k])

        if isnothing(idx)
            Y = [Y[1]]
            break
        elseif idx < k - 3
            Y = Y[1:idx]
            if length(Y) > 2
                Y[2:(end - 1)] .*= 2
            elseif length(Y) == 2
                # No middle elements
            end
            break
        end

        if k_exp == 18
            @warn "Chebyshev approximation did not converge"
        end
    end

    return Y
end

"""
Compute integral of Chebyshev series on [-1, 1].
"""
function _sum_unit_interval(c::Vector{Float64})
    n = length(c)
    n == 0 && return 0.0
    n == 1 && return 2 * c[1]

    # Set odd-index coefficients to zero (they integrate to zero)
    c_copy = copy(c)
    c_copy[2:2:end] .= 0

    # Build integration weights
    weights = zeros(n)
    weights[1] = 2.0
    for i in 3:n
        j = i - 1  # Chebyshev index (0-based)
        if iseven(j)
            weights[i] = 2 / (1 - j^2)
        end
    end

    return dot(weights, c_copy)
end

"""
Compute indefinite integral (cumsum) of Chebyshev series.
Returns CDF coefficients that give F(-1) = 0 and F(1) = 1.
"""
function _simple_cumsum(Y::Vector{Float64})
    Y_rev = reverse(Y)
    n = length(Y_rev)

    c = [0.0; 0.0; Y_rev]
    cout = zeros(n)

    # Compute C_{n+1} ... C_2 (coefficients for degrees n-1 down to 1)
    for i in 1:(n - 1)
        cout[i] = (c[i + 2] - c[i]) / (2 * (n - i + 1))
    end

    # Compute C_1
    cout[n] = c[end] - c[end - 2] / 2

    # Compute C_0 (constant to make F(-1) = 0)
    v = ones(n)
    v[(end - 1):-2:1] .= -1
    c0 = dot(cout, v)
    push!(cout, c0)

    return cout
end

"""
Evaluate Chebyshev series at points using Clenshaw algorithm.
"""
function _clenshaw_evaluate(c::Vector{Float64}, x::Vector{Float64})
    n = length(c)
    n == 0 && return zeros(length(x))

    bk1 = zeros(length(x))
    bk2 = zeros(length(x))
    x2 = 2 .* x

    for k in 1:(n - 1)
        bk = c[k] .+ x2 .* bk1 .- bk2
        bk2 = bk1
        bk1 = bk
    end

    return c[n] .+ 0.5 .* x2 .* bk1 .- bk2
end

"""
Convert Chebyshev coefficients to values at Chebyshev nodes.
"""
function _simple_chebpolyval(c::Vector{Float64})
    c = copy(c)
    lc = length(c)
    lc == 0 && return Float64[]
    lc == 1 && return c

    ii = 2:(lc - 1)
    c[ii] .= 0.5 .* c[ii]

    v = [reverse(c); c[ii]]
    v = real.(ifft(v))

    result = (lc - 1) .* [2 * v[1]; v[ii] .+ v[2 * lc .- ii]; 2 * v[lc]]
    result = reverse(result)

    return result
end

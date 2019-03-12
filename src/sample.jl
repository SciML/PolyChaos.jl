export sample_inversecdf
"""
```
sample_inversecdf(pdf::Function, dom::Tuple{Float64,Float64}, N::Int64)
```
Inverse transform sampling using Chebyshev technology taken from
Olver, Sheehan, and Alex Townsend. "Fast inverse transform sampling in one and two dimensions." arXiv preprint arXiv:1307.1223 (2013).
https://github.com/dlfivefifty/InverseTransformSampling/blob/master/sample.m

samples from a given probabiliy density function.
"""
function sample_inversecdf(pdf::Function, dom::Tuple{Float64,Float64}, N::Int64)::Vector{Float64}
    # scale pdf onto [-1,1]
    f(t) = pdf(scale_pdf(dom, t))

    Y = simple_constructor(f)
    Y = Y ./ sum_unit_interval(Y)

    #cumulative density function
    cout = simple_cumsum(Y)
    cdf = cout

    v = simple_chebpolyval(cdf)
    tol = 100*eps()
    idx1 = findfirst(v_ -> v_ > tol, v)
    idx2 = findlast(v_ -> v_ < 1-tol, v)

    k = length(v)-1
    x = map(n->sin(0.5*pi*n/k),-k:2:k)
    dnew = Tuple(scale_pdf(dom, [x[idx1], x[idx2]]))
    # scale pdf onto [-1,1];
    f_(t) = pdf(scale_pdf(dnew, t))

    Y = simple_constructor(f_)
    Y = Y ./ sum_unit_interval(Y)

    scale_pdf(dnew, (generate_random_samples(Y, N)))
end

function sample_inversecdf(pdf::Function, dom::Tuple{T,TT},N::Int64) where T <: Number where TT <: Number
    sample_inversecdf(pdf,Float64.(dom),N)
end

function scale_pdf(dom::Tuple{Float64, Float64}, t)
    ( ( dom[2] - dom[1] ) .* (t .+ 1) ) ./2 .+ dom[1]
end

function generate_random_samples(Y::Vector{Float64}, N::Int64)::Vector{Float64}
    # cumulative density function
    c = simple_cumsum(Y)
    # generate random samples from uniform distribution
    r = rand(N)
    # bisection method
    a, b = -ones(N), ones(N)
    while norm(b-a, Inf) > 1e-14
        vals = Clenshaw_evaluate(c,(a+b)/2)
        I1 = (vals - r) .<= -1e-14
        I2 = (vals - r) .>= 1e-14
        I3 = .~I1 .& .~I2
        a = I1 .* (a + b) / 2 + I2 .* a + I3 .* (a + b) / 2
        b = I1 .* b + I2 .* (a + b) / 2 + I3 .* (a + b) / 2
    end
    (a + b) / 2
end

function simple_cumsum(Y_::Vector{Float64})
    Y = reverse(Y_)
    n = length(Y)
    c = vcat(zeros(Float64,2), Y)
    cout = zeros(Float64,n)   # initialize vector C_r
    @inbounds cout[1:n-1] = ( c[3:end-1] - c[1:end-3]) ./ (2 * (n:-1:2) )
    @inbounds cout[n] = c[end] - c[end-2] / 2
    v = ones(n)
    @inbounds v[end-1:-2:1] = -ones(trunc(Int, n/2))
    push!(cout, dot(cout,v))
end

function simple_constructor(f::Function)::Vector{Float64}
    Y = Vector{Float64}()
    # simple Chebfun constructor.
    for k = [2^i for i=3:18]
        x = map(n->sin(0.5*pi*n/k),-k:2:k)
        vals = f.(x)
        #Laurent fold in columns.
        Y = real(fft([vals[end:-1:1, :]; vals[2:end-1,:]])/(2*length(x)-2))
        idx = findlast(y->y > 10*log2(k)*eps(), abs.(Y[1:k]))
        if idx < k - 3
            Y = Y[1:idx]
            Y[2:end-1,:] = 2*Y[2:end-1,:]
            break
        end
        k == 2^18 && error("Sampling from inverse CDF failed.")
    end
    return Y
end

function Clenshaw_evaluate(c::Vector{Float64},x::Vector{Float64})::Array{Float64}
    bk1, bk2 = zero(x), zero(x)
    x *= 2
    for k = 1:length(c) - 1
        @inbounds bk = c[k] .+ x .* bk1 .- bk2
        bk2 = bk1
        bk1 = bk
    end
    @inbounds vals = last(c) .+ .5 * x .* bk1 .- bk2
end

function sum_unit_interval(c::Vector{Float64})::Float64
    n = length(c)
    n == 1 && return 2c[1]
    @inbounds  c[2:2:end] .= 0
    @inbounds return dot([2; 0; 2 ./ (1 .- ((2:n-1)).^2) ], c)
end

function simple_chebpolyval(c::Vector{Float64})
    lc = length(c)
    lc == 1 && return c

    ii = 2:lc-1
    c[ii] *= 0.5
    v = [ c[end:-1:1] ; c[ii] ]
    if isreal(c)
        v = real(ifft(v))
    elseif isreal(1i*c)
        v = 1i*real(ifft(imag(v)))
    else
        v = ifft(v)
    end
    v = (lc-1)*[ 2 * v[1] ; v[ii] + v[2 * lc .- ii] ; 2 * v[lc]]
    reverse(v)
end

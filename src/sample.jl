"""
```
sample_withoutOOP(pdf::Function, dom::Tuple{Float64,Float64}, N::Int64)
```
Inverse transform sampling using Chebyshev technology without OOP
Olver, Sheehan, and Alex Townsend. "Fast inverse transform sampling in one and two dimensions." arXiv preprint arXiv:1307.1223 (2013).
https://github.com/dlfivefifty/InverseTransformSampling/blob/master/sample.m

samples from a given probabiliy density function.
"""
function sample_withoutOOP(pdf::Function, dom::Tuple{Float64,Float64}, N::Int64)::Vector{Float64}
    # scale pdf onto [-1,1]
    map(t) = (dom[2]-dom[1]).*(t+1)./2 + dom[1]
    f(x) = pdf(map(x))

    Y = simple_constructor(f)
    out = sum_unit_interval(Y)
    Y = Y./out

    #cumulative density function
    cout = simple_cumsum(Y)
    cdf = cout

    v = simple_chebpolyval(cdf)
    tol = 100*eps()

    idx1 = findfirst(v_ -> v_ > tol, v) 
    idx2 = findlast(v_ -> v_ < 1-tol, v) 

    k = length(v)-1
    x =  transpose(sin.(pi*(-k:2:k)/(2*k)))
    
    dnew = (dom[2] - dom[1]).*([x[idx1] x[idx2]] .+ 1)./2 .+ dom[1]

    #scale pdf onto [-1,1];
    map_(t) = (dnew[2]-dnew[1]).*(t.+1)./2 .+ dnew[1]
    f_(x) = pdf(map_(x))

    Y = simple_constructor(f_)
    out = sum_unit_interval(Y)
    Y = Y./out

    x = map_(generate_random_samples( Y, N ) )

    return x
end


function generate_random_samples(Y::Vector{Float64}, N::Int64)::Vector{Float64}
    #cumulative density function
    cout = simple_cumsum(Y)
    cdf = cout

    # generate random samples from uniform distribution
    r = rand(N)
    c = cdf

    # bisection method
    a = -ones(N)
    b =  ones(N)

    while norm(b-a, Inf) > eps()
        vals = Clenshaw_evaluate(c,(a+b)/2)
        I1 = ((vals-r) .<= -eps())
        I2 = ((vals-r) .>= eps())
        I3 = .~I1 .& .~I2
        a = I1.*(a+b)/2 + I2.*a + I3.*(a+b)/2
        b = I1.*b + I2.*(a+b)/2 + I3.*(a+b)/2
    end
    x = (a+b)/2

    return x
end


function simple_cumsum(Y::Vector{Float64})
    Y = Y[end:-1:1]
    n = length(Y)

    c = [0;0;Y]       # obtain Cheb coeffs c_r
    cout = zeros(n)   # initialize vector C_r

    # compute C_(n+1) ... C_2
    cout[1:n-1] = (c[3:end-1]-c[1:end-3])./(2*(n:-1:2))

    # compute C_1
    cout[n,1] = c[end] - c[end-2]/2;         
    v = ones(n)
    v[end-1:-2:1] = -ones(trunc(Int, n/2))

    # compute C_0
    push!(cout, cout'*v)

    return cout
end

function simple_constructor(f::Function)::Array{Float64,1}
    Y = []

    # simple Chebfun constructor.
    for k = [2^i for i=3:18]
        x = sin.(π*(-k:2:k)/(2*k))
        vals = f.(x)

        #Laurent fold in columns.
        Y = real(fft([vals[end:-1:1, :]; vals[2:end-1,:]])/(2*length(x)-2))
        idx = findlast(y -> y > 10*log2(k)*eps(), Y[1:k]) 

        if idx < k-3
            Y = Y[1:idx]
            Y[2:end-1,:] = 2*Y[2:end-1,:]
            break
        end

        k == 2^18 ? error() : ()
    end

    return Y
end



function Clenshaw_evaluate(c,x)::Array{Float64}
    bk1 = zeros(size(x))
    bk2 = bk1
    x = 2*x

    for k = 1:size(c,1)-1
        bk = c[k] .+ x.*bk1 .- bk2
        bk2 = bk1
        bk1 = bk
    end

    vals = c[end] .+ .5*x.*bk1 .- bk2

    return vals
end


function sum_unit_interval(c)::Float64
    #Integral on the unit interval
    n = length(c)
    if n == 1 return c*2 end

    c[2:2:end] .= 0
    
    return [2 0 2/(1 .- ((2:n-1)).^2)]*c
end

function simple_chebpolyval(c::Array{Float64,1})::Array{Float64}
    # Convert coefficients to values. 
    c = c[:]       # Input should be a column vector
    lc = length(c)

    if lc == 1 return c end

    ii = 2:lc-1
    c[ii] = 0.5*c[ii]
    v = [c[end:-1:1] ; c[ii]]
    if isreal(c)
        v = real(ifft(v))
    elseif isreal(1i*c)
        v = 1i*real(ifft(imag(v)))
    else
        v = ifft(v)
    end

    v = (lc-1)*[ 2*v[1] ; v[ii]+v[2*lc.-ii] ; 2*v[lc]]
    v = v[end:-1:1]

    return v
end

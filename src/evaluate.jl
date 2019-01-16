export evaluate

"""
__Univariate__
```
evaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})
evaluate(n::Int64,x::Float64,a::Vector{Float64},b::Vector{Float64})
evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly)
evaluate(n::Int64,x::Float64,op::OrthoPoly)
```
Evaluate the `n`-th univariate basis polynomial at point(s) `x`
The function is multiply dispatched to facilitate its use with the composite type `OrthoPoly`

!!! note
    - `n` is the degree of the univariate basis polynomial
    - `length(x) = N`, where `N` is the number of points
    - `(a,b)` are the recursion coefficients

__Multivariate__
```
evaluate(n::Vector{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
evaluate(n::Vector{Int64},x::Vector{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
evaluate(n::Vector{Int64},x::Matrix{Float64},op::MultiOrthoPoly)
evaluate(n::Vector{Int64},x::Vector{Float64},op::MultiOrthoPoly)
```
Evaluate the n-th p-variate basis polynomial at point(s) x
The function is multiply dispatched to facilitate its use with the composite type `MultiOrthoPoly`

!!! note
    - `n is a multi-index
    - `length(n) == p`, i.e. a p-variate basis polynomial
    - `size(x) = (N,p)`, where `N` is the number of points
    - `size(a)==size(b)=p`.
"""
function evaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})
    @assert n>=0 "Degree n has to be non-negative (currently n=$n)."
    # if length(a)==0 warn("Length of a is 0.") end
    @assert length(a)==length(b) "Inconsistent number of recurrence coefficients."
    @assert n<=length(a) "Specified degree is $n, but you only provided $(length(a)) coefficients."
    # recurrence relation for orthogonal polynomials
    nx = length(x)
    pminus, p = zeros(nx), ones(nx)
    if n==0
        if nx==1
            return p[1]
        else
            return p
        end
    end
    pplus = (x .- a[1]).*p .- b[1]*pminus
    for k=2:n
        pminus = p
        p = pplus
        pplus = (x .- a[k]).*p .- b[k]*pminus
    end
    nx==1 ? pplus[1] : pplus
end
evaluate(n::Int64,x::Float64,a::Vector{Float64},b::Vector{Float64}) = evaluate(n,[x],a,b)
evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly) = evaluate(n,x,op.α,op.β)
evaluate(n::Int64,x::Float64,op::OrthoPoly) = evaluate(n,[x],op)

function evaluate(n::Vector{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
    @assert length(a)==length(b)
    @assert length(n)==size(x,2)
    Nmulti = length(n)
    Npoints = size(x,1)
    val = ones(Float64,Npoints)
    for i=1:length(n)
        val = val.*evaluate(n[i],x[:,i],a[i],b[i])
    end
    return val
end

function evaluate(n::Vector{Int64},x::Vector{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
    x_ = zeros(1,length(x))
    x_[:] = x[:]
    evaluate(n,x_,a,b)
end

function evaluate(n::Vector{Int64},x::Matrix{Float64},op::MultiOrthoPoly)
    p = length(n)
    @assert p==length(op.uni)
    a,b = [ op.uni[i].α for i=1:p], [ op.uni[i].β for i=1:p]
    evaluate(n,x,a,b)
end

function evaluate(n::Vector{Int64},x::Vector{Float64},op::MultiOrthoPoly)
    x_ = zeros(1,length(x))
    x_[:] = x[:]
    evaluate(n,x_,op)
end

export evaluate

"""
__Univariate__
```
evaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})
evaluate(n::Int64,x::Float64,a::Vector{Float64},b::Vector{Float64})
evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly)
evaluate(n::Int64,x::Float64,op::OrthoPoly)
evaluate(n::Int64,x::Vector{Float64},opq::OrthoPolyQ) = evaluate(n,x,opq.op)
evaluate(n::Int64,x::Float64,opq::OrthoPolyQ) = evaluate(n,[x],opq.op)
```
Evaluate the `n`-th univariate basis polynomial at point(s) `x`
The function is multiply dispatched to facilitate its use with the composite type `OrthoPoly`

If several basis polynomials (stored in `ns`) are to be evaluated at points `x`, then call
```
evaluate(ns::Vector{Int64},x::Vector{Float64},op::OrthoPoly) = evaluate(ns,x,op.α,op.β)
evaluate(ns::Vector{Int64},x::Float64,op::OrthoPoly) = evaluate(ns,[x],op)
evaluate(ns::Vector{Int64},x::Vector{Float64},opq::OrthoPolyQ) = evaluate(ns,x,opq.op)
evaluate(ns::Vector{Int64},x::Float64,opq::OrthoPolyQ) = evaluate(ns,[x],opq.op)
```

If *all* basis polynomials are to be evaluated at points `x`, then call
```
evaluate(x::Vector{Float64},op::OrthoPoly) = evaluate(collect(0:op.deg),x,op)
evaluate(x::Float64,op::OrthoPoly) = evaluate([x],op)
evaluate(x::Vector{Float64},opq::OrthoPolyQ) = evaluate(x,opq.op)
evaluate(x::Float64,opq::OrthoPolyQ) = evaluate([x],opq)
```
which returns an Array of dimensions `(length(x),op.deg+1)`.

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

If several basis polynomials are to be evaluated at points `x`, then call
```
evaluate(ind::Matrix{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
evaluate(ind::Matrix{Int64},x::Matrix{Float64},op::MultiOrthoPoly)
```
where `ind` is a matrix of multi-indices.

If *all* basis polynomials are to be evaluated at points `x`, then call
```
evaluate(x::Matrix{Float64},mop::MultiOrthoPoly) = evaluate(mop.ind,x,mop)
```
which returns an array of dimensions `(mop.dim,size(x,1))`.

!!! note
    - `n` is a multi-index
    - `length(n) == p`, i.e. a p-variate basis polynomial
    - `size(x) = (N,p)`, where `N` is the number of points
    - `size(a)==size(b)=p`.
"""
function evaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})
    @assert n >= 0 "Degree n has to be non-negative (currently n=$n)."
    # if length(a)==0 warn("Length of a is 0.") end
    @assert length(a) == length(b) "Inconsistent number of recurrence coefficients."
    @assert n <= length(a) "Specified degree is $n, but you only provided $(length(a)) coefficients."
    # recurrence relation for orthogonal polynomials
    nx = length(x)
    pminus, p = zeros(nx), ones(nx)
    if n==0
        if nx==1
            return first(p)
        else
            return p
        end
    end
    pplus = (x .- first(a)).*p .- first(b)*pminus
    for k in 2:n
        pminus = p
        p = pplus
        @inbounds pplus = (x .- a[k]).*p .- b[k]*pminus
    end
    nx == 1 ? first(pplus) : pplus
end
evaluate(n::Int64,x::Float64,a::Vector{Float64},b::Vector{Float64}) = evaluate(n,[x],a,b)
evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly) = evaluate(n,x,op.α,op.β)
evaluate(n::Int64,x::Float64,op::OrthoPoly) = evaluate(n,[x],op)
evaluate(n::Int64,x::Vector{Float64},opq::OrthoPolyQ) = evaluate(n,x,opq.op)
evaluate(n::Int64,x::Float64,opq::OrthoPolyQ) = evaluate(n,[x],opq.op)

# univariate + several bases
function evaluate(ns::Vector{Int64},x::Array{Float64},a::Vector{Float64},b::Vector{Float64})
    hcat(map(i->evaluate(i,x,a,b),ns)...)
end
evaluate(ns::Vector{Int64},x::Float64,a::Vector{Float64},b::Vector{Float64}) = evaluate(ns,[x],a,b)

evaluate(ns::Vector{Int64},x::Vector{Float64},op::OrthoPoly) = evaluate(ns,x,op.α,op.β)
evaluate(ns::Vector{Int64},x::Float64,op::OrthoPoly) = evaluate(ns,[x],op)
evaluate(x::Vector{Float64},op::OrthoPoly) = evaluate(collect(0:op.deg),x,op)
evaluate(x::Float64,op::OrthoPoly) = evaluate([x],op)

evaluate(ns::Vector{Int64},x::Vector{Float64},opq::OrthoPolyQ) = evaluate(ns,x,opq.op)
evaluate(ns::Vector{Int64},x::Float64,opq::OrthoPolyQ) = evaluate(ns,[x],opq.op)
evaluate(x::Vector{Float64},opq::OrthoPolyQ) = evaluate(x,opq.op)
evaluate(x::Float64,opq::OrthoPolyQ) = evaluate([x],opq)

# multivariate
function evaluate(n::Vector{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
    @assert length(n) == size(x,2) "number of univariate bases (= $(length(n))) inconsistent with columns points x (= $(size(x,2)))"
    val = ones(Float64,size(x,1))
    for i in 1:length(n)
        @inbounds val = val.*evaluate(n[i],x[:,i],a[i],b[i])
    end
    return val
end
evaluate(n::Vector{Int64},x::Vector{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}}) = evaluate(n,reshape(x,1,length(x)),a,b)
evaluate(n::Vector{Int64},x::Matrix{Float64},op::MultiOrthoPoly) = evaluate(n,x,coeffs(op)...)
evaluate(n::Vector{Int64},x::Vector{Float64},op::MultiOrthoPoly) = evaluate(n,reshape(x,1,length(x)),op)

# using multi-index + multivariate
function evaluate(ind::Matrix{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})
    vals = map(i->evaluate(ind[i,:],x,a,b),Base.OneTo(size(ind,1)))
    hcat(vals...) |> transpose |> Matrix
end
    evaluate(ind::Matrix{Int64},x::Matrix{Float64},op::MultiOrthoPoly) = evaluate(ind,x,coeffs(op)...)
evaluate(x::Matrix{Float64},mop::MultiOrthoPoly) = evaluate(mop.ind,x,mop)

evaluate(ind::Matrix{Int64},x::Vector{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}}) = evaluate(ind,reshape(x,1,length(x)),a,b)
evaluate(ind::Matrix{Int64},x::Vector{Float64},op::MultiOrthoPoly) = evaluate(ind,reshape(x,1,length(x)),coeffs(op)...)
evaluate(x::Vector{Float64},mop::MultiOrthoPoly) = evaluate(mop.ind,reshape(x,1,length(x)),mop)

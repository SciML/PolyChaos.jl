export evaluate

"""
__Univariate__
```
evaluate(n::Int,x::Array{<:Real},a::AbstractVector{<:Real},b::AbstractVector{<:Real})
evaluate(n::Int,x::Real,a::AbstractVector{<:Real},b::AbstractVector{<:Real})
evaluate(n::Int,x::AbstractVector{<:Real},op::AbstractOrthoPoly)
evaluate(n::Int,x::Real,op::AbstractOrthoPoly)
```
Evaluate the `n`-th univariate basis polynomial at point(s) `x`
The function is multiply dispatched to facilitate its use with the composite type `AbstractOrthoPoly`

If several basis polynomials (stored in `ns`) are to be evaluated at points `x`, then call
```
evaluate(ns::AbstractVector{<:Int},x::AbstractVector{<:Real},op::AbstractOrthoPoly) = evaluate(ns,x,op.α,op.β)
evaluate(ns::AbstractVector{<:Int},x::Real,op::AbstractOrthoPoly) = evaluate(ns,[x],op)
```

If *all* basis polynomials are to be evaluated at points `x`, then call
```
evaluate(x::AbstractVector{<:Real},op::AbstractOrthoPoly) = evaluate(collect(0:op.deg),x,op)
evaluate(x::Real,op::AbstractOrthoPoly) = evaluate([x],op)
```
which returns an Array of dimensions `(length(x),op.deg+1)`.

!!! note
    - `n` is the degree of the univariate basis polynomial
    - `length(x) = N`, where `N` is the number of points
    - `(a,b)` are the recursion coefficients

__Multivariate__
```
evaluate(n::AbstractVector{<:Int},x::AbstractMatrix{<:Real},a::Vector{<:AbstractVector{<:Real}},b::Vector{<:AbstractVector{<:Real}})
evaluate(n::AbstractVector{<:Int},x::AbstractVector{<:Real},a::Vector{<:AbstractVector{<:Real}},b::Vector{<:AbstractVector{<:Real}})
evaluate(n::AbstractVector{<:Int},x::AbstractMatrix{<:Real},op::MultiOrthoPoly)
evaluate(n::AbstractVector{<:Int},x::AbstractVector{<:Real},op::MultiOrthoPoly)
```
Evaluate the n-th p-variate basis polynomial at point(s) x
The function is multiply dispatched to facilitate its use with the composite type `MultiOrthoPoly`

If several basis polynomials are to be evaluated at points `x`, then call
```
evaluate(ind::AbstractMatrix{<:Int},x::AbstractMatrix{<:Real},a::Vector{<:AbstractVector{<:Real}},b::Vector{<:AbstractVector{<:Real}})
evaluate(ind::AbstractMatrix{<:Int},x::AbstractMatrix{<:Real},op::MultiOrthoPoly)
```
where `ind` is a matrix of multi-indices.

If *all* basis polynomials are to be evaluated at points `x`, then call
```
evaluate(x::AbstractMatrix{<:Real},mop::MultiOrthoPoly) = evaluate(mop.ind,x,mop)
```
which returns an array of dimensions `(mop.dim,size(x,1))`.

!!! note
    - `n` is a multi-index
    - `length(n) == p`, i.e. a p-variate basis polynomial
    - `size(x) = (N,p)`, where `N` is the number of points
    - `size(a)==size(b)=p`.
"""
function evaluate(n::Int, x::AbstractArray{<:Real}, a::AbstractVector{<:Real},
                  b::AbstractVector{<:Real})
    @assert n>=0 "Degree n has to be non-negative (currently n=$n)."
    # if length(a)==0 warn("Length of a is 0.") end
    @assert length(a)==length(b) "Inconsistent number of recurrence coefficients."
    @assert n<=length(a) "Specified degree is $n, but you only provided $(length(a)) coefficients."
    # recurrence relation for orthogonal polynomials
    nx = length(x)
    pminus, p = zeros(nx), ones(nx)
    if n == 0
        if nx == 1
            return first(p)
        else
            return p
        end
    end
    pplus = (x .- first(a)) .* p .- first(b) * pminus
    for k in 2:n
        pminus = p
        p = pplus
        @inbounds pplus = (x .- a[k]) .* p .- b[k] * pminus
    end
    nx == 1 ? first(pplus) : pplus
end
function evaluate(n::Int, x::Real, a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    evaluate(n, [x], a, b)
end
function evaluate(n::Int, x::AbstractVector{<:Real}, op::AbstractOrthoPoly)
    evaluate(n, x, op.α, op.β)
end
evaluate(n::Int, x::Real, op::AbstractOrthoPoly) = evaluate(n, [x], op)

# univariate + several bases
function evaluate(ns, x::AbstractArray{<:Real}, a::AbstractVector{<:Real},
                  b::AbstractVector{<:Real})
    hcat(map(i -> evaluate(i, x, a, b), ns)...)
end
function evaluate(ns, x::Real, a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    evaluate(ns, [x], a, b)
end

evaluate(ns, x::AbstractVector{<:Real}, op::AbstractOrthoPoly) = evaluate(ns, x, op.α, op.β)
evaluate(ns, x::Real, op::AbstractOrthoPoly) = evaluate(ns, [x], op)
function evaluate(x::AbstractVector{<:Real}, op::AbstractOrthoPoly)
    evaluate(collect(0:(op.deg)), x, op)
end
evaluate(x::Real, op::AbstractOrthoPoly) = evaluate([x], op)

# multivariate
function evaluate(n::AbstractVector{<:Int}, x::AbstractMatrix{<:Real},
                  a::AbstractVector{<:AbstractVector{<:Real}},
                  b::AbstractVector{<:AbstractVector{<:Real}})
    @assert length(n)==size(x, 2) "number of univariate bases (= $(length(n))) inconsistent with columns points x (= $(size(x,2)))"
    val = ones(Float64, size(x, 1))
    for i in 1:length(n)
        @inbounds val = val .* evaluate(n[i], x[:, i], a[i], b[i])
    end
    return val
end
function evaluate(n::AbstractVector{<:Int}, x::AbstractVector{<:Real},
                  a::AbstractVector{<:AbstractVector{<:Real}},
                  b::AbstractVector{<:AbstractVector{<:Real}})
    evaluate(n, reshape(x, 1, length(x)), a, b)
end
function evaluate(n::AbstractVector{<:Int}, x::AbstractMatrix{<:Real}, op::MultiOrthoPoly)
    evaluate(n, x, coeffs(op)...)
end
function evaluate(n::AbstractVector{<:Int}, x::AbstractVector{<:Real}, op::MultiOrthoPoly)
    evaluate(n, reshape(x, 1, length(x)), op)
end

# using multi-index + multivariate
function evaluate(ind::AbstractMatrix{<:Int}, x::AbstractMatrix{<:Real},
                  a::AbstractVector{<:AbstractVector{<:Real}},
                  b::AbstractVector{<:AbstractVector{<:Real}})
    vals = map(i -> evaluate(ind[i, :], x, a, b), Base.OneTo(size(ind, 1)))
    hcat(vals...) |> transpose |> Matrix
end

function evaluate(ind::AbstractMatrix{<:Int}, x::AbstractMatrix{<:Real}, op::MultiOrthoPoly)
    evaluate(ind, x, coeffs(op)...)
end
evaluate(x::AbstractMatrix{<:Real}, mop::MultiOrthoPoly) = evaluate(mop.ind, x, mop)

function evaluate(ind::AbstractMatrix{<:Int}, x::AbstractVector{<:Real},
                  a::AbstractVector{<:AbstractVector{<:Real}},
                  b::AbstractVector{<:AbstractVector{<:Real}})
    evaluate(ind, reshape(x, 1, length(x)), a, b)
end
function evaluate(ind::AbstractMatrix{<:Int}, x::AbstractVector{<:Real}, op::MultiOrthoPoly)
    evaluate(ind, reshape(x, 1, length(x)), coeffs(op)...)
end
function evaluate(x::AbstractVector{<:Real}, mop::MultiOrthoPoly)
    evaluate(mop.ind, reshape(x, 1, length(x)), mop)
end

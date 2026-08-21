export evaluate

"""
    evaluate(n, x, α, β)
    evaluate(n, x, op::AbstractOrthoPoly)
    evaluate(x, op::AbstractOrthoPoly)
    evaluate(n, x, op::MultiOrthoPoly)
    evaluate(x, mop::MultiOrthoPoly)

Evaluate monic orthogonal basis polynomials using their three-term recurrence.
Univariate points can be scalars or vectors. For a `MultiOrthoPoly`, each row
of a point matrix is one multivariate point and each row of an index matrix is
one multi-index.

# Arguments

- `n`: polynomial degree for a univariate basis, or a multi-index for one
  multivariate basis polynomial.
- `ns`: degrees of several univariate basis polynomials.
- `ind`: matrix whose rows are the multi-indices of several multivariate basis
  polynomials.
- `x`: evaluation point(s). A multivariate point matrix has one point per row
  and one variable per column.
- `α`, `β`: recurrence coefficient vectors, or one vector per variable.
- `op`, `mop`: univariate or multivariate orthogonal-polynomial basis.

# Returns

Values of the requested basis polynomials. Evaluating all polynomials of `op`
returns an array with one row per point and one column per degree.

# Examples

```jldoctest
julia> using PolyChaos

julia> evaluate(1, 0.5, LegendreOrthoPoly(2))
0.5
```
"""
function evaluate(
        n::Int, x::AbstractArray{<:Real}, a::AbstractVector{<:Real},
        b::AbstractVector{<:Real}
    )
    @assert n >= 0 "Degree n has to be non-negative (currently n=$n)."
    # if length(a)==0 warn("Length of a is 0.") end
    @assert length(a) == length(b) "Inconsistent number of recurrence coefficients."
    @assert n <= length(a) "Specified degree is $n, but you only provided $(length(a)) coefficients."
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
    return nx == 1 ? first(pplus) : pplus
end
function evaluate(n::Int, x::Real, a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    return evaluate(n, [x], a, b)
end
function evaluate(n::Int, x::AbstractVector{<:Real}, op::AbstractOrthoPoly)
    return evaluate(n, x, op.α, op.β)
end
evaluate(n::Int, x::Real, op::AbstractOrthoPoly) = evaluate(n, [x], op)

# univariate + several bases
function evaluate(
        ns, x::AbstractArray{<:Real}, a::AbstractVector{<:Real},
        b::AbstractVector{<:Real}
    )
    return hcat(map(i -> evaluate(i, x, a, b), ns)...)
end
function evaluate(ns, x::Real, a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    return evaluate(ns, [x], a, b)
end

evaluate(ns, x::AbstractVector{<:Real}, op::AbstractOrthoPoly) = evaluate(ns, x, op.α, op.β)
evaluate(ns, x::Real, op::AbstractOrthoPoly) = evaluate(ns, [x], op)
function evaluate(x::AbstractVector{<:Real}, op::AbstractOrthoPoly)
    return evaluate(collect(0:(op.deg)), x, op)
end
evaluate(x::Real, op::AbstractOrthoPoly) = evaluate([x], op)

# multivariate
function evaluate(
        n::AbstractVector{<:Int}, x::AbstractMatrix{<:Real},
        a::AbstractVector{<:AbstractVector{<:Real}},
        b::AbstractVector{<:AbstractVector{<:Real}}
    )
    @assert length(n) == size(x, 2) "number of univariate bases (= $(length(n))) inconsistent with columns points x (= $(size(x, 2)))"
    val = ones(Float64, size(x, 1))
    for i in 1:length(n)
        @inbounds val = val .* evaluate(n[i], x[:, i], a[i], b[i])
    end
    return val
end
function evaluate(
        n::AbstractVector{<:Int}, x::AbstractVector{<:Real},
        a::AbstractVector{<:AbstractVector{<:Real}},
        b::AbstractVector{<:AbstractVector{<:Real}}
    )
    return evaluate(n, reshape(x, 1, length(x)), a, b)
end
function evaluate(n::AbstractVector{<:Int}, x::AbstractMatrix{<:Real}, op::MultiOrthoPoly)
    return evaluate(n, x, coeffs(op)...)
end
function evaluate(n::AbstractVector{<:Int}, x::AbstractVector{<:Real}, op::MultiOrthoPoly)
    return evaluate(n, reshape(x, 1, length(x)), op)
end

# using multi-index + multivariate
function evaluate(
        ind::AbstractMatrix{<:Int}, x::AbstractMatrix{<:Real},
        a::AbstractVector{<:AbstractVector{<:Real}},
        b::AbstractVector{<:AbstractVector{<:Real}}
    )
    vals = map(i -> evaluate(ind[i, :], x, a, b), Base.OneTo(size(ind, 1)))
    return hcat(vals...) |> transpose |> Matrix
end

function evaluate(ind::AbstractMatrix{<:Int}, x::AbstractMatrix{<:Real}, op::MultiOrthoPoly)
    return evaluate(ind, x, coeffs(op)...)
end
evaluate(x::AbstractMatrix{<:Real}, mop::MultiOrthoPoly) = evaluate(mop.ind, x, mop)

function evaluate(
        ind::AbstractMatrix{<:Int}, x::AbstractVector{<:Real},
        a::AbstractVector{<:AbstractVector{<:Real}},
        b::AbstractVector{<:AbstractVector{<:Real}}
    )
    return evaluate(ind, reshape(x, 1, length(x)), a, b)
end
function evaluate(ind::AbstractMatrix{<:Int}, x::AbstractVector{<:Real}, op::MultiOrthoPoly)
    return evaluate(ind, reshape(x, 1, length(x)), coeffs(op)...)
end
function evaluate(x::AbstractVector{<:Real}, mop::MultiOrthoPoly)
    return evaluate(mop.ind, reshape(x, 1, length(x)), mop)
end

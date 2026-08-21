export coeffs,
    nw,
    dim,
    deg,
    multi2uni,
    getentry,
    issymmetric,
    integrate

"""
    dim(basis)

Return the number of polynomials represented by an orthogonal basis.

# Arguments

- `basis`: a univariate [`AbstractOrthoPoly`](@ref) or a
  [`MultiOrthoPoly`](@ref).

For univariate bases this is `deg(basis) + 1`; for multivariate bases it is the
number of rows in the total-degree multi-index matrix.
"""
dim(op::AbstractOrthoPoly) = op.deg + 1
dim(mop::MultiOrthoPoly) = mop.dim

"""
    deg(basis)

Return the maximum represented polynomial degree.

# Arguments

- `basis`: a univariate [`AbstractOrthoPoly`](@ref) or a
  [`MultiOrthoPoly`](@ref).
"""
deg(op::AbstractOrthoPoly) = op.deg
deg(mop::MultiOrthoPoly) = mop.deg

"""
    nw(q::AbstractQuad)
    nw(op::AbstractOrthoPoly)
    nw(ops::AbstractVector)
    nw(mop::MultiOrthoPoly)

Return quadrature nodes and weights in matrix form. `nw(EmptyQuad())` returns
an empty `0 x 2` matrix.

# Arguments

- `q`, `op`, `ops`, `mop`: a quadrature rule, basis, collection of bases, or
  multivariate basis.

# Returns

A node/weight matrix for a univariate input, or one node and weight array per
component for a collection or multivariate basis.

# Examples

```jldoctest
julia> using PolyChaos

julia> size(nw(LegendreOrthoPoly(2)), 2)
2
```
"""
nw(::EmptyQuad) = Array{Float64}(undef, 0, 2)

function nw(quad::AbstractQuad)
    return [quad.nodes quad.weights]
end

nw(op::AbstractOrthoPoly) = nw(op.quad)

function nw(quads::Vector{<:AbstractQuad})
    nodes = [quad.nodes for quad in quads]
    weights = [quad.weights for quad in quads]
    return nodes, weights
end

function nw(ops::AbstractVector)
    quad = [op.quad for op in ops]
    return nw(quad)
end

nw(mop::MultiOrthoPoly) = nw(mop.uni)

"""
    coeffs(op::AbstractOrthoPoly)
    coeffs(ops::AbstractVector)
    coeffs(mop::MultiOrthoPoly)

Return the monic recurrence coefficients associated with an orthogonal basis.

# Arguments

- `op`, `ops`, `mop`: a basis, collection of univariate bases, or multivariate
  basis.

# Returns

A two-column `α`/`β` matrix for a univariate basis, or the corresponding
coefficient arrays for collections and multivariate bases.
"""
function coeffs(op::AbstractOrthoPoly)
    return [op.α op.β]
end

function coeffs(op::AbstractVector)
    a = [p.α for p in op]
    b = [p.β for p in op]
    return a, b
end

coeffs(mop::MultiOrthoPoly) = coeffs(mop.uni)

"""
    integrate(f::Function, nodes::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    integrate(f::Function, q::AbstractQuad)
    integrate(f::Function, op::AbstractOrthoPoly)

Integrate `f` using the supplied quadrature rule.

# Arguments

- `f`: scalar-valued integrand.
- `nodes`, `weights`: paired quadrature data.
- `q`, `op`: a quadrature rule or basis with an attached rule.

# Returns

The weighted sum of `f` evaluated at the quadrature nodes.

For example ``\\int_0^1 6x^5 = 1`` can be solved as follows:

```@repl
julia> opq = Uniform01OrthoPoly(3) # a quadrature rule is added by default

julia> integrate(x -> 6x^5, opq)
0.9999999999999993
```

!!! note


  - function ``f`` is assumed to return a scalar.
  - interval of integration is "hidden" in `nodes`.
"""
function integrate(
        f::Function, nodes::AbstractVector{<:Real},
        weights::AbstractVector{<:Real}
    )
    return dot(weights, f.(nodes))
end

function integrate(f::Function, quad::AbstractQuad)
    quad isa EmptyQuad &&
        throw(DomainError(quad, "supplied an empty quadrature"))
    return integrate(f, quad.nodes, quad.weights)
end

integrate(f::Function, op::AbstractOrthoPoly) = integrate(f, op.quad)

"""
```
integrate(f::Function, mop::MultiOrthoPoly)
```

Integrate a multivariate function `f` using tensor product quadrature from a
`MultiOrthoPoly`. The function `f` should accept the same number of arguments
as there are univariate orthogonal polynomials in `mop`.

For product measures, this computes the integral by evaluating `f` at all
combinations of quadrature nodes and weighting by the product of the
corresponding univariate weights.

# Arguments

- `f`: function accepting one scalar argument per coordinate.
- `mop`: multivariate basis with attached quadrature rules.

# Returns

The tensor-product quadrature estimate of the integral.

# Examples
```julia
op1 = GaussOrthoPoly(3)
op2 = Uniform01OrthoPoly(5)
mop = MultiOrthoPoly([op1, op2], 3)

# Integrate f(x,y) = x*y over the product measure
integrate((x, y) -> x * y, mop)
```
"""
function integrate(f::Function, mop::MultiOrthoPoly)
    nodes, weights = nw(mop)
    p = length(nodes)
    any(isempty, nodes) && throw(
        DomainError(
            mop,
            "one or more univariate orthogonal polynomials have empty quadrature; use addQuadrature=true"
        )
    )
    result = 0.0
    for idx in Iterators.product([eachindex(n) for n in nodes]...)
        node_vals = [nodes[d][idx[d]] for d in 1:p]
        weight_prod = prod(weights[d][idx[d]] for d in 1:p)
        result += weight_prod * f(node_vals...)
    end
    return result
end

"""
    issymmetric(m::AbstractMeasure)
    issymmetric(op::AbstractOrthoPoly)

Return whether the measure underlying `m` or `op` is symmetric about zero.

# Arguments

- `m`: measure implementing the [`AbstractMeasure`](@ref) contract.
- `op`: basis with an underlying measure.

# Returns

`true` when the measure is symmetric and `false` otherwise.
"""
issymmetric(m::AbstractMeasure) = m.symmetric
issymmetric(op::AbstractOrthoPoly) = issymmetric(op.measure)

"""
    multi2uni(a, ind)

Convert scalar-product basis indices to coordinate-wise univariate indices.

# Arguments

- `a`: nonnegative row indices into `ind`, using PolyChaos's zero-based basis
  convention.
- `ind`: total-degree multi-index matrix.

# Returns

A matrix whose column `j` is the univariate multi-index represented by `a[j]`.
"""
function multi2uni(a::AbstractVector{<:Int}, ind::AbstractMatrix{<:Int})
    minimum(a) < 0 && throw(DomainError(a, "no negative degrees allowed"))
    l, p = size(ind) # p-variate basis
    m = length(a) # dimension of scalar product
    l -= 1 # (l+1)-dimensional basis
    maximum(a) > l && throw(
        DomainError(
            a,
            "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"
        )
    )
    A = zeros(Int64, p, m)
    for (i, a_element) in enumerate(a)
        A[:, i] = ind[a_element + 1, :]
    end
    return A
end

"""
    getentry(a, T, ind, dim)

Return a sparse scalar-product tensor entry.

# Arguments

- `a`: zero-based polynomial indices; its length must equal `dim`.
- `T`: sparse tensor storage returned by [`computeTensorizedSP`](@ref).
- `ind`: total-degree multi-index matrix used to construct `T`.
- `dim`: scalar-product order.

`a` is sorted in place in descending order before lookup.
"""
function getentry(
        a::AbstractVector{<:Int}, T::SparseVector{<:Real, <:Int},
        ind::AbstractMatrix{<:Int}, dim::Int
    )
    m = length(a)
    l = size(ind, 1) - 1
    minimum(a) < 0 && throw(DomainError(a, "no negative degrees allowed"))
    maximum(a) > l && throw(
        DomainError(
            a,
            "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"
        )
    )
    m != dim && throw(
        DomainError(
            m,
            "length $m of provided index $a is inconsistent with dimension $(dim) of multi-index"
        )
    )
    # a .+= 1
    sort!(a; rev = true)

    sp_ind = 1 + evalpoly(l, a)
    return T[sp_ind]
end

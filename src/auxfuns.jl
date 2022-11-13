export coeffs,
       nw,
       dim,
       deg,
       multi2uni,
       getentry,
       issymmetric,
       integrate

dim(op::AbstractOrthoPoly) = op.deg + 1
dim(mop::MultiOrthoPoly) = mop.dim

deg(op::AbstractOrthoPoly) = op.deg
deg(mop::MultiOrthoPoly) = mop.deg

"""
```
nw(q::EmptyQuad)
nw(q::AbstractQuad)
nw(opq::AbstractOrthoPoly)
nw(opq::AbstractVector)
nw(mop::MultiOrthoPoly)
```
returns nodes and weights in matrix form
"""
nw(quad::typeof(EmptyQuad())) = Array{Float64}(undef, 0, 2)

function nw(quad::AbstractQuad)
    [quad.nodes quad.weights]
end

nw(op::AbstractOrthoPoly) = nw(op.quad)

function nw(quads::Vector{<:AbstractQuad})
    nodes = [quad.nodes for quad in quads]
    weights = [quad.weights for quad in quads]
    return nodes, weights
end

function nw(ops::AbstractVector)
    quad = [op.quad for op in ops]
    nw(quad)
end

nw(mop::MultiOrthoPoly) = nw(mop.uni)

"""
```
coeffs(op::AbstractOrthoPoly)
coeffs(op::AbstractVector)
coeffs(mop::MultiOrthoPoly)
```
returns recurrence coefficients of in matrix form
"""
function coeffs(op::AbstractOrthoPoly)
    [op.α op.β]
end

function coeffs(op::AbstractVector)
    a = [p.α for p in op]
    b = [p.β for p in op]
    return a, b
end

coeffs(mop::MultiOrthoPoly) = coeffs(mop.uni)

"""
```
integrate(f::Function,nodes::AbstractVector{<:Real},weights::AbstractVector{<:Real})
integrate(f::Function,q::AbstractQuad)
integrate(f::Function,opq::AbstractOrthoPoly)
```
integrate function `f` using quadrature rule specified via `nodes`, `weights`.
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
function integrate(f::Function, nodes::AbstractVector{<:Real},
                   weights::AbstractVector{<:Real})
    dot(weights, f.(nodes))
end

function integrate(f::Function, quad::AbstractQuad)
    typeof(quad) == typeof(EmptyQuad()) &&
        throw(DomainError(quad, "supplied an empty quadrature"))
    integrate(f, quad.nodes, quad.weights)
end

integrate(f::Function, op::AbstractOrthoPoly) = integrate(f, op.quad)

"""
```
issymmetric(m::AbstractMeasure)
issymmetric(op::AbstractOrthoPoly)
```
Is the measure symmetric (around any point in the domain)?
"""
issymmetric(m::AbstractMeasure) = m.symmetric
issymmetric(op::AbstractOrthoPoly) = issymmetric(op.measure)

function multi2uni(a::AbstractVector{<:Int}, ind::AbstractMatrix{<:Int})
    minimum(a) < 0 && throw(DomainError(a, "no negative degrees allowed"))
    l, p = size(ind) # p-variate basis
    m = length(a) # dimension of scalar product
    l -= 1 # (l+1)-dimensional basis
    maximum(a) > l && throw(DomainError(a,
                      "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"))
    A = zeros(Int64, p, m)
    for (i, a_element) in enumerate(a)
        A[:, i] = ind[a_element + 1, :]
    end
    return A
end

function getentry(a::AbstractVector{<:Int}, T::SparseVector{<:Real, <:Int},
                  ind::AbstractMatrix{<:Int}, dim::Int)
    m = length(a)
    l = size(ind, 1) - 1
    minimum(a) < 0 && throw(DomainError(a, "no negative degrees allowed"))
    maximum(a) > l && throw(DomainError(a,
                      "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"))
    m != dim && throw(DomainError(m,
                      "length $m of provided index $a is inconsistent with dimension $(dim) of multi-index"))
    # a .+= 1
    sort!(a; rev = true)

    sp_ind = 1 + evalpoly(l, a)
    return T[sp_ind]
end

export Tensor

"""
    Tensor(dim, basis)

Construct a sparse tensor of scalar products of order `dim` for an
orthogonal-polynomial basis.

# Arguments

- `dim`: positive scalar-product order.
- `basis`: a univariate [`AbstractOrthoPoly`](@ref) or multivariate
  [`MultiOrthoPoly`](@ref) with attached quadrature rules.

# Fields

- `dim`: scalar-product order.
- `T`: sparse vector of tensor entries.
- `get`: index-to-entry accessor.
- `op`: source orthogonal-polynomial basis.

# Examples

```jldoctest
julia> using PolyChaos

julia> Tensor(2, LegendreOrthoPoly(2; Nrec = 4)).dim
2
```
"""
struct Tensor{OP} <: AbstractTensor{OP}
    dim::Int          # "dimension"
    T::SparseVector{Float64, Int}
    get::Function
    op::OP
end
function Tensor(dim::Int, mop::MultiOrthoPoly)
    tensorEntries = computeTensorizedSP(dim, mop)
    getfun(ind) = getentry(ind, tensorEntries, mop.ind, dim)
    return Tensor(dim, tensorEntries, getfun, mop)
end
function Tensor(dim::Int, opq::AbstractOrthoPoly)
    tensorEntries = computeTensorizedSP(dim, opq)
    getfun(ind) = getentry(ind, tensorEntries, calculateMultiIndices(1, opq.deg), dim)
    return Tensor(dim, tensorEntries, getfun, opq)
end

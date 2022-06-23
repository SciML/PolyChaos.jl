export computeTensorizedSP

function computeTensorizedSP(m::Int,
                             α::AbstractVector{<:AbstractVector{<:Real}},
                             β::AbstractVector{<:AbstractVector{<:Real}},
                             nodes::AbstractVector{<:AbstractVector{<:Real}},
                             weights::AbstractVector{<:AbstractVector{<:Real}},
                             ind::AbstractMatrix{<:Int};
                             issymmetric::BitArray = falses(length(α)))
    m < 1 && throw(DomainError(m, "`dimension` has to be positive"))
    l, p = size(ind)
    l -= 1
    T = spzeros(reduce(+, [l^i for i in 0:m]))
    for tensor_ind in collect(with_replacement_combinations(0:l, m))
        index = 1 + reduce(+, [idx * l^(m - i) for (i, idx) in enumerate(tensor_ind)])
        T[index] = computeSP(tensor_ind, α, β, nodes, weights, ind; issymmetric = issymmetric)
    end
    return T
end
function computeTensorizedSP(m::Int, op::AbstractVector, ind::AbstractMatrix{<:Int})
    α, β = coeffs(op)
    nodes, weights = nw(op)
    return computeTensorizedSP(m, α, β, nodes, weights, ind; issymmetric = issymmetric.(op))
end

function computeTensorizedSP(m::Int, mop::MultiOrthoPoly)
    any([typeof(op.quad) == typeof(EmptyQuad()) for op in mop.uni]) &&
        throw(InconsistencyError("at least one quadrature rule missing"))
    return computeTensorizedSP(m, mop.uni, mop.ind)
end

function computeTensorizedSP(m::Int, op::AbstractOrthoPoly)
    typeof(op.quad) == typeof(EmptyQuad()) &&
        throw(InconsistencyError("no quadrature rule provided"))
    return computeTensorizedSP(m, [op], calculateMultiIndices(1, deg(op)))
end

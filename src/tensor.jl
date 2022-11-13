export computeTensorizedSP

function computeTensorizedSP(m::Integer,
                             α::AbstractVector{<:AbstractVector{<:Real}},
                             β::AbstractVector{<:AbstractVector{<:Real}},
                             nodes::AbstractVector{<:AbstractVector{<:Real}},
                             weights::AbstractVector{<:AbstractVector{<:Real}},
                             ind::AbstractMatrix{<:Integer};
                             issymmetric::BitArray = falses(length(α)))
    m < 1 && throw(DomainError(m, "`dimension` has to be positive"))
    l, p = size(ind)
    l -= 1
    s = 1
    for _ in 1:m
        s = muladd(s, l, 1)
    end
    T = spzeros(s)
    for tensor_ind in with_replacement_combinations(0:l, m)
        index = 1 + evalpoly(l, reverse!(tensor_ind))
        T[index] = computeSP(tensor_ind, α, β, nodes, weights, ind;
                             issymmetric = issymmetric)
    end
    return T
end
function computeTensorizedSP(m::Integer, op::AbstractVector, ind::AbstractMatrix{<:Integer})
    α, β = coeffs(op)
    nodes, weights = nw(op)
    computeTensorizedSP(m, α, β, nodes, weights, ind; issymmetric = issymmetric.(op))
end

function computeTensorizedSP(m::Integer, mop::MultiOrthoPoly)
    any([typeof(op.quad) == typeof(EmptyQuad()) for op in mop.uni]) &&
        throw(InconsistencyError("at least one quadrature rule missing"))
    computeTensorizedSP(m, mop.uni, mop.ind)
end

function computeTensorizedSP(m::Integer, op::AbstractOrthoPoly)
    typeof(op.quad) == typeof(EmptyQuad()) &&
        throw(InconsistencyError("no quadrature rule provided"))
    computeTensorizedSP(m, [op], calculateMultiIndices(1, deg(op)))
end

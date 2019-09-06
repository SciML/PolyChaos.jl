export computeTensorizedSP

function computeTensorizedSP(m::Int,
                             α::Vector{<:Vector{<:Real}},β::Vector{<:Vector{<:Real}},
                             nodes::Vector{<:Vector{<:Real}},weights::Vector{<:Vector{<:Real}},
                             ind::Matrix{<:Int};
                             issymmetric::BitArray=falses(length(α)))
    m < 1 && throw(DomainError(m, "`dimension` has to be positive"))
    l,p = size(ind)
    l -= 1
    T = spzeros(reduce(+,[l^i for i in 0:m]))
    for tensor_ind in collect(with_replacement_combinations(0:l, m))
        index = 1 + reduce(+,[idx*l^(m-i) for (i,idx) in enumerate(tensor_ind)])
        T[index] = computeSP(tensor_ind, α, β, nodes, weights, ind; issymmetric=issymmetric)
    end
    return T
end
function computeTensorizedSP(m::Int, op::Vector{<:AbstractOrthoPoly}, ind::Matrix{<:Int})
    α, β = coeffs(op)
    nodes, weights = nw(op)
    computeTensorizedSP(m, α, β, nodes, weights, ind; issymmetric=issymmetric.(op))
end

function computeTensorizedSP(m::Int, mop::MultiOrthoPoly)
    any([typeof(op.quad) == EmptyQuad for op in mop.uni]) && throw(InconsistencyError("at least one quadrature rule missing"))
    computeTensorizedSP(m, mop.uni, mop.ind)
end

function computeTensorizedSP(m::Int, op::AbstractOrthoPoly)
    typeof(op.quad) == EmptyQuad && throw(InconsistencyError("no quadrature rule provided"))
    computeTensorizedSP(m, [op], calculateMultiIndices(1, deg(op)))
end

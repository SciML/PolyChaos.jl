function computeTensorizedSP(m::Int64,
                             α::Vector{Vector{Float64}},β::Vector{Vector{Float64}},
                             nodes::Vector{Vector{Float64}},weights::Vector{Vector{Float64}},
                             ind::Matrix{Int64};
                             issymmetric::BitArray=falses(length(α)))
    @assert m>=1 "`dimension` has to be positive"
    l,p = size(ind)
    l-=1
    T = spzeros(reduce(+,[l^i for i in 0:m]))
    # @show length(collect(with_replacement_combinations(0:l,m)))
    for tensor_ind in collect(with_replacement_combinations(0:l,m))
        index = 1 + reduce(+,[idx*l^(m-i) for (i,idx) in enumerate(tensor_ind)])
        T[index] = computeSP(tensor_ind,α,β,nodes,weights,ind;issymmetric=issymmetric)
    end
    return T
end
function computeTensorizedSP(m::Int64,op::Vector{OrthoPoly},q::Vector{Quad},ind::Matrix{Int64})
    @assert length(op)==length(q) "Inconsistent lengh of polynomials and quadrature rules"
    @assert issymmetric.(op)==issymmetric.(q)
    α,β = coeffs(op)
    n,w = nw(q)
    computeTensorizedSP(m,α,β,n,w,ind;issymmetric=issymmetric.(op))
end
function computeTensorizedSP(m::Int64,mOP::MultiOrthoPoly)
    @assert typeof(mOP.uni)==Vector{OrthoPolyQ} "no quadrature rules provided"
    ops = [ mOP.uni[i].op for i=1:length(mOP.uni) ]
    quads = [ mOP.uni[i].quad for i=1:length(mOP.uni) ]
    computeTensorizedSP(m,ops,quads,mOP.ind)
end

function computeTensorizedSP(m::Int64,opq::OrthoPolyQ)
    computeTensorizedSP(m,[opq.op],[opq.quad],calculateMultiIndices(1,deg(opq)) )
end

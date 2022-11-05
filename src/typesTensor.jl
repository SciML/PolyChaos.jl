export Tensor

struct Tensor{OP} <: AbstractTensor{OP}
    dim::Int          # "dimension"
    T::SparseVector{Float64, Int}
    get::Function
    op::OP
end
function Tensor(dim::Int, mop::MultiOrthoPoly)
    tensorEntries = computeTensorizedSP(dim, mop)
    getfun(ind) = getentry(ind, tensorEntries, mop.ind, dim)
    Tensor(dim, tensorEntries, getfun, mop)
end
function Tensor(dim::Int, opq::AbstractOrthoPoly)
    tensorEntries = computeTensorizedSP(dim, opq)
    getfun(ind) = getentry(ind, tensorEntries, calculateMultiIndices(1, opq.deg), dim)
    Tensor(dim, tensorEntries, getfun, opq)
end

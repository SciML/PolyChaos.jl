export Tensor 

struct Tensor <: AbstractTensor
dim::Int          # "dimension"
T::SparseVector{Float64,Int}
get::Function
op::AbstractOrthoPoly

    # inner constructors
    function Tensor(dim::Int,mop::MultiOrthoPoly)
        tensorEntries = computeTensorizedSP(dim, mop)
        getfun(ind) = getentry(ind, tensorEntries, mop.ind, dim)
        new(dim, tensorEntries, getfun, mop)
    end
    function Tensor(dim::Int,opq::AbstractOrthoPoly)
        tensorEntries = computeTensorizedSP(dim, opq)
        getfun(ind) = getentry(ind, tensorEntries, calculateMultiIndices(1, opq.deg), dim)
        new(dim, tensorEntries, getfun, opq)
    end
end
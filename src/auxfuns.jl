export  coeffs,
        nw,
        dim,
        deg,
        multi2uni,
        getentry,
        issymmetric,
        integrate

dim(op::OrthoPoly)::Int64 = op.deg+1
dim(opq::OrthoPolyQ) = dim(opq.op)
dim(mop::MultiOrthoPoly) = mop.dim

deg(op::OrthoPoly)::Int64 = op.deg
deg(opq::OrthoPolyQ) = deg(opq.op)
deg(mop::MultiOrthoPoly) = mop.deg

"""
```
nw(q::Quad)
nw(opq::OrthoPolyQ)
nw(opq::Vector{OrthoPolyQ})
nw(mOP::MultiOrthoPoly)
```
returns nodes and weights in matrix form
"""
function nw(q::Quad)
    [q.nodes q.weights]
end
nw(opq::OrthoPolyQ) = nw(opq.quad)
function nw(qq::Vector{Quad})
    n = [ q.nodes for q in qq]
    w = [ q.weights for q in qq]
    return n,w
end
function nw(opq::Vector{OrthoPolyQ})
    q = [ p.quad for p in opq ]
    nw(q)
end

nw(mOP::MultiOrthoPoly) = nw(mOP.uni)

"""
```
coeffs(op::OrthoPoly)
coeffs(opq::OrthoPolyQ)
coeffs(op::Vector{OrthoPoly})
coeffs(opq::Vector{OrthoPolyQ})
coeffs(mop::MultiOrthoPoly)
```
returns recurrence coefficients of in matrix form
"""
function coeffs(op::OrthoPoly)
    [op.α op.β]
end
coeffs(opq::OrthoPolyQ) = coeffs(opq.op)

function coeffs(op::Vector{OrthoPoly})
    a = [ p.α for p in op]
    b = [ p.β for p in op]
    return a,b
end
function coeffs(opq::Vector{OrthoPolyQ})
    a = [ p.op.α for p in opq]
    b = [ p.op.β for p in opq]
    return a,b
end
coeffs(mop::MultiOrthoPoly) = coeffs(mop.uni)


"""
```
integrate(f::Function,nodes::Vector{Float64},weights::Vector{Float64})
integrate(f::Function,q::Quad)
integrate(f::Function,opq::OrthogonalPolyQ)
```
integrate function `f` using quadrature rule specified via `nodes`, `weights`

!!! note
- function ``f`` is assumed to return a scalar
- interval of integration is "hidden" in ``nodes``
"""
function integrate(f::Function,nodes::Vector{Float64},weights::Vector{Float64})::Float64
    dot(weights,f.(nodes))
end
integrate(f::Function,q::Quad)::Float64 = integrate(f,q.nodes,q.weights)
integrate(f::Function,opq::OrthoPolyQ)::Float64 = integrate(f,opq.quad)

"""
```
issymmetric(m::Measure)::Bool
issymmetric(op::OrthoPoly)::Bool
issymmetric(q::Quad)::Bool
issymmetric(opq::OrthoPolyQ)::Bool
```
is measure symmetric (around any point in the domain)?
"""
function issymmetric(m::Measure)::Bool
    m.symmetric
end
issymmetric(op::OrthoPoly)::Bool = issymmetric(op.meas)
issymmetric(q::Quad)::Bool = issymmetric(q.meas)
function issymmetric(opq::OrthoPolyQ)::Bool
    @assert issymmetric(opq.op)==issymmetric(opq.quad) "inconsistent symmetries"
    issymmetric(opq.op)
end

function multi2uni(a::Vector{Int64},ind::Matrix{Int64})
    @assert minimum(a)>=0 "no negative degrees allowed"
    l,p = size(ind) # p-variate basis
    m = length(a) # dimension of scalar product
    l=l-1 # (l+1)-dimensional basis
    @assert maximum(a)<=l "not enough elements in multi-index"
    A = zeros(Int64,p,m)
    for (i,a_) in enumerate(a)
        A[:,i] = ind[a_+1,:]
    end
    return A
end

function getentry(a::Vector{Int64},T::SparseVector{Float64,Int64},ind::Matrix{Int64},dim::Int64)
    m::Int64 = length(a)
    l::Int64 = size(ind,1)-1
    @assert minimum(a)>=0 "no negative degrees allowed"
    @assert maximum(a)<=l "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"
    @assert m==dim "length $m of provided index $a is inconsistent with dimension $(dim) of multi-index "
    # a .+= 1
    sort!(a)
    sp_ind::Int64 = 1 + reduce(+,[idx*l^(m-i) for (i,idx) in enumerate(a)])
    return T[sp_ind]::Float64
end

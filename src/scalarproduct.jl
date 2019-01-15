function computeSP(a::Vector{Int64},
                   α::Vector{Vector{Float64}},β::Vector{Vector{Float64}},
                   nodes::Vector{Vector{Float64}},weights::Vector{Vector{Float64}},
                   ind::Matrix{Int64};
                   issymmetric::BitArray=falses(length(α)))::Float64
    @assert minimum(a)>=0 "no negative degrees allowed"
    l,p = size(ind); # p-variate basis
    p==1 ? computeSP(a,α[1],β[1],nodes[1],weights[1];issymmetric=issymmetric[1]) : ()
    l-=1
    m = length(a) # m-dimensional scalar product
    @assert maximum(a)<=l "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"
    @assert length(α)==length(β)==length(nodes)==length(weights)==length(issymmetric)==p "inconsistent number of rec. coefficients and/or nodes/weights"
    # remove zero entries from index
    zero_inds = findall(x->x==0,a)
    a = a[setdiff(1:length(a),zero_inds)]
    if  length(a)==0
        return prod(β[i][1] for i=1:p)
    elseif length(a)==1
        return 0.
    else
        inds_uni = multi2uni(a,ind)
        # display(inds_uni)
        val = 1.
        for i=1:p
            # compute univariate scalar product
            v = computeSP(inds_uni[i,:],α[i],β[i],nodes[i],weights[i];issymmetric=issymmetric[i])
            v==0. ? (return 0.) : (val*=v)
        end
    end
    return val
end

function computeSP(a::Vector{Int64},op::Vector{OrthoPoly},q::Vector{Quad},ind::Matrix{Int64})
    @assert length(op)==length(q) "Inconsistent lengh of polynomials and quadrature rules"
    @assert issymmetric.(op)==issymmetric.(q)
    α,β = coeffs(op)
    n,w = nw(q)
    computeSP(a,α,β,n,w,ind;issymmetric=issymmetric.(op))
end

function computeSP(a::Vector{Int64},mOP::MultiOrthoPoly)
    @assert typeof(mOP.uni)==Vector{OrthoPolyQ} "no quadrature rules provided"
    ops = [ mOP.uni[i].op for i=1:length(mOP.uni) ]
    quads = [ mOP.uni[i].quad for i=1:length(mOP.uni) ]
    computeSP(a,ops,quads,mOP.ind)
end

"""
__Univariate__
```
computeSP(a::Vector{Int64},α::Vector{Float64},β::Vector{Float64},nodes::Vector{Float64},weights::Vector{Float64};issymmetric::Bool=false)
computeSP(a::Vector{Int64},op::OrthoPoly,q::Quad;issymmetric=issymmetric(op))
computeSP(a::Vector{Int64},opq::OrthoPolyQ)
```
__Multivariate__
```
computeSP( a::Vector{Int64},
           α::Vector{Vector{Float64}},β::Vector{Vector{Float64}},
           nodes::Vector{Vector{Float64}},weights::Vector{Vector{Float64}},
           ind::Matrix{Int64};
           issymmetric::BitArray=falses(length(α)))
computeSP(a::Vector{Int64},op::Vector{OrthoPoly},q::Vector{Quad},ind::Matrix{Int64})
computeSP(a::Vector{Int64},mOP::MultiOrthoPoly)
```

Computes the scalar product
```math
\\langle \\phi_{a_1},\\phi_{a_2},\\cdots,\\phi_{a_n} \\rangle,
```
where `n = length(a)`.
This requires to provide the recurrence coefficients `(α,β)` and the quadrature rule
`(nodes,weights)`, as well as the multi-index `ind`.
If provided via the keyword `issymmetric`, symmetry of the weight function is exploited.
All computations of the multivariate scalar products resort back to efficient computations
of the univariate scalar products.
Mathematically, this follows from Fubini's theorem.

The function is multiply dispatched to facilitate its use with `OrthoPolyQ` or a
suitable combination of `OrthoPoly` and its quadrature rule `Quad`.

!!! note
    - Zero entries of ``a`` are removed automatically to simplify computations, which follows from
    ```math
    \\langle \\phi_i, \\phi_j, \\phi_0,\\cdots,\\phi_0 \\rangle = \\langle \\phi_i, \\phi_j \\rangle,
    ```
    because `\\phi_0 = 1`.
    - It is checked whether enough quadrature points are supplied to solve the integral exactly.
"""
function computeSP(a_::Vector{Int64},α::Vector{Float64},β::Vector{Float64},nodes::Vector{Float64},weights::Vector{Float64};issymmetric::Bool=false)
    @assert minimum(a_)>=0 "no negative degrees allowed"
    # remove zero entries
    # this works because ``<Φ_i,Φ_j,Φ_0,...,Φ_0> = <Φ_i,Φ_j>``
    # hence, avoiding unnecessary operations that introduced numerical gibberish
    zero_inds = findall(x->x==0,a_)
    a = a_[setdiff(1:length(a_),zero_inds)]
    # simplification in case the measure is symmetric w.r.t t=0
    # exploit symmetry of the density, see Theorem 1.17 in W. Gautschi's book
    # and exploit the fact that
    issymmetric && isodd(sum(a)) ? (return 0.) : ()
    # Gauss quadrature rules have exactness 2N-1
    @assert Int(ceil(0.5*(sum(a)+1)))<=length(nodes) "not enough nodes to integrate exactly ($(length(nodes)) provided, where $(Int(ceil(0.5*(sum(a)+1)))) are needed)."
    if length(a)==0
        return β[1]
    elseif length(a)==1
        return 0.
    elseif length(a)==2
        a[1]==a[2] ? (return computeSP2(a[1],β)[end]) : return 0.
    else
        f = ones(Float64,length(nodes))
        for i=1:length(a)
            f = f.*evaluate(a[i],nodes,α,β)
        end
        return dot(weights,f)
    end
end
function computeSP(a::Vector{Int64},op::OrthoPoly,q::Quad)
    @assert issymmetric(op)==issymmetric(q) "inconsistent symmetries"
    computeSP(a,op.α,op.β,q.nodes,q.weights;issymmetric=issymmetric(op))
end
computeSP(a::Vector{Int64},opq::OrthoPolyQ) = computeSP(a,opq.op,opq.quad)

"""
    computeSP2(n::Int64,β::Vector{Float64})
    computeSP2(n::Int64,op::OrthoPoly) = computeSP2(n,op.β)
    computeSP2(op::OrthoPoly) = computeSP2(op.deg,op.β)
    computeSP2(opq::OrthoPolyQ) = computeSP2(opq.op)

Computes the `n` *regular* scalar products aka 2-norms of the orthogonal polynomials, namely
```math
\\|ϕ_i\\|^2 = \\langle \\phi_i,\\phi_i\\rangle \\quad \\forall i \\in \\{ 0,\\dots,n \\}.
```
Notice that only the values of `β` of the recurrence coefficients `(α,β)` are required.
The computation is based on equation (1.3.7) from Gautschi, W. "Orthogonal Polynomials: Computation and Approximation".
Whenever there exists an analytic expressions for `β`, this function should be used.

The function is multiply dispatched to facilitate its use with the composite types
`OrthoPoly` and `OrthoPolyQ`.
"""
function computeSP2(n::Int64,β::Vector{Float64})
    @assert n>=0 "can only compute scalar products for non-negative degrees"
    @assert length(β)>=n+1
    n==0 ? (return β[1]) : ()
    s = ones(Float64,n+1)
    s[1] = β[1]  # order 0computeSP2
    for i=2:n+1
        s[i] = s[i-1]*β[i]
    end
    return s
end
computeSP2(n::Int64,op::OrthoPoly) = computeSP2(n,op.β)
computeSP2(op::OrthoPoly) = computeSP2(op.deg,op.β)
computeSP2(opq::OrthoPolyQ) = computeSP2(opq.op)


# function computeSP2(n::Int64,α::Vector{Float64},β::Vector{Float64},nodes::Vector{Float64},weights::Vector{Float64})
#     @assert n>=0 "degree has to be non-negative"
#     @assert n<=length(α) "not enough recurrence coefficients"
#     Nquad, Nrec = length(nodes), length(α)
#     n==0 ? (return β[1]) : (dot(weights,evaluate(n,nodes,α,β).^2))
# end
# computeSP2(n::Int64,opq::OrthoPolyQ) = computeSP2(n,opq.op.α,opq.op.β,opq.quad.nodes,opq.quad.weights)

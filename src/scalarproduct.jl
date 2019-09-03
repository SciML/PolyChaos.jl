export  computeSP,
        computeSP2

function computeSP(a_::Vector{<:Int},
                   α::Vector{Vector{<:Real}},β::Vector{Vector{<:Real}},
                   nodes::Vector{Vector{<:Real}},weights::Vector{Vector{<:Real}},
                   ind::Matrix{<:Int};
                   issymmetric::BitArray=falses(length(α)),
                   zerotol::Float64=1e-10)
    minimum(a_) < 0 && throw(DomainError(minimum(a_),"no negative degrees allowed"))
    l, p = size(ind) # p-variate basis
    p == 1 && computeSP(a_,α[1],β[1],nodes[1],weights[1];issymmetric=issymmetric[1])
    l -= 1
    maximum(a_) > l && throw(DomainError(maximum(a_), "not enough elements in multi-index (requested: $(maximum(a)), max: $l)"))
    !(length(α)==length(β)==length(nodes)==length(weights)==length(issymmetric)==p) && throm(InconsistencyError("inconsistent number of rec. coefficients and/or nodes/weights"))

    a = Vector{Int64}()

    for aa in a_
        !iszero(aa) && push!(a,aa)
    end

    if iszero(length(a))
        return prod(β[i][1] for i in 1:p)
    elseif isone(length(a))
        return 0.
    else
        inds_uni = multi2uni(a,ind)
        val = 1.
        for i in 1:p
            # compute univariate scalar product
            @inbounds v = computeSP(inds_uni[i,:],α[i],β[i],nodes[i],weights[i];issymmetric=issymmetric[i])
            isapprox(v,0,atol=zerotol) ? (return 0.) : (val*=v)
        end
    end
    return val
end

function computeSP(a::Vector{<:Int}, op::Vector{<:AbstractOrthoPoly}, ind::Matrix{<:Int})
    α, β = coeffs(op)
    nodes, weights = nw(op)
    computeSP(a, α, β, nodes, weights, ind; issymmetric=issymmetric.(op))
end

computeSP(a::Vector{<:Int},mop::MultiOrthoPoly) = computeSP(a, mop.uni, [ op.quad for op in mop.uni ], mop.ind)

"""
__Univariate__
```
computeSP(a::Vector{<:Int},α::Vector{<:Real},β::Vector{<:Real},nodes::Vector{<:Real},weights::Vector{<:Real};issymmetric::Bool=false)
computeSP(a::Vector{<:Int},op::AbstractOrthoPoly;issymmetric=issymmetric(op))
```
__Multivariate__
```
computeSP( a::Vector{<:Int},
           α::Vector{Vector{<:Real}},β::Vector{Vector{<:Real}},
           nodes::Vector{Vector{<:Real}},weights::Vector{Vector{<:Real}},
           ind::Matrix{<:Int};
           issymmetric::BitArray=falses(length(α)))
computeSP(a::Vector{<:Int},op::Vector{<:AbstractOrthoPoly},ind::Matrix{<:Int})
computeSP(a::Vector{<:Int},mOP::MultiOrthoPoly)
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
function computeSP(a_::Vector{<:Int},α::Vector{<:Real},β::Vector{<:Real},nodes::Vector{<:Real},weights::Vector{<:Real};issymmetric::Bool=false)
    minimum(a_) < 0 && throw(DomainError(minimum(a_), "no negative degrees allowed"))
    # this works because ``<Φ_i,Φ_j,Φ_0,...,Φ_0> = <Φ_i,Φ_j>``
    # hence, avoiding unnecessary operations that introduce numerical gibberish
    # a = Vector{Int64}()
    # for aa in a_
    #     !iszero(aa) && push!(a,aa)
    # end
    a = filter(x -> x > 0, a_)
    # simplification in case the measure is symmetric w.r.t t=0
    # exploit symmetry of the density, see Theorem 1.17 in W. Gautschi's book
    # and exploit the fact that
    (issymmetric && isodd(sum(a))) && return 0.
    # Gauss quadrature rules have exactness 2N-1
    !(Int(ceil(0.5*(sum(a)+1)))<=length(nodes)) && throw(InconsistencyError("not enough nodes to integrate exactly ($(length(nodes)) provided, where $(Int(ceil(0.5*(sum(a)+1)))) are needed)")) 

    res = 
    if iszero(length(a))
        first(β)
    elseif isone(length(a))
        0.
    elseif length(a) == 2
        @inbounds a[1] == a[2] ? computeSP2(first(a),β)[end] : 0.
    else
        f = ones(Float64,length(nodes))
        @simd for i = 1:length(a)
            @inbounds f = f.*evaluate(a[i],nodes,α,β)
        end
        dot(weights,f)
    end
end

function computeSP(a::Vector{<:Int},op::AbstractOrthoPoly)
    computeSP(a,op.α,op.β,op.quad.nodes,op.quad.weights;issymmetric=issymmetric(op))
end

"""
    computeSP2(n::Int,β::Vector{<:Real})
    computeSP2(n::Int,op::AbstractOrthoPoly) = computeSP2(n,op.β)
    computeSP2(op::AbstractOrthoPoly) = computeSP2(op.deg,op.β)

Computes the `n` *regular* scalar products aka 2-norms of the orthogonal polynomials, namely
```math
\\|ϕ_i\\|^2 = \\langle \\phi_i,\\phi_i\\rangle \\quad \\forall i \\in \\{ 0,\\dots,n \\}.
```
Notice that only the values of `β` of the recurrence coefficients `(α,β)` are required.
The computation is based on equation (1.3.7) from Gautschi, W. "Orthogonal Polynomials: Computation and Approximation".
Whenever there exists an analytic expressions for `β`, this function should be used.

The function is multiply dispatched to facilitate its use with `AbstractOrthoPoly`.
"""
function computeSP2(n::Int,β::Vector{<:Real})
    n < 0 && throw(DomainError(n, "can only compute scalar products for non-negative degrees"))
    length(β) < n + 1 && throw(InconsistencyError("inconsistent length of β"))
    n == 0 && return β[1]
    s = ones(Float64,n+1)
    @inbounds s[1] = β[1]  # order 0computeSP2
    for i in 2:n+1
        @inbounds s[i] = s[i-1]*β[i]
    end
    return s
end
computeSP2(n::Int,op::AbstractOrthoPoly) = computeSP2(n,op.β)
computeSP2(op::AbstractOrthoPoly) = computeSP2(op.deg,op.β)


# function computeSP2(n::Int64,α::Vector{<:Real},β::Vector{<:Real},nodes::Vector{<:Real},weights::Vector{<:Real})
#     @assert n>=0 "degree has to be non-negative"
#     @assert n<=length(α) "not enough recurrence coefficients"
#     Nquad, Nrec = length(nodes), length(α)
#     n==0 ? (return β[1]) : (dot(weights,evaluate(n,nodes,α,β).^2))
# end
# computeSP2(n::Int64,opq::OrthoPolyQ) = computeSP2(n,opq.op.α,opq.op.β,opq.quad.nodes,opq.quad.weights)

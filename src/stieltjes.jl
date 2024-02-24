export stieltjes, lanczos, mcdiscretization

function removeZeroWeights(n::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    nw = [n w]
    inds = findall(w -> abs(w) >= eps(), nw[:, 2])
    # nw = sortrows(nw, by=x->x[1])
    nw = nw[inds, :]
    nw = sortslices(nw, dims = 1, by = x -> x[1])
    nw[:, 1], nw[:, 2]
end

"""
    stieltjes(N::Int,nodes_::AbstractVector{<:Real},weights_::AbstractVector{<:Real};removezeroweights::Bool=true)

Stieltjes procedure---Given the nodes and weights, the function
generates the first`N` recurrence coefficients of the
corresponding discrete orthogonal polynomials.

Set the Boolean `removezeroweights` to `true` if zero weights should be removed.
"""
function stieltjes(
        N::Int, nodes_::AbstractVector{<:Real}, weights_::AbstractVector{<:Real};
        removezeroweights::Bool = true)
    tiny = 10 * floatmin()
    huge = 0.1 * floatmax()
    α, β = zeros(Float64, N), zeros(Float64, N)
    nodes, weights = removezeroweights ? removeZeroWeights(nodes_, weights_) :
                     (nodes_, weights_)
    Ncap = length(nodes)
    @assert N > 0&&N <= Ncap "N is out of range."
    s0::Float64 = sum(weights)
    @inbounds α[1] = dot(nodes, weights) / s0
    @inbounds β[1] = s0
    N == 1 && return α, β
    p0, p1, p2 = zeros(Float64, Ncap), zeros(Float64, Ncap), ones(Float64, Ncap)
    for k in 1:(N - 1)
        p0 = p1
        p1 = p2
        @inbounds p2 = (nodes .- α[k]) .* p1 - β[k] * p0
        s1 = dot(weights, p2 .^ 2)
        s2 = dot(nodes, weights .* p2 .^ 2)
        # s1 = sum( weights[ν]*p2[ν]^2 for ν=1:Ncap )
        # s2 = sum( weights[ν]*nodes[ν]*p2[ν]^2 for ν=1:Ncap )
        abs(s1) < tiny && throw(DomainError(tiny,
            "Underflow in stieltjes() for k=$k; try using `removeZeroWeights`"))
        !(maximum(abs.(p2)) <= huge && abs(s2) <= huge) &&
            throw(DomainError(huge, "Overflow in stieltjes for k=$k"))
        @inbounds α[k + 1] = s2 / s1
        @inbounds β[k + 1] = s1 / s0
        s0 = s1
    end
    α, β
end

"""
    lanczos(N::Int,nodes::AbstractVector{<:Real},weights::AbstractVector{<:Real};removezeroweights::Bool=true)

Lanczos procedure---given the nodes and weights, the function
generates the first `N` recurrence coefficients of the
corresponding discrete orthogonal polynomials.

Set the Boolean `removezeroweights` to `true` if zero weights should be removed.

The script is adapted from the routine RKPW in
W.B. Gragg and W.J. Harrod, *The numerically stable
reconstruction of Jacobi matrices from spectral data*,
Numer. Math. 44 (1984), 317-335.
"""
function lanczos(N::Int, nodes::AbstractVector{<:Real}, weights::AbstractVector{<:Real};
        removezeroweights::Bool = true)
    !(length(nodes) == length(weights) > 0) &&
        throw(InconsistencyError("inconsistent number of nodes and weights"))
    nodes, weights = removezeroweights ? removeZeroWeights(nodes, weights) :
                     (nodes, weights)
    Ncap = length(nodes)
    (N <= 0 || N > Ncap) && throw(DomainError(N, "out of range"))
    p0 = copy(nodes)
    p1 = zeros(Float64, Ncap)
    @inbounds p1[1] = first(weights)
    for n in 1:(Ncap - 1)
        @inbounds pn = weights[n + 1]
        gam, sig, t = 1.0, 0.0, 0.0
        @inbounds xlam = nodes[n + 1]
        for k in 1:(n + 1)
            @inbounds ρ = p1[k] + pn
            tmp = gam * ρ
            tsig = sig
            gam, sig = ρ > 0 ? (p1[k] / ρ, pn / ρ) : (1.0, 0.0)
            @inbounds tk = sig * (p0[k] - xlam) - gam * t
            @inbounds p0[k] -= (tk - t)
            t = tk
            @inbounds pn = sig <= 0 ? tsig * p1[k] : (t^2) / sig
            tsig = sig
            @inbounds p1[k] = tmp
        end
    end
    @inbounds p0[Base.OneTo(N)], p1[Base.OneTo(N)]
end

"""
    mcdiscretization(N::Int,quads::Vector{},discretemeasure::AbstractMatrix{<:Real}=zeros(0,2);discretization::Function=stieltjes,Nmax::Integer=300,ε::Float64=1e-8,gaussquad::Bool=false)

This routine returns ``N`` recurrence coefficients of the polynomials that are
orthogonal relative to a weight function ``w`` that
is decomposed as a sum of ``m`` weights ``w_i`` with domains ``[a_i,b_i]`` for ``i=1,\\dots,m``,

```math
w(t) = \\sum_{i}^{m} w_i(t) \\quad \\text{with } \\operatorname{dom}(w_i) = [a_i, b_i].
```

For each weight ``w_i`` and its domain ``[a_i, b_i]`` the function `mcdiscretization()`
expects a quadrature rule of the form
nodes::AbstractVector{<:Real}, weights::AbstractVector{<:Real} = my_quad_i(N::Int)
all of which are stacked in the parameter `quad`
quad = [ my_quad_1, ..., my_quad_m ]
If the weight function has a discrete part (specified by `discretemeasure`)
it is added on to the discretized continuous weight function.

The function `mcdiscretization()` performs a sequence of discretizations of the
given weight ``w(t)``, each discretization being
followed by an application of the Stieltjes or Lanczos procedure
(keyword `discretization in [stieltjes, lanczos]`) to produce approximations to the desired recurrence
coefficients.
The function applies to each subinterval ``i`` an `N`-point quadrature rule (the ``i``th entry of `quad`)
to discretize the weight function ``w_i`` on that subinterval.
If the procedure converges to within a prescribed accuracy
`ε` before `N` reaches its maximum allowed value `Nmax`.
If the function does not converge, the function prompts an error message.

The keyword `gaussquad` should be set to `true` if Gauss quadrature rules
are available *for all* ``m`` weights ``w_i(t)`` with ``i = 1, \\dots, m``.

For further information, please see W. Gautschi "Orthogonal Polynomials: Approximation
and Computation", Section 2.2.4.
"""
function mcdiscretization(N::Int, quads::AbstractVector,
        discretemeasure::AbstractMatrix{<:Real} = zeros(0, 2);
        discretization::Function = stieltjes, Nmax::Integer = 300,
        ε::Float64 = 1e-8, gaussquad::Bool = false,
        removezeroweights::Bool = false)
    !(Nmax > 0 && Nmax > N) && throw(DomainError(Nmax, "invalid choice of Nmax=$Nmax"))
    ε <= 0 && throw(DomainError(ε, "invalid choice of ε"))
    discretization ∉ [stieltjes, lanczos] &&
        throw(DomainError(discretization, "unknown discretization"))
    length(quads) <= 0 && throw(InconsistencyError("no quadrature rule specified"))
    δ = !gaussquad ? 1 : 2
    # δ::Int=1
    # gaussquad ? δ=2 : ()
    mc, mp = size(quads, 1), size(discretemeasure, 1)
    Δ, kount = 1, -1
    α, β, b = zeros(Float64, N), zeros(Float64, N), ones(Float64, N)
    Mi = floor(Int64, (2 * N - 1) / δ)
    while any(abs.(β - b) .> ε * abs.(β))
        b = copy(β)
        kount += 1
        if kount > 1
            Δ = 2^(floor(Int64, kount / 5)) * N
        end
        # kount>1 ? Δ=2^(floor(Int64,kount/5))*N : ()
        Mi += Δ
        Mi > Nmax && error("Mi=$Mi exceeds Nmax=$Nmax")
        Ntot = mc * Mi
        xx, ww = zeros(Float64, Ntot), zeros(Float64, Ntot)
        for i in 1:mc
            nn = (i - 1) * Mi
            @inbounds x, w = quads[i](Mi)
            @inbounds xx[(nn + 1):(nn + Mi)], ww[(nn + 1):(nn + Mi)] = x, w
        end
        if mp > 0
            @inbounds xx[(Ntot + 1):(Ntot + mp)], ww[(Ntot + 1):(Ntot + mp)] = discretemeasure[
                :,
                1],
            discretemeasure[:,
                2]
        end
        α, β = discretization(N, xx, ww; removezeroweights = removezeroweights)
    end
    # printstyled("\nSuccess: ",color=:green)
    # print("converged after $kount iteration(s).\n\n")
    return α, β
end

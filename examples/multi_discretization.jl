using PolyChaos

function myquad(N::Int64, rec::Function, γ::Float64)
    α, β = rec(N)
    n, w = gauss(N, α, β)
    return n, γ * w
end

N = 10
γ = 0.5
quads = [n -> myquad(n, rm_chebyshev1, γ); n -> myquad(n, rm_legendre, 1 - γ)]

α, β = mcdiscretization(N, quads; Nmax = 300, ε = 1.0e-7, gaussquad = true)

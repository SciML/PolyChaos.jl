using PolyChaos, Test, LinearAlgebra
import QuadGK: quadgk

function integration_test(op::AbstractOrthoPoly, dim::Int, name::String; tol::Float64 = 1e-6,
                          dom::Tuple{<:Real,<:Real} = op.measure.dom)
    N, w, ind = op.deg, op.measure.w, zeros(Int64, dim)
    @testset "$name" begin
        for ind_ in Iterators.product([collect(1:N) for i in 1:dim]...)
            ind[:] .= ind_[:]
            s1 = computeSP(ind, op)
            f(t) = prod(evaluate(ind[i], t, op) for i in 1:dim) * w(t)
            s2 = quadgk(f, dom[1], dom[2])[1]
            if abs(s1) <= 1e-3
                @test isapprox(s1, s2; atol = tol)
            else
                @test isapprox(abs(s1 / s2), 1; atol = tol)
            end
        end
    end
end

N, Nrec, dim, tol = 3, 1000, 3, 1e-4

op = GaussOrthoPoly(N; Nrec = Nrec)
set_gaussian = integration_test(op, dim, "gaussian"; tol = tol, dom = (-10.0, 10.0))

α, β = 1.32, 3.23
op = Beta01OrthoPoly(N, α, β; Nrec = Nrec)
set_beta = integration_test(op, dim, "beta01"; tol = tol)

op = GammaOrthoPoly(N, 1, 1; Nrec = Nrec)
set_gamma = integration_test(op, dim, "gamma"; tol = tol)

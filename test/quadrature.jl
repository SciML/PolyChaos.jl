using PolyChaos, Test, LinearAlgebra, IterTools
import QuadGK: quadgk

function integration_test(opq::OrthoPolyQ,dim::Int64,name::String;tol::Float64=1e-6,
                          dom::Tuple{Float64,Float64}=opq.op.meas.dom)
    op::OrthoPoly = opq.op
    N = op.deg
    w::Function = op.meas.w
    ind = zeros(Int64,dim)
    @testset "$name" begin
        for ind_ in IterTools.product([collect(1:N) for i=1:dim]...)
            ind[:] .= ind_[:]
            s1 = computeSP(ind,opq)
            f(t) = prod(evaluate(ind[i],t,op) for i=1:dim)*w(t)
            # f(t) = evaluate(ind[1],t,op)*evaluate(ind[2],t,op)*evaluate(ind[3],t,op)*w(t)
            s2 = quadgk(f,dom[1],dom[2])[1]
            if abs(s1)<=1e-3
                @test isapprox(s1,s2;atol=tol)
            else
                @test isapprox(abs(s1/s2),1;atol=tol)
            end
        end
    end
end
##

N, Nrec, dim, tol = 4, 1000, 3, 1e-4

op = OrthoPoly("gaussian",N;Nrec=Nrec);
opq = OrthoPolyQ(op);
integration_test(opq,dim,"gaussian";tol=tol,dom=(-10.,10.))

α, β = 1.32, 3.23
op = OrthoPoly("beta01",N,Dict(:shape_a=>α,:shape_b=>β);Nrec=Nrec);
opq = OrthoPolyQ(op);
set_beta = integration_test(opq,dim,"beta01";tol=tol)

op = OrthoPoly("gamma",N,Dict(:shape=>1.0,:rate=>1.);Nrec=Nrec);
opq = OrthoPolyQ(op);
set_gamma = integration_test(opq,dim,"gamma";tol=tol,dom=(0.,Inf))

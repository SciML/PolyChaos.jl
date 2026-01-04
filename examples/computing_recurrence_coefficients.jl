# [ include("../src/"*s) for s in readdir("../src") ]
using PolyChaos
# using LinearAlgebra
# import FFTW
# import SpecialFunctions
# numerically compute recurrence coefficients for (almost) Gaussian density w(t)
N = 10
w(t) = exp(-t^2)
lb, ub = -Inf, Inf
@time α, β = rm_compute(w, lb, ub; Nquad = 200, Npoly = N)
# analytical solution
α_ana = zeros(N)
β_ana = [√π; [0.5 * k for k in 1:(N - 1)]]
# compare
display("Deviation for α: $(α - α_ana)")
display("Deviation for β: $(β - β_ana)")
## do the same for Chebyshev polynomials #4
v(t) = sqrt(1 - t) / sqrt(1 + t)
lb, ub = -1 + 1.0e-8, 1 - 1.0e-8
@time α, β = rm_compute(v, lb, ub; Nquad = 2000, Npoly = N)
# analytical solution
α_ana = [-0.5; zeros(N - 1)]
β_ana = [π; [0.25 for k in 1:(N - 1)]]
# compare
display("Deviation for α: $(α - α_ana)")
display("Deviation for β: $(β - β_ana)")

# Todos for PolyChaos.jl

## High priority
  - moments of order greater than 2 in `polynomial_chaos.jl` (needs `IterTools` for an elegant solution)
  - sampling (`polynomial_chaos.jl`)
    - rejection sampling
    - inverse CDF
  - plotting
  - gamma distribution needs to be double checked → construction of basis for feasible rates other than 1
  - code generation for random linear ODEs
  - ~~migrate to `Julia 1.1`~~
  - extend code tests
  - documentation
  - remove dependencies on `FastGaussQuadrature.jl`, for example: for few nodes, the golub welsch algorithm is used which we have implemented ourselves.

## Low priority
  - orthonormal polynomials
  - arbitrary polynomials

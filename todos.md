# Todos for PolyChaos.jl

## High priority
  - code:
    - multi-discretization from Gautschi ([see here](https://www.cs.purdue.edu/archives/2002/wxg/codes/mcdis.m))
    - moments of order greater than 2 in `polynomial_chaos.jl` (needs `IterTools` for an elegant solution)
    - sampling (`polynomial_chaos.jl`)
      - rejection sampling
      - inverse CDF
    - plotting by using recipes from `Plot.jl`
    - gamma distribution needs to be double checked → construction of basis for feasible rates other than 1
    - code generation for random linear ODEs
    - ~~migrate to `Julia 1.1`~~
    - ~~extend code tests~~
    - remove dependencies on `FastGaussQuadrature.jl`, for example: for few nodes, the golub welsch algorithm is used which we have implemented ourselves.
  - documentation:
      - have automated way to convert from notebook to `@repl` ⟶ see `docs/src/conversion.jl`
      - quickstart
      - examples:
        - quadrature rules (fejer, fejer2, clenshaw-curtis)
        - discretization procedures (stieltjes, lanczos), see `test/discretization.jl`
        - multi-discretization
        - `quadgp()`
        - optimizing Rosenbrock/quadratic function
        - DC-OPF with stochastic uncertainties


## Low priority
  - replace `PolyChaos` by `PolyChaos.jl`
  - orthonormal polynomials
  - arbitrary polynomials

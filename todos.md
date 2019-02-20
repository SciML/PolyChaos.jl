# Todos for PolyChaos.jl

## High priority
  - code:
    - moments of order greater than 2 in `polynomial_chaos.jl` (needs `IterTools` for an elegant solution)
    - sampling (`polynomial_chaos.jl`)
      - rejection sampling
      - inverse CDF --> transpile code from [here](https://github.com/dlfivefifty/InverseTransformSampling/blob/master/sample.m)
    - plotting by using recipes from `Plot.jl`
    - gamma distribution needs to be double checked → construction of basis for feasible rates other than 1
    - check code for efficiency following [guidelines](https://docs.julialang.org/en/v1/manual/performance-tips/): use `@inbounds` and `StaticArrays.jl`, see [here](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-annotations-1)
    - add [macros](https://docs.julialang.org/en/v1/manual/metaprogramming/#man-macros-1) to generate `OrthoPoly` and the rest
    - rewrite `evaluate.jl` without for loops, for example recursively or using `reduce`?
    - code generation for random linear ODEs
    - extend code testing
      - ade `codedev`
  - documentation:
    - make conversion to Latex work
    - add `quickstart.md`
    - examples:
      - quadrature rules (fejer, fejer2, clenshaw-curtis, gauss, radau, lobatto)
      - discretization procedures (stieltjes, lanczos), see `test/discretization.jl`
      - `quadgp()`
      - optimizing Rosenbrock/quadratic function
      - DC-OPF with stochastic uncertainties


## Low priority
  - orthonormal polynomials
  - arbitrary polynomials

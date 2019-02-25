```@setup mysetup
using PolyChaos
```
# Overview
PolyChaos is a collection of numerical routines for orthogonal polynomials written in the [Julia](https://julialang.org/) programming language.
Starting from some non-negative weight (aka an absolutely continuous nonnegative measure), PolyChaos allows
- to compute the coefficients for the monic three-term recurrence relation,
- to evaluate the orthogonal polynomials at arbitrary points,
- to compute the quadrature rule,
- to compute tensors of scalar products,
- to do all of the above in a multivariate setting (aka product measures).

If the weight function is a probability density function, PolyChaos further provides routines to compute [polynomial chaos expansions](https://en.wikipedia.org/wiki/Polynomial_chaos) (PCEs) of random variables with this very density function.
These routines allow
- to compute affine PCE coefficients for arbitrary densities,
- to compute moments,
- to compute the tensors of scalar products.

PolyChaos contains several *canonical* orthogonal polynomials such as Jacobi or Hermite polynomials.
For these, closed-form expressions and state-of-the art quadrature rules are used whenever possible.
However, a cornerstone principle of PolyChaos is to provide all the functionality for user-specific, arbitrary weights.

!!! note
    What PolyChaos is not (at least currently):
    - a self-contained introduction to orthogonal polynomials, quadrature rules and/or polynomial chaos expansions. We assume the user brings some experience to the table. However, over time we will focus on strengthening the tutorial charater of the package.
    - a symbolic toolbox
    - a replacement for [FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl)

## Installation
The package requires `Julia 1.0` or newer.
In `Julia` switch to the package manager
```julia
julia> ]
(v1.0) pkg> add PolyChaos
```
This will install PolyChaos and its dependencies.
Once that is done, load the package:
```julia
julia> using PolyChaos
```
That's it.

Let's take a look at a simple example.
We would like to solve the integral
```math
\int_0^1 6 x^5 \mathrm{d}x.
```
Exploiting the underlying uniform measure, the integration can be done exactly with a 3-point quadrature rule.
```@example mysetup
opq = OrthoPolyQ("uniform01",3)
integrate(x->6x^5,opq)
```

To get going with PolyChaos check out the tutorials such as the one on [numerical integration](@ref NumericalIntegration).
In case you are unfamiliar with orthogonal polynomials, perhaps [this background information](@ref MathematicalBackground) is of help.

## References
The code base of `PolyChaos` is partially based on Walter Gautschi's [Matlab suite of programs for generating orthogonal polynomials and related quadrature rules](https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html), with much of the theory presented in his book *Orthogonal Polynomials: Computation and Approximation* published in 2004 by the Oxford University Press.

For the theory of polynomial chaos expansion we mainly consulted T. J. Sullivan. *Introduction to Uncertainty Quantification*. Springer International Publishing Switzerland. 2015.

## Citing
Currently, there is no publication about `PolyChaos`.
Meanwhile, in case you find `PolyChaos` useful, feel free to get in touch, or simply participate in Github's gamification. ;)

## Collaboration
We are always looking for contributors.
If you are interested, just get in touch: tillmann [dot] muehlpfordt [at] kit [dot] edu.
Or just fork and/or star the repository.
Much appreciated.

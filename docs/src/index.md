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

A very specific application of orthogonal polynomials is the theory of [polynomial chaos expansions](https://en.wikipedia.org/wiki/Polynomial_chaos)---for which `PolyChaos` provides rich functionalities.
Loosely speaking, polynomial chaos is to random variables what Fourier series expansion is to periodic signals: a Hilbert space technique that allows to represent an infinite-dimensional mathematical object in terms of finitely many coefficients.

!!! note
    What PolyChaos is not (at least currently):
    - a self-contained introduction to orthogonal polynomials, quadrature rules and/or polynomial chaos expansions. We assume the user brings some experience to the table. However, over time we will focus on strengthening the tutorial charater of the package.
    - a symbolic toolbox
    - a replacement for [FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl)

## References
The code base of `PolyChaos` is partially based on Walter Gautschi's [Matlab suite of programs for generating orthogonal polynomials and related quadrature rules](https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html), with much of the theory presented in his book *Orthogonal Polynomials: Computation and Approximation* published in 2004 by the Oxford University Press.

For the theory of polynomial chaos expansion we mainly consulted T. J. Sullivan. *Introduction to Uncertainty Quantification*. Springer International Publishing Switzerland. 2015.

## Citing
Currently, there is no publication about `PolyChaos`.
Meanwhile, in case you find `PolyChaos` useful, feel free to get in touch, or simply participate in Github's gamification. ;)

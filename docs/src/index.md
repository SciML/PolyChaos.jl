# PolyChaos

PolyChaos is a collection of numerical routines for orthogonal polynomials written in the [Julia](https://julialang.org/) programming language.
Starting from some non-negative weight (aka an absolutely continuous nonnegative measure), PolyChaos allows
- to compute the coefficients for the monic three-term recurrence relation,
- to evaluate the orthogonal polynomials at arbitrary points,
- to compute the quadrature rule,
- to compute tensors of scalar products,
- to do all of the above in a multivariate setting (aka product measures).

If the weight function is a probability density function, PolyChaos further provides routines to compute polynomial chaos expansions (PCEs) of random variables with this very density function.
These routines allow
- to compute affine PCE coefficients for arbitrary densities,
- to compute moments of arbitrary degree,
- to compute the tensors of scalar products.

PolyChaos contains several *canonical* orthogonal polynomials such as Jacobi or Hermite polynomials.
For these, closed-form expressions and state-of-the art quadrature rules are used whenever possible.
However, a cornerstone principle of PolyChaos is to provide all the functionality for user-specific, arbitrary weights.

PolyChaos comes with several *canonical* orthogonal polynomials, for example Jacobi polynomials,


What PolyChaos is not (at least currently):
- symbolic toolbox
- replacement for [FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl)


A very specific application of orthogonal polynomials is the theory of polynomial chaos expansion---for which PolyChaos provides rich functionalities.
Loosely speaking, polynomial chaos is to random variables what Fourier series expansion is to periodic signals: a Hilbert space technique that allows to represent an infinite-dimensional mathematical object in terms of finitely many coefficients.

that combines orthogonal polynomials, quadrature rules and polynomial chaos expansions.

```@contents
Depth = 3
```
## Types


## Functions

### Recurrence Coefficients for Monic Orthogonal Polynomials
The functions below provide analytic expressions for the recurrence coefficients of common orthogonal polynomials.
All of these provide *monic orthogonal polynomials* relative to the weights.

!!! note
    The number `N` of recurrence coefficients has to be positive for all functions below.
```@docs
r_scale(c::Float64,a::Vector{Float64},b::Vector{Float64})
rm_compute(weight::Function,lb::Float64,ub::Float64;Npoly::Int64=4,Nquad::Int64=10,quadrature::Function=clenshaw_curtis)
rm_logistic(N::Int)
rm_hermite(N::Int,mu::Float64)
rm_hermite_prob(N::Int)
rm_laguerre(N::Int,a::Float64)
rm_legendre(N::Int)
rm_legendre01(N::Int)
rm_jacobi(N::Int,a::Float64,b::Float64)
rm_jacobi01(N::Int,a::Float64,b::Float64)
rm_meixner_pollaczek(N::Int,lambda::Float64,phi::Float64)
```

## Scalar Products
```@docs
computeSP2
computeSP
```

## Evaluating Orthogonal Polynomials
```@docs
evaluate
```

## Quadrature Rules
```@docs
fejer
fejer2
clenshaw_curtis
quadpts_beta01
quadpts_gamma
quadpts_gaussian
quadpts_logistic
quadpts_uniform01
quadgp
```


## Index

```@index
```

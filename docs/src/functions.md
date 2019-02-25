# Functions

!!! note
    The core interface of (so we hope) all essential functions are *not* dependent on specialized types such as `OrthoPoly`/`OrthoPolyQ`.
    Having said that, for exactly those essential functions there exist overloaded functions that accept specialized types such as `OrthoPoly`/`OrthoPolyQ` as arguments.

    Too abstract?
    For example, the function `evaluate()` that evaluates a polynomial of degree `n` at points `x` has the core interface
    ```
        evaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})
    ```
    where `a` and `b` are the vectors of recurrence coefficients.
    For simplicity, there also exists the interface
    ```
        evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly)
    ```
    which is defined as
    ```
        evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly) = evaluate(n,x,op.α,op.β)
    ```
    So fret not upon the encounter of multiply-dispatched versions of the same thing. It's there to simplify your life (so we hope).

    The idea of this approach is to make it simpler for others to copy and paste code snippets and use them in their own work.

List of all functions in `PolyChaos`.

```@index
```

## Recurrence Coefficients for Monic Orthogonal Polynomials
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
stieltjes
lanczos
mcdiscretization
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
quadgp
gauss
radau
radau_jacobi
radau_laguerre
lobatto
lobatto_jacobi
```

## Polynomial Chaos
```@docs
mean
var
std
sampleMeasure
evaluatePCE
samplePCE
calculateAffinePCE
convert2affinePCE
```

## Auxiliary Functions
```@docs
nw
coeffs
integrate
issymmetric
```

# Functions

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

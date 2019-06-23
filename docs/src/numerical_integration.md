# [Numerical Integration](@id NumericalIntegration)
The goal of this tutorial is to solve an integral using Gauss quadrature,
```math
I := \int_{0}^{1} f(t) \mathrm{d} t \approx \sum_{k=1}^n w_k f(t_k),
```
where we choose $f(t) = \sin(t)$, and $n = 5$.

Make sure to check out [this tutorial](@ref QuadratureRules) too.

### Variant 0
```jldoctest mylabel
julia> using PolyChaos

julia> n = 5;

julia> f(t) = sin(t);

julia> opq = OrthoPolyQ("uniform01",n);

julia> I0 = integrate(f,opq)
0.4596976941320484

julia> print("Numerical error: $(abs(1-cos(1)-I0))")
Numerical error: 1.8818280267396403e-13
```
with negligible numerical errors.

### Variant 1
Let us  now solve the same problem, while elaborating what is going on under the hood.
At first, we load the package by calling
```@repl
using PolyChaos
```
Now we define a measure, specifically the uniform measure $\mathrm{d}\lambda(t) = w(t) \mathrm{d} t$ with the weight $w$ defined as
```math
  w: \mathcal{W} = [0,1] \rightarrow \mathbb{R}, \quad w(t) = 1.
```
This measure can be defined using the composite type `Measure`:
```jldoctest mylabel
julia> m = Measure("uniform01")
Measure("uniform01", PolyChaos.w_uniform01, (0.0, 1.0), true, Dict{Any,Any}())
```
Next, we need to compute the quadrature rule relative to the uniform measure.
To do this we use the composite type `Quad`.

```jldoctest mylabel
julia> q1 = Quad(n-1,m)
Quad("quadgp", 4, [1.0, 0.853553, 0.5, 0.146447, 0.0], [0.0333333, 0.266667, 0.4, 0.266667, 0.0333333], Measure("uniform01", PolyChaos.w_uniform01, (0.0, 1.0), true, Dict{Any,Any}()))

julia> nw(q1)
5×2 Array{Float64,2}:
 1.0       0.0333333
 0.853553  0.266667
 0.5       0.4
 0.146447  0.266667
 0.0       0.0333333
```
This creates a quadrature rule `q` named `"myq"` with `n-1` nodes and weights relative to the measure `m`.
The function `nw()` prints the nodes and weights.
To solve the integral we call `integrate()`
```jldoctest mylabel
julia> I1 = integrate(f,q1)
0.4596977927043755

julia> print("Numerical error: $(abs(1-cos(1)-I1))")
Numerical error: 9.857251526135258e-8
```

### Variant 2
There is another variant to solve the integral, which computes the quadrature rule based on the recurrence coefficients of the polynomials that are orthogonal relative to the measure `m`.
First, we compute the orthogonal polynomials using the composite type `OrthoPoly`.
```jldoctest mylabel
julia> op = OrthoPoly("uniform01",n)
OrthoPoly("uniform01", 5, [0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [1.0, 0.0833333, 0.0666667, 0.0642857, 0.0634921, 0.0631313], Measure("uniform01", PolyChaos.w_uniform01, (0.0, 1.0), true, Dict{Any,Any}()))

julia> coeffs(op)
6×2 Array{Float64,2}:
 0.5  1.0
 0.5  0.0833333
 0.5  0.0666667
 0.5  0.0642857
 0.5  0.0634921
 0.5  0.0631313
```
The resulting system of orthogonal polynomials is characterized by its recursion coefficients $(\alpha, \beta)$, which can be extracted with the function `coeffs()`.

Now, the quadrature rule can be constructed based on `op`, and the integral be solved.
```jldoctest mylabel
julia> q2 = Quad(n,op)
Quad("golubwelsch", 5, [0.0469101, 0.230765, 0.5, 0.769235, 0.95309], [0.118463, 0.239314, 0.284444, 0.239314, 0.118463], Measure("uniform01", PolyChaos.w_uniform01, (0.0, 1.0), true, Dict{Any,Any}()))

julia> nw(q2)
5×2 Array{Float64,2}:
 0.0469101  0.118463
 0.230765   0.239314
 0.5        0.284444
 0.769235   0.239314
 0.95309    0.118463

julia> I2 = integrate(f,q2)
0.4596976941320484

julia> print("Numerical error: $(abs(1-cos(1)-I2))")
Numerical error: 1.8818280267396403e-13
```

### Comparison
We see that the different variants provide slightly different results:
```jldoctest mylabel
julia> 1-cos(1) .- [I0 I1 I2]
1×3 Array{Float64,2}:
 -1.88183e-13  -9.85725e-8  -1.88183e-13
```
with `I0` and `I2` being the same and more accurate than `I1`.
The increased accuracy is based on the fact that for `I0` and `I2` the quadrature rules are based on the recursion coefficients of the underlying orthogonal polynomials.
The quadrature for `I1` is based on an general-purpose method that can be significantly less accurate, see also [the next tutorial](@ref QuadratureRules).

# [Numerical Integration](@id NumericalIntegration)
```@setup usePolyChaos
using PolyChaos
n = 5; f(t) = sin(t)
opq = OrthoPolyQ("uniform01",n-1);
I0 = integrate(f,opq);
m = Measure("uniform01");
q = Quad(n-1,m);
I1 = integrate(f,q)
op = OrthoPoly("uniform01",n-1)
q = Quad(n,op)
I2 = integrate(f,q)
```
The goal of this tutorial is to solve an integral using Gauss quadrature,
```math
I := \int_{0}^{1} f(t) \mathrm{d} t \approx \sum_{k=1}^n w_k f(t_k),
```
where we choose $f(t) = \sin(t)$, and $n = 5$.

### Variant 0
```@repl
using PolyChaos
n = 5;
f(t) = sin(t);
opq = OrthoPolyQ("uniform01",n-1);
I0 = integrate(f,opq)
print("Numerical error: $(abs(1-cos(1)-I0))")
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
```@repl usePolyChaos
m = Measure("uniform01");
```
Next, we need to compute the quadrature rule relative to the uniform measure.
To do this we use the composite type `Quad`.
```@repl usePolyChaos
q1 = Quad(n-1,m);
nw(q)
```
This creates a quadrature rule `q` named `"myq"` with `n-1` nodes and weights relative to the measure `m`.
The function `nw()` prints the nodes and weights.
To solve the integral we call `integrate()`
```@repl usePolyChaos
I1 = integrate(f,q1)
print("Numerical error: $(abs(1-cos(1)-I1))")
```

### Variant 2
There is another variant to solve the integral, which computes the quadrature rule based on the recurrence coefficients of the polynomials that are orthogonal relative to the measure `m`.
First, we compute the orthogonal polynomials using the composite type `OrthoPoly`.
```@repl usePolyChaos
op = OrthoPoly("uniform01",n-1);
coeffs(op)
```
The resulting system of orthogonal polynomials is characterized by its recursion coefficients $(\alpha, \beta)$, which can be extracted with the function `coeffs()`.

Now, the quadrature rule can be constructed based on `op`, and the integral be solved.
```@repl usePolyChaos
q2 = Quad(n,op)
nw(q)
I2 = integrate(f,q2)
print("Numerical error: $(abs(1-cos(1)-I2))")
```

### Comparison
We see that the different variants provide slightly different results:
```@repl usePolyChaos
1-cos(1) .- [I0 I1 I2]
```
with `I0` and `I2` being the same and more accurate than `I1`.
The increased accuracy is based on the fact that for `I0` and `I2` the quadrature rules are based on the recursion coefficients of the underlying orthogonal polynomials.
The quadrature for `I1` is based on an general-purpose method that can be significantly less accurate.

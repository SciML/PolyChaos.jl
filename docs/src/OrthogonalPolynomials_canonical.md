# [Univariate Monic Orthogonal Polynomials](@id UnivariateMonicOrthogonalPolynomials)
Univariate monic orthogonal polynomials make up the core building block of the package.
These are real polynomials $\{ \pi_k \}_{k \geq 0}$, which are univariate $\pi_k: \mathbb{R} \rightarrow \mathbb{R}$ and orthogonal relative to a nonnegative weight function $w: \mathbb{R} \rightarrow \mathbb{R}_{\geq 0}$, and which have a leading coefficient equal to one:
```math
\begin{aligned}
\pi_k(t) &= t^k + a_{k-1} t^{k-1} + \dots + a_1 t + a_0 \quad \forall k = 0, 1, \dots \\
\langle \pi_k, \pi_l \rangle &= \int_{\mathbb{R}} \pi_k(t) \pi_l(t) w(t) \mathrm{d}t =
\begin{cases}
0 & k \neq l, \text{ and }k,l \geq 0 \\
\| \pi_k \|^2 > 0 & k = l \geq 0
\end{cases}
\end{aligned}
```

These univariate monic orthogonal polynomials satisfy the paramount three-term recurrence relation
```math
\begin{aligned}
\pi_{k+1}(t) &= (t - \alpha_k) \pi_k(t) - \beta_k \pi_{k-1}(t), \quad k= 0, 1, \dots, \\
\pi_o(t) &= 1, \\
\pi_{-1}(t) &= 0.
\end{aligned}
```

Hence, every system of $n$ univariate monic orthogonal polynomials $\{ \pi_k \}_{k=0}^n$ is isomorphic to its recurrence coefficients $\{ \alpha_k, \beta_k \}_{k=0}^n$.


## Classical Orthogonal Polynomials

The so-called *classical* orthogonal polynomials are polynomials named after famous mathematicians who each discovered a special family of orthogonal polynomials, for example [Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials) or [Jacobi polynomials](https://en.wikipedia.org/wiki/Jacobi_polynomials).
For *classical* orthogonal polynomials there exist closed-form expressions of---among others---the recurrence coefficients.
Also quadrature rules for *classical* orthogonal polynomials are well-studied (with dedicated packages such as [FastGaussQuadrature.jl](https://github.com/ajt60gaibb/FastGaussQuadrature.jl).
However, more often than not these *classical* orthogonal polynomials are neither monic nor orthogonal, hence not normalized in any sense.
For example, there is a distinction between the [*probabilists'* Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition) and the [*physicists'* Hermite polynomials](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition).
The difference is in the weight function $w(t)$ relative to which the polynomials are orthogonal:
```math
\begin{aligned}
&\text{Probabilists':} &&&w(t) = \frac{1}{\sqrt{2 \pi}} \, \exp \left( - \frac{t^2}{2} \right) \\
&\text{Physicists':} &&&w(t) =  \exp \left( - t^2 \right).
\end{aligned}
```

To streamline the computations, all *classical* orthogonal polynomials are converted to __monic__ orthogonal polynomials (for which, of course, the closed-form expressions persist).
Currently, the following weight functions (hence *classical* orthogonal polynomials) are supported:

| Name | Weight $w(t)$| Parameters | Support| *Classical* polynomial |
| --- | --- | --- | --- | --- |
| `hermite` | $ \exp \left( - t^2 \right)$ | - | $(-\infty, \infty)$ | Hermite |
| `genhermite` | $ \lvert t \rvert^{2 \mu}\exp \left( - t^2 \right)$ | $\mu > -\frac{1}{2}$ | $(-\infty, \infty)$ | Generalized Hermite |
| `legendre` | $1$ | - | $(-1,1)$ | Legendre
| `jacobi` | $(1-t)^{\alpha} (1+t)^{\beta}$ | $\alpha, \beta > -1$ | $(-1,1)$ | Jacobi |
| `laguerre` | $\exp(-t)$ | - | $(0,\infty)$ | Laguerre |
| `genlaguerre` | $t^{\alpha}\exp(-t)$ | $\alpha>-1$ | $(0,\infty)$ | Generalized Laguerre |
| `meixnerpollaczek` | $\frac{1}{2 \pi} \exp((2\phi-\pi)t) \lvert\Gamma(\lambda + \mathrm{i}t)\rvert^2$ |$\lambda > 0, 0<\phi<\pi$ | $(-\infty,\infty)$ | Meixner-Pollaczek


Additionally, the following weight functions that are equivalent to probability density functions are supported:

| Name | Weight $w(t)$| Parameters | Support| *Classical* polynomial |
| --- | --- | --- | --- | --- |
| `gaussian` | $\frac{1}{\sqrt{2 \pi}} \, \exp \left( - \frac{t^2}{2} \right)$ | - | $(-\infty, \infty)$ | Probabilists' Hermite |
| `uniform01` | $1$ | - | $(0,1)$ |  Legendre
| `beta01` | $\frac{1}{B(\alpha,\beta)} \, t^{\alpha-1} (1-t)^{\beta-1}$ |$\alpha, \beta > 0$ | $(0,1)$ | Jacobi |
| `gamma` | $\frac{\beta^\alpha}{\Gamma(\alpha)} t^{\alpha-1} \exp(-\beta t)$ | $\alpha, \beta > 0$ | $(0,\infty)$ | Laguerre |
| `logistic` | $\frac{\exp(-t)}{(1+\exp(-t))^2}$ | - | $(-\infty,\infty)$ | -

To generate the orthogonal polynomials up to maximum degree `deg`, simply call


```julia
using PolyChaos
deg = 4
op = OrthoPoly("gaussian",deg)
```

This generates `op`as an `OrthoPoly` type with the underlying Gaussian measure `op.meas`.
The recurrence coefficients are accessible via `coeffs()`.


```julia
coeffs(op)
```

By default, the constructor for `OrthoPoly` generates `deg+1` recurrence coefficients.
Sometimes, some other number `Nrec` may be required.
This is why `Nrec` is a keyword for the constructor `OrthoPoly`.


```julia
N = 100
op_ = OrthoPoly("logistic",deg;Nrec=N)
```

Let's check whether we truly have more coefficients:


```julia
size(coeffs(op_),1)==N
```

## Arbitrary Weights
If you are given a weight function $w$ that does not belong to the Table above, it is still possible to generate the respective univariate monic orthogonal polynomials.
First, we define the measure by specifying a name, the weight, the support, symmetry, and parameters


```julia
supp = (-1,1)
function w(t)
    supp[1]<=t<=supp[2] ? (1. + t) : error("$t not in support")
end
my_meas = Measure("my_meas",w,supp,false,Dict())
```

Notice: it is advisable to define the weight such that an error is thrown for arguments outside of the support.

Now, we want to construct the univariate monic orthogonal polynomials up to degree `deg` relative to `m`.
The constructor is


```julia
my_op = OrthoPoly("my_op",deg,my_meas;Nquad=200)
```

By default, the recurrence coefficients are computed using the [Stieltjes procuedure](https://warwick.ac.uk/fac/sci/maths/research/grants/equip/grouplunch/1985Gautschi.pdf) with [Clenshaw-Curtis](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature) quadrature (with `Nquad` nodes and weights).
Hence, the choice of `Nquad` influences accuracy.


## [Multivariate Monic Orthogonal Polynomials](@id MultivariateMonicOrthogonalPolynomials)
Suppose we have $p$ systems of univariate monic orthogonal polynomials,
```math
\{ \pi_k^{(1)} \}_{k\geq 0}, \: \{ \pi_k^{(2)} \}_{k\geq 0}, \dots, \{ \pi_k^{(p)} \}_{k\geq 0},
```
each system being orthogonal relative to the weights $w^{(1)}, w^{(2)}, \dots, w^{(p)}$ with supports $\mathcal{W}^{(1)}, \mathcal{W}^{(2)}, \dots, \mathcal{W}^{(p)}$.
Also, let $d^{(i)}$ be the maximum degree of the $i$-th system of univariate orthogonal polynomials.
We would like to construct a $p$-variate monic basis $\{ \psi_k \}_{k \geq 0}$ with $\psi: \mathbb{R}^p \rightarrow \mathbb{R}$ of degree at most $0 \leq d \leq \min_{i=1,\dots,k}\{ d^{(i)}\}$.
Further, this basis shall be orthogonal relative to the product measure $w: \mathcal{W} = \mathcal{W}^{(1)} \otimes \mathcal{W}^{(2)} \mathcal{W}^{(1)} \cdots \otimes \mathcal{W}^{(p)} \rightarrow \mathbb{R}_{\geq 0}$ given by
```math
w(t) = \prod_{i=1}^{p} w^{(i)}(t_i),
```
hence satisfies
```math
\langle \psi_k, \psi_l \rangle = \int_{\mathcal{W}} \psi_k(t) \psi_l(t) w(t) \mathrm{d} t =
\begin{cases}
0 & k \neq l, \text{ and }k,l \geq 0 \\
\| \psi_k \|^2 > 0 & k = l \geq 0
\end{cases}
```

For this, there exists the composite struct `MultiOrthoPoly`.
Let's consider an example where we mix *classical* orthogonal polynomials with an arbitrary weight.


```julia
deg = [3,5,6,4]
d = minimum(deg)

op1 = OrthoPoly("gaussian",deg[1])
op2 = OrthoPoly("uniform01",deg[2])
op3 = OrthoPoly("beta01",deg[3],Dict(:shape_a=>2,:shape_b=>1.2))
ops = [op1,op2,op3,my_op]
mop = MultiOrthoPoly(ops,d)
```

The total number of  basis polynomials is stored in the field `dim`.
The univariate basis polynomials making up the multivariate basis are stored in the field `uni`.



```julia
mop.uni
```

The field `ind` contains the multi-index, i.e. row $i$ stores what combination of univariate polynomials makes up the $i$-th multivariate polynomial.
For example,


```julia
i = 11
mop.ind[i+1,:]
```

translates mathematically to
```math
\psi_{11}(t) = \pi_0^{(1)}(t_1) \pi_1^{(2)}(t_2) \pi_0^{(3)}(t_3) \pi_1^{(4)}(t_4).
```

Notice that there is an offset by one, because the basis counting starts at 0, but Julia is 1-indexed.
The underlying measure of `mop` is now of type `MultiMeasure`, and stored in the field `meas`


```julia
mop.meas
```

The weight $w$ can be evaluated as expected


```julia
mop.meas.w(0.5*ones(length(ops)))
```

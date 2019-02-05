# Type Hierarchy

Let's look at the types `PolyChaos` provides.
The high-level perspective looks as such:
![type hierarchy](assets/TypeHierarchy_transparent.png)

The left hand side covers types related to univariate measures; the right hand side covers multivariate measures.
An arrow beginning at one type and ending at another type means that the beginning type is a field of the ending type.
For example, the type `OrthoPoly` has a field of type `Measure`; the type `OrthoPolyQ` has a field of type `OrthoPoly` and a field of type `Quad`, and so on.
Let's begin with the univariate case.

!!! note
    If you are unfamiliar with the mathematical background of orthogonal polynomials, please consult [this tutorial](@ref MathematicalBackground).


## Measure
It all begins with a measure, more specifically absolutely continuous measures.
What are the fields of such a type `measure`?


| Field | Meaning |
| --- | --- |
| `name::String` | Name of measure|
| `w::Function` | Weight function $w: \Omega \rightarrow \mathbb{R}$ |
| `dom::Tuple{Float64,Float64}` | Domain $ \Omega$|
| `symmetric::Bool` | Is $w$ symmetric relative to some $m \in \Omega$, hence $w(m-x) = w(m+x)$ for all $x \in \Omega$? |
| `pars::Dict` | Additional parameters (e.g. shape parameters for Beta distribution |


They are a `name`, a weight function $w: \Omega \rightarrow \mathbb{R}$ with domain $\Omega$ (`dom`).
If the weight function is symmetric relative to some $m \in \Omega$, the field `symmetric` should be set to `true`.
Symmetry relative to $m$ means that
```math
\forall x \in \Omega: \quad w(m-x) = w(m+x).
```
For example, the Gaussian probability density
```math
w(x) = \frac{1}{\sqrt{2\pi}} \mathrm{e}^{-x^2/2}
```
is symmetric relative to the origin $m=0$.
If the weight function has any parameters, then they are stored in the dictionary `pars`.
For example, the probability density of the Beta distribution on $\Omega = [0,1]$ has two positive shape parameters $\alpha, \beta > 0$
```math
w(x) = \frac{1}{B(\alpha,\beta)} x^{\alpha-1} (1-x)^{\beta-1}.
```
[This tutorial shows the above in action.](@ref UnivariateMonicOrthogonalPolynomials)

## OrthoPoly
Given an absolutely continuous measure we are wondering what are the monic polynomials $\phi_i: \Omega \rightarrow \mathbb{R}$ that are orthogonal relative to this very measure?
Mathematically this reads
```math
\langle \phi_i, \phi_j \rangle = \int_{\Omega} \phi_i(t) \phi_j(t) w(t) \mathrm{d}t =
\begin{cases}
> 0, & i=j \\
= 0, & i\neq j.
\end{cases}
```
They can be constructed using the type `OrthoPoly`, which has the fields

| Name | Meaning |
| --- | --- |
| `name::String`| Name|
| `deg::Int64`         | Maximum degree|
| `α::Vector{Float64}` | Vector of recurrence coefficients α|
| `β::Vector{Float64}` | Vector of recurrence coefficients β |
| `meas::Measure` | Underlying measure|


The purpose of `name` is obvious.
The integer `deg` stands for the maxium degree of the polynomials.
Rather than storing the polynomials $\phi_i$ themselves we store the recurrence coefficients `α`, `β` that characterize the system of orthogonal polynomials.
These recurrence coefficients are the single most important piece of information for the orthogonal polynomials.
For several common measures, there exist analytic formulae.
These are built-in to `PolyChaos` and should be used whenever possible.

[This tutorial shows the above in action.](@ref UnivariateMonicOrthogonalPolynomials)

## Quad
Quadrature rules are intricately related to orthogonal polynomials.
An $n$-point quadrature rule is a pair of so-called nodes $t_k$ and weights $w_k$ for $k=1,\dots,n$ that allow to solve integrals relative to the measure
```math
\int_\Omega f(t) w(t) \mathrm{d} t \approx \sum_{k=1}^n w_k f(t_k).
```
If the integrand $f$ is polynomial, then the specific Gauss quadrature rules possess the remarkable property that an $n$-point quadrature rule can integrate polynomial integrands $f$ of degree at most $2n-1$ *exactly*; no approximation error is made.

The fields of `Quad` are

| Name | Meaning |
| --- | --- |
| `name::String`        | Name |
| `Nquad::Int64`          | Number $n$ of quadrature points |
| `nodes::Vector{Float64}` | Nodes |
| `weights::Vector{Float64}` | Weights |
| `meas::Measure` | Underlying measure |

with obvious meanings.

[This tutorial shows the above in action.](@ref NumericalIntegration)

 ## OrthoPolyQ
As you would expect from the figure at the top, the type `OrthoPolyQ` is an amalgamation of `OrthoPoly` and `Quad`.
It has just those two fields

| Name | Meaning |
| --- | --- |
| `op::OrthoPoly`        | Orthogonal polynomials |
| `quad::Quad` | Quadrature rule |

Clearly, the underlying measures have to be the same.

[This tutorial shows the above in action.](@ref ComputationOfScalarProducts)

[Make sure to check out this tutorial too.](@ref NumericalIntegration)

## MultiMeasure
So far, everything was univariate, the weight of the measure was mapping real numbers to real numbers.
`PolyChaos` can handle product measures too.
Let's assume the weight function is a product of two independent Gaussian densities
```math
w: \mathbb{R} \times \mathbb{R} \rightarrow \mathbb{R}, \quad w(x) = \frac{1}{\sqrt{2\pi}} \mathrm{e}^{x_1^2/2} \frac{1}{\sqrt{2\pi}} \mathrm{e}^{x_2^2/2}.
```
The type `MultiMeasure` serves this purpose, with its fields

| Name | Meaning |
| --- | --- |
| `name::Vector{String}` | Name |
| `w::Function` | Weight function of product measure |
| `w_uni::Vector{Function}` | Weight functions of underlying univariate measures |
| `dom::Vector{Tuple{Float64,Float64}}` | Domain |
| `symmetric::Vector{Bool}` | Symmetry properties |
| `pars::Vector{Dict}` | Additioanl parameters |

All fields from `Measure` appear in vectorized versions (except for the weight $w$, which is the weight of the product measure)
The only *new* field is `w_uni`, which stacks the univariate weight functions.

[This tutorial shows the above in action.](@ref MultivariateMonicOrthogonalPolynomials)

## MultiOrthoPoly
Just as we did in the univariate case, we use `MultiMeasure` as a building block for multivariate orthogonal polynomials.
The type `MultiOrthoPoly` combines product measures with the respective orthogonal polynomials and their quadrature rules.
Its fields are

| Name | Meaning |
| --- | --- |
| `name::Vector{String}`| Names |
| `deg::Int64`| Maximum degree |
| `dim::Int64`| Dimension of basis |
| `ind::Matrix{Int64}`| Multi-index |
| `meas::MultiMeasure`| Underlying product measure |
| `uni::Union{Vector{OrthoPoly},Vector{OrthoPolyQ}}`| Underlying univariate orthogonal polynomials |

The names of the univariate bases are stored in `names`; the maximum degree of the basis is `deg`; the overall dimension of the multivariate basis is `dim`; the multi-index `ind` maps the $j$-th multivariate basis to the elements of the univariate bases; the product measure is stored in `meas`; finally, all univariate bases are collected in `uni`.

[This tutorial shows the above in action.](@ref MultivariateMonicOrthogonalPolynomials)

## Tensor
The last type we need to address is `Tensor`.
It is used to store the results of scalar products.
Its fields are

| Name | Meaning |
| --- | --- |
| `dim::Int64`| *Dimension* $m$ of tensor $\langle \phi_{i_1} \phi_{i_2} \cdots \phi_{i_{m-1}}, \phi_{i_m} \rangle$ |
| `T::SparseVector{Float64,Int64}`| Entries of tensor |
| `get::Function`| Function to get entries from `T` |
| `op::Union{OrthoPolyQ,MultiOrthoPoly}`| Underlying univariate orthogonal polynomials |

The *dimension* $m$ of the tensor is the number of terms that appear in the scalar product.
Let's assume we set $m = 3$, hence have $\langle \phi_{i} \phi_{j}, \phi_{k} \rangle$, then the concrete entry is obtained as `Tensor.get([i,j,k])`.

[This tutorial shows the above in action.](@ref ComputationOfScalarProducts)

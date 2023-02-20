```@setup mysetup
using PolyChaos
d = 6
myops = Dict()
polynames = ["gaussian", "beta01", "uniform01", "logistic"]
gaussian = GaussOrthoPoly(d);
myops["gaussian"] = gaussian
α, β = 1.3, 2.2
beta01 = Beta01OrthoPoly(d,α,β);
myops["beta01"] = beta01
uniform01 = Uniform01OrthoPoly(d);
myops["uniform01"] = uniform01
logistic = LogisticOrthoPoly(d);
myops["logistic"] = logistic;
myops
points, degrees = randn(10), 0:2:d
[ evaluate(degree,points,gaussian) for degree in degrees ]
μ, σ = 2., 0.2
pce_gaussian = convert2affinePCE(μ,σ,gaussian)
a, b = -0.3, 1.2
convert2affinePCE(a,b,uniform01)
pce_uniform = convert2affinePCE(μ,σ,uniform01;kind="μσ")
convert2affinePCE(a,b,beta01)
pce_beta = convert2affinePCE(μ,σ,beta01; kind="μσ")
a1, a2 = μ, sqrt(3)*σ/pi
pce_logistic = convert2affinePCE(a1,a2,logistic)
mean(pce_gaussian,myops["gaussian"]), std(pce_gaussian,myops["gaussian"])
mean(pce_uniform,myops["uniform01"]), std(pce_uniform,myops["uniform01"])
mean(pce_beta,myops["beta01"]), std(pce_beta,myops["beta01"])
mean(pce_logistic,myops["logistic"]), std(pce_logistic,myops["logistic"])
using Statistics
N = 1000
ξ_gaussian = sampleMeasure(N,myops["gaussian"])
samples_gaussian = evaluatePCE(pce_gaussian,ξ_gaussian,myops["gaussian"])
ξ_uniform = sampleMeasure(N,myops["uniform01"])
samples_uniform = evaluatePCE(pce_uniform,ξ_uniform,myops["uniform01"])
ξ_beta = sampleMeasure(N,myops["beta01"])
samples_beta = evaluatePCE(pce_beta,ξ_beta,myops["beta01"])
ξ_logistic = sampleMeasure(N,myops["logistic"])
samples_logistic = evaluatePCE(pce_logistic,ξ_logistic,myops["logistic"])

```

# [Common Random Variables](@id CommonRandomVariables)

Polynomial chaos expansion (PCE) is a Hilbert space technique for random variables with finite variance.
Mathematically equivalent to Fourier series expansions for periodic signals, PCE allows characterizing a random variable in terms of its PCE coefficients (aka Fourier coefficients).
That is, the PCE of a random variable $\mathsf{x}$ is given by

```math
\mathsf{x} = \sum_{i=0}^L x_i \phi_i,
```

where $x_i$ are the so-called PCE coefficients, and $\phi_i$ are the orthogonal polynomials that are orthogonal relative to the probability density function of $\mathsf{x}$.

This tutorial walks you through the PCE of common random variables, namely Gaussian (`gaussian`), Beta (`beta01`), Uniform(`uniform01`), Logistic (`logistic`), and shows how they are implemented in `PolyChaos`.

## Construction of Basis

```@example mysetup
using PolyChaos
```

The orthogonal polynomials are constructed using the `OrthoPoly`-type (here of degree at most `d`). For canonical measures, special constructors are implemented:

```@example mysetup
d = 6

myops = Dict()
polynames = ["gaussian", "beta01", "uniform01", "logistic"]

# gaussian
gaussian = GaussOrthoPoly(d);
myops["gaussian"] = gaussian

# beta01
α, β = 1.3, 2.2
beta01 = Beta01OrthoPoly(d, α, β);
myops["beta01"] = beta01

# uniform01
uniform01 = Uniform01OrthoPoly(d);
myops["uniform01"] = uniform01

# logistic
logistic = LogisticOrthoPoly(d);
myops["logistic"] = logistic;

myops
```

For example, let's evaluate the Gaussian basis polynomials at some points

```@example mysetup
points, degrees = randn(10), 0:2:d

[evaluate(degree, points, gaussian) for degree in degrees]
```

## Finding PCE Coefficients

Having constructed the orthogonal bases, the question remains how to find the PCE coefficients for the common random variables.
Every random variable can be characterized exactly by two PCE coefficients.
For a Gaussian random variable, this is familiar: the mean and the variance suffice to describe a Gaussian random variable entirely.
The same is true for any random variable of finite variance given the right basis.
The function `convert2affinePCE` provides the first two PCE coefficients (hence the name affine) for the common random variables.

### Gaussian

Given the Gaussian random variable $\mathsf{x} \sim \mathcal{N}(\mu, \sigma^2)$ with $\sigma > 0$, the affine PCE coefficients are

```@example mysetup
# Gaussian
μ, σ = 2.0, 0.2
pce_gaussian = convert2affinePCE(μ, σ, gaussian)
```

### Uniform

Given the uniform random variable $\mathsf{x} \sim \mathcal{U}(a, b)$ with finite support $a<b$, the affine PCE coefficients are

```@example mysetup
a, b = -0.3, 1.2
convert2affinePCE(a, b, uniform01)
```

Instead, if the expected value and standard deviation are known, the affine PCE coefficients of the uniform random variable are

```@example mysetup
pce_uniform = convert2affinePCE(μ, σ, uniform01; kind = "μσ")
# notice that the zero-order coefficient IS equivalent to the expected value μ
```

### Beta

Given the Beta random variable $\mathsf{x} \sim \mathcal{B}(a, b, \alpha, \beta)$ with finite support $a<b$ and shape parameters $\alpha, \beta > 0$, the affine PCE coefficients are

```@example mysetup
convert2affinePCE(a, b, beta01)
```

Instead, if the expected value and standard deviation are known, the affine PCE coefficients of the uniform random variable are

```@example mysetup
pce_beta = convert2affinePCE(μ, σ, beta01; kind = "μσ")
```

### Logistic

Given the logistic random variable $\mathsf{x} \sim \mathcal{L}(a_1,a_2)$ where $a_2>0$ with the probability density function

```math
\rho(t) = \frac{1}{4 a_2} \, \operatorname{sech}^2 \left(\frac{t-a_1}{2a_2}\right)
```

the affine PCE coefficients of the uniform random variable are

```@example mysetup
a1, a2 = μ, sqrt(3) * σ / pi
pce_logistic = convert2affinePCE(a1, a2, logistic)
```

## Moments

It is a key feature of PCE to compute moments from the PCE coefficients alone; no sampling is required.

### Gaussian

```@example mysetup
mean(pce_gaussian, myops["gaussian"]), std(pce_gaussian, myops["gaussian"])
```

### Uniform

```@example mysetup
mean(pce_uniform, myops["uniform01"]), std(pce_uniform, myops["uniform01"])
```

### Beta

```@example mysetup
mean(pce_beta, myops["beta01"]), std(pce_beta, myops["beta01"])
```

### Logistic

```@example mysetup
mean(pce_logistic, myops["logistic"]), std(pce_logistic, myops["logistic"])
```

## Sampling

Having found the PCE coefficients, it may be useful to sample the random variables.
That means, find $N$ realizations of the random variable that obey the random variable's probability density function.
This is done in two steps:

 1. Draw $N$ samples from the measure (`sampleMeasure()`), and then
 2. Evaluate the basis polynomials and multiply times the PCE coefficients, i.e. $\sum_{i=0}^L x_i \phi_i(\xi_j)$ where $\xi_j$ is the $j$-th sample from the measure (`evaluatePCE()`).

Both steps are combined in the function `samplePCE()`.

### Gaussian

```@example mysetup
using Statistics
N = 1000
ξ_gaussian = sampleMeasure(N, myops["gaussian"])
samples_gaussian = evaluatePCE(pce_gaussian, ξ_gaussian, myops["gaussian"])
# samplePCE(N,pce_gaussian,myops["gaussian"])
```

### Uniform

```@example mysetup
ξ_uniform = sampleMeasure(N, myops["uniform01"])
samples_uniform = evaluatePCE(pce_uniform, ξ_uniform, myops["uniform01"])
# samples_uniform = samplePCE(N,pce_uniform,myops["uniform01"])
```

### Beta

```@example mysetup
ξ_beta = sampleMeasure(N, myops["beta01"])
samples_beta = evaluatePCE(pce_beta, ξ_beta, myops["beta01"])
# samples_beta = samplePCE(N,pce_beta,myops["beta01"])
```

### Logistic

```@example mysetup
ξ_logistic = sampleMeasure(N, myops["logistic"])
samples_logistic = evaluatePCE(pce_logistic, ξ_logistic, myops["logistic"])
# samples_logistic = samplePCE(N,pce_logistic,myops["logistic"])
```

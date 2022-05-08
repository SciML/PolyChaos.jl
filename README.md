[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://polychaos.sciml.ai/stable)
[![Build Status](https://travis-ci.org/timueh/PolyChaos.jl.svg?branch=master)](https://travis-ci.org/timueh/PolyChaos.jl)
[![codecov](https://codecov.io/gh/timueh/PolyChaos.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/timueh/PolyChaos.jl)
[![Coverage Status](https://coveralls.io/repos/github/timueh/PolyChaos.jl/badge.svg?branch=master)](https://coveralls.io/github/timueh/PolyChaos.jl?branch=master)
[![Code DOI](https://zenodo.org/badge/165908711.svg)](https://zenodo.org/badge/latestdoi/165908711)
[![Paper@arXiv](https://img.shields.io/badge/arXiv-2004.03970-green.svg)](https://arxiv.org/abs/2004.03970)



# PolyChaos -- Orthogonal Polynomials, Quadrature, and Polynomial Chaos

A Julia package to construct orthogonal polynomials, their quadrature rules, and use it with polynomial chaos expansions.

The package requires `Julia 1.3` or newer.
In `Julia` switch to the package manager
```julia
julia> ]
(v1.3) pkg> add PolyChaos
```
This will install PolyChaos and its dependencies.
Once that is done, load the package:
```julia
julia> using PolyChaos
```
That's it.

Let's take a look at a simple example.
We would like to solve the integral

![equation](https://latex.codecogs.com/gif.latex?%5Cint_0%5E1%206%20x%5E5%20%5Cmathrm%7Bd%7Dx.)

Exploiting the underlying uniform measure, the integration can be done exactly with a 3-point quadrature rule.
```@example mysetup
opq = Uniform01OrthoPoly(3)
integrate(x -> 6x^5, opq)
```


For more information please visit the [documentation](https://polychaos.sciml.ai/stable).

## Citing

If you like `PolyChaos.jl`, consider citing our paper

```
@ARTICLE{2020arXiv200403970M,
       author = {{M{\"u}hlpfordt}, Tillmann and {Zahn}, Frederik and {Hagenmeyer}, Veit and {Faulwasser}, Timm},
        title = "{PolyChaos.jl -- A Julia Package for Polynomial Chaos in Systems and Control}",
      journal = {arXiv e-prints},
     keywords = {Electrical Engineering and Systems Science - Systems and Control, Mathematics - Numerical Analysis, Mathematics - Optimization and Control},
         year = 2020,
        month = apr,
          eid = {arXiv:2004.03970},
        pages = {arXiv:2004.03970},
archivePrefix = {arXiv},
       eprint = {2004.03970},
 primaryClass = {eess.SY},
}
```

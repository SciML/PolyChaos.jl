# [Computation of Scalar Products](@id ComputationOfScalarProducts)
By now, we are able to construct orthogonal polynomials, and to construct quadrature rules for a given nonnegative weight function, respectively.
Now we combine both ideas to solve integrals involving the orthogonal polynomials
```math
\langle \phi_{i_1} \phi_{i_2} \cdots \phi_{i_{m-1}}, \phi_{i_m} \rangle
= \int \phi_{i_1}(t) \phi_{i_2}(t) \cdots \phi_{i_{m-1}}(t) \phi_{i_m}(t) w(t) \mathrm{d} t,
```
both for the univariate and multivariate case.
The integrand is a polynomial (possibly multivariate) that can be solved exactly with the appropriate Gauss quadrature rules.

!!! note
    To simplify notation we drop the integration interval.
    It is clear from the context.


## Univariate Polynomials
### Classical Polynomials
Let's begin with a univariate basis for some *classical* orthogonal polynomial


```julia
using PolyChaos
deg, n = 4, 20
s_α, s_β = 2.1, 3.2
op = OrthoPoly("beta01",deg,Dict(:shape_a=>s_α,:shape_b=>s_β);Nrec=n)
```

To add the corresponding quadrature rule there is the composite struct `OrthoPolyQ` whose simplest constructor reads


```julia
opq = OrthoPolyQ(op,n)
```

By default, an $n$-point Gauss quadrature rule is create relative to the underlying measure `op.meas`, where $n$ is the number of recurrence coefficients stored in `op.α` and `op.β`.
The type `OrthoPolyQ` has just two fields: an `OrthoPoly`, and a `Quad`.

To compute the squared norms
```math
\| \phi_k \|^2 = \langle \phi_k, \phi_k  \rangle
= \int \phi_k(t) \phi_k(t) w(t) \mathrm{d} t
```

of the basis we call `computeSP2()`


```julia
normsq = computeSP2(opq)
```

For the general case
```math
\langle \phi_{i_1} \phi_{i_2} \cdots \phi_{i_{m-1}}, \phi_{i_m} \rangle
= \int \phi_{i_1}(t) \phi_{i_2}(t) \cdots \phi_{i_{m-1}}(t) \phi_{i_m}(t) w(t) \mathrm{d} t,
```
there exists a type `Tensor` that requires only two arguments: the *dimension* $m \geq 1$, and an `OrthoPolyQ`


```julia
m = 3
t = Tensor(3,opq)
```

To get the desired entries, `Tensor`comes with a `get()` function that is called for some index $a \in \mathbb{N}_0^m$ that has the entries $a = [i_1, i_2, \dots, i_m]$.
For example



```julia
t.get([1,2,3])
```

Or using comprehension


```julia
T = [ t.get([i1,i2,i3]) for i1=0:dim(opq)-1,i2=0:dim(opq)-1,i3=0:dim(opq)-1]
```

Notice that we can cross-check the results.


```julia
using LinearAlgebra
@show normsq == LinearAlgebra.diag(T[:,:,1])
@show normsq == LinearAlgebra.diag(T[:,1,:])
@show normsq == LinearAlgebra.diag(T[1,:,:])
```

Also, `normsq` can be computed analogously in `Tensor` format


```julia
t2 = Tensor(2,opq)
@show normsq == [ t2.get([i,i]) for i=0:dim(opq)-1]
```

### Arbitrary Weights
Of course, the type `OrthoPolyQ` can be constructed for arbitrary weights $w(t)$.
In this case we have to compute the orthogonal basis and the respective quadrature rule.
Let's re-work the above example by hand.


```julia
using SpecialFunctions
supp = (0,1)
function w(t)
    supp[1]<=t<=supp[2] ? (t^(s_α-1)*(1-t)^(s_β-1)/SpecialFunctions.beta(s_α,s_β)) : error("$t not in support")
end
my_meas = Measure("my_meas",w,supp,false,Dict())
my_op = OrthoPoly("my_op",deg,my_meas;Nrec=n)
my_quad = Quad(n,my_op)
my_opq = OrthoPolyQ(my_op,my_quad)
```

Now we can compute the squared norms $\| \phi_k \|^2$


```julia
my_normsq = computeSP2(my_opq)
```

And the tensor


```julia
my_t = Tensor(m,my_opq)
my_T = [ my_t.get([i1,i2,i3]) for i1=0:dim(opq)-1,i2=0:dim(opq)-1,i3=0:dim(opq)-1]
```

Let's compare the results:


```julia
@show abs.(normsq-my_normsq)
@show norm(T-my_T)
```

!!! note
    The possibility to create quadrature rules for arbitrary weights should be reserved to cases different from *classical* ones.

## Multivariate Polynomials
For multivariate polynomials the syntax for `Tensor` is very much alike, except that we are dealing with the type `MultiOrthoPoly` now.


```julia
mop = MultiOrthoPoly([opq,my_opq],deg)
```


```julia
mt2 = Tensor(2,mop)
mt3 = Tensor(3,mop)
mT2 = [ mt2.get([i,i]) for i=0:dim(mop)-1 ]
```

Notice that `mT2` carries the elements of the 2-dimensional tensors for the univariate bases `opq` and `my_opq`.
The encoding is given by the multi-index `mop.ind`


```julia
mop.ind
```

To cross-check the results we can distribute the multi-index back to its univariate indices with the help of `findUnivariateIndices`.


```julia
ind_opq = findUnivariateIndices(1,mop.ind)
ind_my_opq = findUnivariateIndices(2,mop.ind)
```


```julia
@show mT2[ind_opq] - normsq
@show mT2[ind_my_opq] - my_normsq;
```

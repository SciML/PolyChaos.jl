# Developer API

The following interfaces are intended for packages that extend PolyChaos. They
are versioned contracts, but ordinary applications should construct the concrete
measure, basis, quadrature, and tensor types documented in the user guide.

```@docs
AbstractMeasure
AbstractCanonicalMeasure
AbstractQuad
AbstractOrthoPoly
AbstractCanonicalOrthoPoly
AbstractTensor
```

## Generic Interface Guarantees

The generic methods are tested against independent subtypes rather than only
PolyChaos's concrete types. An `AbstractMeasure` subtype supplies `w`, `dom`,
and `symmetric`; an `AbstractQuad` subtype supplies paired `nodes` and
`weights`; and a univariate `AbstractOrthoPoly` subtype supplies `deg`, `α`,
`β`, `measure`, and `quad`. This supports `issymmetric`, `nw`, `integrate`,
`dim`, `deg`, `coeffs`, `evaluate`, `computeSP2`, and `Quad` without a concrete
PolyChaos basis type.

An `AbstractTensor` subtype supplies `dim`, sparse storage `T`, an index lookup
function `get`, and its basis `op`, which supports generic PCE variance
calculation.

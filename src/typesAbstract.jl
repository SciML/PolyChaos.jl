export AbstractMeasure,
    AbstractCanonicalMeasure,
    AbstractQuad,
    AbstractOrthoPoly,
    AbstractCanonicalOrthoPoly,
    AbstractTensor

"""
    AbstractQuad{T<:Real}

Developer extension interface for a quadrature rule with real node and weight
element type `T`.

# Required properties

Concrete subtypes must provide `nodes::AbstractVector{T}` and
`weights::AbstractVector{T}` fields of equal length. `nw(quad)` and
`integrate(f, quad)` use those fields directly. A rule with no nodes should use
[`EmptyQuad`](@ref) instead of a custom empty subtype.

# Developer API

This interface is versioned for PolyChaos extensions. It is not a general
quadrature interface: external code should use [`nw`](@ref) or
[`integrate`](@ref) rather than relying on concrete field layout.
"""
abstract type AbstractQuad{T <: Real} end

"""
    AbstractMeasure

Developer extension interface for a univariate measure used to define an
orthogonal-polynomial family.

# Required properties

Concrete subtypes must expose `w::Function`, `dom::Tuple{<:Real,<:Real}`, and
`symmetric::Bool`. The weight `w(t)` must be nonnegative on `dom`; `symmetric`
must be true only when the measure is symmetric about zero. Generic
implementations of [`issymmetric`](@ref), [`Quad`](@ref), and [`OrthoPoly`](@ref)
read these properties.

# Developer API

Use [`Measure`](@ref) for an ordinary user-defined measure. Subtyping is for
package developers that also provide an appropriate recurrence or
discretization path.
"""
abstract type AbstractMeasure end

"""
    AbstractCanonicalMeasure <: AbstractMeasure

Developer extension interface for a measure with a known canonical
orthogonal-polynomial family.

# Rules

Subtypes satisfy the [`AbstractMeasure`](@ref) contract and should be paired
with an `OrthoPoly(::YourMeasure, deg; Nrec, addQuadrature)` method. This lets
the generic `OrthoPoly(measure, deg)` entry point construct the canonical
basis without numerical discretization.
"""
abstract type AbstractCanonicalMeasure <: AbstractMeasure end

"""
    AbstractOrthoPoly{M<:AbstractMeasure,Q<:AbstractQuad}

Developer extension interface for a univariate or multivariate
orthogonal-polynomial basis.

# Required properties

Univariate subtypes must provide `deg::Int`, recurrence coefficient vectors
`α` and `β`, `measure::M`, and `quad::Q`. `length(α) == length(β) >= deg + 1`
is required. [`dim`](@ref), [`deg`](@ref), [`coeffs`](@ref), [`nw`](@ref),
[`evaluate`](@ref), and [`integrate`](@ref) use this contract.

Multivariate subtypes additionally provide `uni`, `ind`, and `dim`; see
[`MultiOrthoPoly`](@ref) for the concrete layout.

# Developer API

This is a versioned extension interface. User code should construct one of the
documented concrete bases rather than subtype it directly.
"""
abstract type AbstractOrthoPoly{M <: AbstractMeasure, Q <: AbstractQuad} end

"""
    AbstractCanonicalOrthoPoly{V,M,Q} <: AbstractOrthoPoly{M,Q}

Developer extension interface for a basis associated with an
[`AbstractCanonicalMeasure`](@ref).

# Rules

Subtypes follow the [`AbstractOrthoPoly`](@ref) field contract. `V` is the
concrete vector type used for both recurrence coefficient vectors. Canonical
bases should retain enough recurrence coefficients to satisfy `Nrec >= deg +
1` and use [`EmptyQuad`](@ref) only when `addQuadrature = false` was requested.
"""
abstract type AbstractCanonicalOrthoPoly{V <: AbstractVector{<:Real}, M, Q} <:
AbstractOrthoPoly{M, Q} end

"""
    AbstractTensor{T<:AbstractOrthoPoly}

Developer extension interface for sparse scalar-product tensors associated
with an orthogonal-polynomial basis.

# Required properties

Concrete subtypes must provide `dim::Int`, `T::SparseVector`, `get::Function`,
and `op::T`. The `get` function accepts a basis-index tuple and returns the
corresponding scalar product. Generic variance calculations read `T` directly.

# Developer API

Use [`Tensor`](@ref) for normal construction. Subtyping is reserved for
extensions that preserve the stated storage contract.
"""
abstract type AbstractTensor{T <: AbstractOrthoPoly} end

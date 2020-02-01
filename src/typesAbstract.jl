export  AbstractMeasure,
        AbstractCanonicalMeasure,
        AbstractQuad,
        AbstractOrthoPoly,
        AbstractCanonicalOrthoPoly,
        AbstractTensor


abstract type AbstractQuad{T<:Real} end

abstract type AbstractMeasure end
abstract type AbstractCanonicalMeasure <: AbstractMeasure end

abstract type AbstractOrthoPoly{M<:AbstractMeasure, Q<:AbstractQuad} end
abstract type AbstractCanonicalOrthoPoly{V<:AbstractVector{<:Real}, M, Q} <: AbstractOrthoPoly{M, Q} end

abstract type AbstractTensor{T<:AbstractOrthoPoly} end

export  AbstractMeasure,
        AbstractCanonicalMeasure,
        AbstractQuad,
        AbstractOrthoPoly,
        AbstractCanonicalOrthoPoly,
        AbstractTensor


abstract type AbstractQuad end

abstract type AbstractMeasure end
abstract type AbstractCanonicalMeasure <: AbstractMeasure end

abstract type AbstractOrthoPoly end
abstract type AbstractCanonicalOrthoPoly <: AbstractOrthoPoly end

abstract type AbstractTensor end

using LinearAlgebra, PolyChaos, SparseArrays, Statistics, Test

struct GenericMeasure <: AbstractMeasure
    w::Function
    dom::Tuple{Float64, Float64}
    symmetric::Bool
end

struct GenericQuad <: AbstractQuad{Float64}
    nodes::Vector{Float64}
    weights::Vector{Float64}
end

struct GenericBasis <: AbstractOrthoPoly{GenericMeasure, GenericQuad}
    deg::Int
    α::Vector{Float64}
    β::Vector{Float64}
    measure::GenericMeasure
    quad::GenericQuad
end

struct GenericTensor <: AbstractTensor{GenericBasis}
    dim::Int
    T::SparseVector{Float64, Int}
    get::Function
    op::GenericBasis
end

@testset "Developer abstract interfaces" begin
    measure = GenericMeasure(x -> 0.5, (-1.0, 1.0), true)
    quad = GenericQuad([-1.0, 1.0], [0.5, 0.5])
    basis = GenericBasis(2, [0.0, 0.0, 0.0], [1.0, 1 / 3, 4 / 15], measure, quad)

    @test PolyChaos.issymmetric(measure)
    @test PolyChaos.issymmetric !== LinearAlgebra.issymmetric
    @test nw(quad) == [-1.0 0.5; 1.0 0.5]
    @test integrate(x -> x^2, quad) == 1.0
    @test Quad(2, measure) isa Quad

    @test deg(basis) == 2
    @test dim(basis) == 3
    @test coeffs(basis) == [0.0 1.0; 0.0 1 / 3; 0.0 4 / 15]
    @test nw(basis) == nw(quad)
    @test evaluate(1, 0.25, basis) == 0.25
    @test computeSP([1, 1], basis) == 1 / 3
    @test computeSP2(basis) == [basis.β[1], basis.β[1] * basis.β[2], prod(basis.β)]
    @test Quad(2, basis) isa Quad

    tensor = GenericTensor(
        2,
        sparsevec([1], [1.0], 3),
        index -> index == [1, 1] ? 1 / 3 : 1.0,
        basis,
    )
    @test PolyChaos.var([1.0, 2.0], tensor) == 4 / 3
    @test PolyChaos.mean !== Statistics.mean
    @test PolyChaos.var !== Statistics.var
    @test PolyChaos.std !== Statistics.std
    @test PolyChaos.mean([1.0, 3.0]) == 2.0
    @test PolyChaos.var([1.0, 3.0]) == 2.0
    @test PolyChaos.std([1.0, 3.0]) == sqrt(2.0)
end

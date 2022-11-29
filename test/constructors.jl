using PolyChaos, Test

deg = 4
Nrec = 10

function compare_quad(A::Quad, B::Quad)
    (A.Nquad == B.Nquad) && (A.name == B.name) && (A.nodes == B.nodes) &&
        (A.weights == B.weights)
end

function compare_quad(A::EmptyQuad, B::EmptyQuad)
    true
end

function compare_orthopoly(A, B)
    (A.deg == B.deg) && compare_quad(A.quad, B.quad) && (A.measure == B.measure) &&
        (A.α == B.α) && (A.β == B.β)
end

@testset "LegendreOrthoPoly constructors" begin
    meas = LegendreMeasure()
    for quad in [true, false]
        polyA = LegendreOrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "JacobiOrthoPoly constructors" begin
    shape_a, shape_b = rand(), rand()
    meas = JacobiMeasure(shape_a, shape_b)
    for quad in [true, false]
        polyA = JacobiOrthoPoly(deg, shape_a, shape_b, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "LaguerreOrthoPoly constructors" begin
    meas = LaguerreMeasure()
    for quad in [true, false]
        polyA = LaguerreOrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "genLaguerreOrthoPoly constructors" begin
    shape = rand()
    meas = genLaguerreMeasure(shape)
    for quad in [true, false]
        polyA = genLaguerreOrthoPoly(deg, shape, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "HermiteOrthoPoly constructors" begin
    meas = HermiteMeasure()
    for quad in [true, false]
        polyA = HermiteOrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "genHermiteOrthoPoly constructors" begin
    mu = rand()
    meas = genHermiteMeasure(mu)
    for quad in [true, false]
        polyA = genHermiteOrthoPoly(deg, mu, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "MeixnerPollaczekOrthoPoly constructors" begin
    lambda, phi = rand(), rand()
    meas = MeixnerPollaczekMeasure(lambda, phi)
    for quad in [true, false]
        polyA = MeixnerPollaczekOrthoPoly(deg, lambda, phi, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "GaussOrthoPoly constructors" begin
    meas = GaussMeasure()
    for quad in [true, false]
        polyA = GaussOrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "Uniform01OrthoPoly constructors" begin
    meas = Uniform01Measure()
    for quad in [true, false]
        polyA = Uniform01OrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "Uniform_11OrthoPoly constructors" begin
    meas = Uniform_11Measure()
    for quad in [true, false]
        polyA = Uniform_11OrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "Beta01OrthoPoly constructors" begin
    shape_a, shape_b = rand(), rand()
    meas = Beta01Measure(shape_a, shape_b)
    for quad in [true, false]
        polyA = Beta01OrthoPoly(deg, shape_a, shape_b, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "GammaOrthoPoly constructors" begin
    shape, rate = rand(), 1.0
    meas = GammaMeasure(shape, rate)
    for quad in [true, false]
        polyA = GammaOrthoPoly(deg, shape, rate, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

@testset "LogisticOrthoPoly constructors" begin
    meas = LogisticMeasure()
    for quad in [true, false]
        polyA = LogisticOrthoPoly(deg, Nrec = 10, addQuadrature = quad)
        polyB = OrthoPoly(meas, deg; Nrec = 10, addQuadrature = quad)
        @test compare_orthopoly(polyA, polyB)
    end
end

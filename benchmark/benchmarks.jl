using PolyChaos, BenchmarkTools
using StableRNGs

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

deg = 8
op = Uniform01OrthoPoly(deg)
op_g = GaussOrthoPoly(deg)
op_l = genLaguerreOrthoPoly(deg, 1.2)

x = 0.4
xs = rand(rng, 500)
inds = collect(0:deg)
cf = rand(rng, deg + 1)

# =============================================================================
# Orthogonal polynomial construction (three-term recurrence coeffs)
# =============================================================================

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["uniform01"] = @benchmarkable Uniform01OrthoPoly($deg)
SUITE["construct"]["gauss"] = @benchmarkable GaussOrthoPoly($deg)
SUITE["construct"]["laguerre"] = @benchmarkable genLaguerreOrthoPoly(
    $deg, 1.2
)
SUITE["construct"]["hermite"] = @benchmarkable GaussOrthoPoly(
    $deg; addQuadrature = false
)

# =============================================================================
# Evaluation / quadrature
# =============================================================================

SUITE["evaluate"] = @benchmarkable evaluate($inds, $x, $op)
SUITE["evaluate_batch"] = @benchmarkable evaluate($inds, $xs, $op)

SUITE["quad"] = BenchmarkGroup()
SUITE["quad"]["fejer"] = @benchmarkable fejer(20)
SUITE["quad"]["integrate"] = @benchmarkable integrate(x -> sin(π * x), $op)
SUITE["quad"]["computeSP"] = @benchmarkable computeSP([1, 2], $op)

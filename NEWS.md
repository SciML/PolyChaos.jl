# PolyChaos Release Notes

## Version 0.2.1
- increased code coverage
- examples work with new type hierarchy
- updated dependencies

## Version 0.2.0
- new type hierarchy that involves four abstract types:
    - `AbstractMeasure`
    - `AbstractOrthoPoly`
    - `AbstractQuad`
    - `AbstractTensor`
- canonical measures and orthogonal polynomials are now subtypes of `AbstractCanonicalMeasure` and `AbstractOrthoPoly`
- code hygiene at several fronts

## Version 0.1.3
- removed default `show` function for types specific to `PolyChaos.jl` as this was causing issues with the Atom editor
- improved `showpoly`, and added `showbasis`
- updated documentation
- updated package dependencies
- added example for chance-constrained DC optimal power flow to documentation [see 004be7c](https://github.com/timueh/PolyChaos.jl/commit/004be7c4581ce035dc033da3810f4266555d9206)
- computation of moments from PCE coefficients is now compatible with `JuMP` (was restricted to `Vector{Float64}` before), see [see f99030](https://github.com/timueh/PolyChaos.jl/commit/f9903093e6556c4f37dabd82f03595323b59cc96)

## Version 0.1.2 (Mar 13, 2019)
- implemented show function for orthogonal polynomials ([thanks @pfitzseb](https://discourse.julialang.org/t/how-to-define-verbose-output-for-a-polynomial/21317/5) )
- improved Golub-Welsch algorithm ([adapted from `GaussQuadrature.jl`](https://github.com/billmclean/GaussQuadrature.jl/blob/master/src/GaussQuadrature.jl) which itself is based on [this Fortran code](https://www.netlib.org/cgi-bin/netlibfiles.pl?filename=/go/gaussq.f))
- improved code coverage, i.e. added tests
- code hygiene

## Version 0.1.1 (Feb 27, 2019)
- improved performance of recurrence coefficients, quadrature rules, evaluation of orthogonal polynomials
- added considerable number of tests for recurrence coefficients (see #4) and quadrature rules (see #5)
- added code coverage
- added GNU General Public License v3.0
- removed dependency on `IterTools`

## Version 0.1.0 (Feb 18, 2019)
- initial version

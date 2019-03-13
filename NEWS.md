# PolyChaos Release Notes

## Version 0.1.2
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

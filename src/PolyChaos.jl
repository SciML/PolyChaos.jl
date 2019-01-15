module PolyChaos

using SpecialFunctions, FastGaussQuadrature, SparseArrays, Distributions
import LinearAlgebra: I, dot, SymTridiagonal, eigen, issymmetric
import FFTW: ifft
import Combinatorics: with_replacement_combinations
import Base: show
import AdaptiveRejectionSampling: RejectionSampler, run_sampler!
import Statistics: mean, std, var

export

calculateMultiIndices,

r_scale,

rm_logistic, rm_laguerre, rm_logisticsum, rm_hermite, rm_hermite_prob, rm_laguerre, rm_jacobi,
rm_jacobi01, rm_hahn, rm_meixner_pollaczek, rm_legendre, rm_legendre01, rm_compute, clenshaw_curtis,

OrthoPoly, MultiOrthoPoly, Quad, OrthoPolyQ, Measure, MultiMeasure, Tensor,

stieltjes, lanczos, golubwelsch, quadpts,

evaluate,

calculateMultiIndices,

w_gaussian, w_uniform01, w_logistic, build_w_beta, build_w_gamma,
build_w_hermite,

coeffs, nw,

quadgp, fejer, fejer2, clenshaw_curtis, quadpts_beta01, quadpts_gamma, quadpts_gaussian,
quadpts_logistic, quadpts_uniform01,

computeSP2, computeSP, computeTensorizedSP,

integrate, issymmetric,

multi2uni,

getentry, findUnivariateIndices,

dim, deg

# add local functions / scripts
include("types.jl")
include("multiIndices.jl")
include("evaluate.jl")
include("quadrature_rules.jl")
include("recurrence_coefficients_monic.jl")
include("stieltjes.jl")
include("densities.jl")
include("auxfuns.jl")
include("scalarproduct.jl")
include("tensor.jl")
include("show.jl")
include("polynomial_chaos.jl")


end

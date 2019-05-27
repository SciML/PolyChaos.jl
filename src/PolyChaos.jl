module PolyChaos

using SpecialFunctions, SparseArrays, Distributions
import LinearAlgebra: I, dot, SymTridiagonal, eigen, issymmetric, pinv
import FFTW: ifft
import Combinatorics: with_replacement_combinations
import Base: show
import AdaptiveRejectionSampling: RejectionSampler, run_sampler!
import Statistics: mean, std, var
import GaussQuadrature: special_eigenproblem!

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
include("sparse_basis.jl")

end

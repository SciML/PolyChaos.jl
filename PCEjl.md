[_back to main.md_](../main.md)
---
## Roadmap for ``PCE.jl``
- The toolbox is in a very good state. It's stable and I was able to add stuff to it in a decent amount of time.
- current drawbacks:
  - type for every uncertainty
  - code does not represent the mathematical picture behind:
    - every monic orthogonal basis is isomorphic to its recurrence relation coefficients
    - recurrence relation is the heart of orthogonal polynomials (which it isn't in our package)
    - no easy addition of features from Gautschi's code snippets possible
  - can we thin out the code to increase performance?
  - can we reduce number of dependencies?
  - with type ``Basis_t`` everything is computed which appears to make the code slow; often times we only need the scalar products. E.g. the basis -- if one is interested in the basis, then there should be a command for this (``print_basis(PCE_type)``)
  - computation of Galerkin tensor → it's computationally extremely demanding; can we define some clever look-up table/hash table?
  - plotting done using PyPlot which is slow; can we switch to ``Plots.jl`` and use recipes?
  - no testing done for ``Julia 0.7``/``Julia 1.0``
  - monic vs "classic" orthogonal polynomials that are neither monic nor orthonormal
- possible collaborators/people we could reach out to:
  - Alex Townsend: the guy who's behind ``FastGaussQuadrature`` (and also involved with ``chebfun``)

# Type structure
  ```julia
  struct OrthoPoly             # =univariate MONIC orthogonal polynomial
    name::String
    deg::Int64          # maximum degree
    w::Function
    supp::Tuple{Float64,Float64}
    α::Vector{Float64}  # recurrence coefficients
    β::Vector{Float64}  # recurrence coefficients

    # --> add constructors 'canonical' uncertainties
    if name=="Gaussian" then a,b =
    elseif name=="Beta" then a,b =..
    else
      r_compute()
    end
  end

  # Struct that contains pre-computed nodes and weights
  struct OrthoPolyQ
    op::OrthoPoly
    nodes::Vector{Float64}
    weights::Vector{Float64}

    # the constructor should call the function quadpts()
    # two parameters are important:
    Nquad::Int => number of nodes, weights
    Npoly::Int => number of rec. coefficients used to computed nodes, weights
  end

  struct MultiOrthoPoly
    name::Vector{String}
    deg::Int64
    w::Function
    supp::Vector{Tuple{Float64,Float64}}
    ind::Matrix{Float64} # multi-index
    uni::Vector{OrthoPolyQ,OrthoPoly} # vector of orthogonal polynomials
  end

  struct Tensor
    dim::Int64          # "dimension"
    T::spzeros()
    op::Union{OrthoPoly,MultiOrthoPoly}

    # constructor has to take into account op
  end

  function get(a::Vector{Int64},T::Tensor)
    @assert length(a)==T.dim
    # ...
  end

  function quadpts(op::OrthoPoly;Nquad::Int,Npoly::Int)
    # computes the nodes and weights for an orthogonal polynomial system
    # Nquad::Int => number of nodes, weights
    # Npoly::Int => number of rec. coefficients used to computed nodes,
    if op.name=="Gaussian" gauss_hermite(Nquad)
    elseif op.name=="Beta" gauss_jacobi(Nquad)
    else
      golub_welsch()
    end
      return nodes, weights
  end

  function evaluate_sp(a::Vector{Int64},op::OrthoPoly)
    this code exists
  end

  function evaluate(op::OrthoPoly,x::Float64)
  end

  function evaluate(mop::MultiOrthoPoly,x::Vector{Float64})
    # use multi-index and then evaluate every univariate OrthoPoly
    # according to evaluate(op::OrthoPoly,x::Float64)

  function integrate(f::Function,nodes,weights)
    returns dot(weights,f.(nodes))
  end
  function integrate(f::Function,op::OrthogonalPolyQ)
    returns integrate(f,op.nodes,op.weights)
  end
  # integrate has to be dispatched for all special cases, e.g. integrating basis polynomials for Hermites and so on...

  function orthonormalize(op::OrthoPoly)
    # returns orthonormalized OrthoPoly
    # need to scale recursion coefficients
  end

  ```

# Meeting Nov 29, 2018
- ``DynamicPolynomials`` and ``Polynomials`` should be removed
- to display the basis we use text formatting
- Re-write the tensor to avoid the definition of large arrays --> it's all there in ``PCE.jl``, but it's not returned

---
[_back to main.md_](../main.md)

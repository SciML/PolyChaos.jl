var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": "using PolyChaos"
},

{
    "location": "#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "PolyChaos is a collection of numerical routines for orthogonal polynomials written in the Julia programming language. Starting from some non-negative weight (aka an absolutely continuous nonnegative measure), PolyChaos allowsto compute the coefficients for the monic three-term recurrence relation,\nto evaluate the orthogonal polynomials at arbitrary points,\nto compute the quadrature rule,\nto compute tensors of scalar products,\nto do all of the above in a multivariate setting (aka product measures).If the weight function is a probability density function, PolyChaos further provides routines to compute polynomial chaos expansions (PCEs) of random variables with this very density function. These routines allowto compute affine PCE coefficients for arbitrary densities,\nto compute moments,\nto compute the tensors of scalar products.PolyChaos contains several canonical orthogonal polynomials such as Jacobi or Hermite polynomials. For these, closed-form expressions and state-of-the art quadrature rules are used whenever possible. However, a cornerstone principle of PolyChaos is to provide all the functionality for user-specific, arbitrary weights.note: Note\nWhat PolyChaos is not (at least currently):a self-contained introduction to orthogonal polynomials, quadrature rules and/or polynomial chaos expansions. We assume the user brings some experience to the table. However, over time we will focus on strengthening the tutorial charater of the package.\na symbolic toolbox\na replacement for FastGaussQuadrature.jl"
},

{
    "location": "#Installation-1",
    "page": "Overview",
    "title": "Installation",
    "category": "section",
    "text": "The package requires Julia 1.0 or newer. In Julia switch to the package managerjulia> ]\n(v1.0) pkg> add PolyChaosThis will install PolyChaos and its dependencies. Once that is done, load the package:julia> using PolyChaosThat\'s it.Let\'s take a look at a simple example. We would like to solve the integralint_0^1 6 x^5 mathrmdxExploiting the underlying uniform measure, the integration can be done exactly with a 3-point quadrature rule.opq = OrthoPolyQ(\"uniform01\",3)\nintegrate(x->6x^5,opq)To get going with PolyChaos check out the tutorials such as the one on numerical integration. In case you are unfamiliar with orthogonal polynomials, perhaps this background information is of help."
},

{
    "location": "#References-1",
    "page": "Overview",
    "title": "References",
    "category": "section",
    "text": "The code base of PolyChaos is partially based on Walter Gautschi\'s Matlab suite of programs for generating orthogonal polynomials and related quadrature rules, with much of the theory presented in his book Orthogonal Polynomials: Computation and Approximation published in 2004 by the Oxford University Press.For the theory of polynomial chaos expansion we mainly consulted T. J. Sullivan. Introduction to Uncertainty Quantification. Springer International Publishing Switzerland. 2015."
},

{
    "location": "#Contributing-1",
    "page": "Overview",
    "title": "Contributing",
    "category": "section",
    "text": "We are always looking for contributors. If you are interested, just get in touch: tillmann [dot] muehlpfordt [at] kit [dot] edu.Or just fork and/or star the repository:Julia\'s package manager works nicely with Github: simply install the hosted package via Pkg.clone and the repository\'s URL. A so-called fork is created withPkg.clone(\"https://github.com/timueh/PolyChaos.jl\")The fork will replace the original package.CallPkg.dir(\"PolyChaos\")to figure out where the package was cloned to. Go to that location and figure out what branch you are on via git branch."
},

{
    "location": "#Citing-1",
    "page": "Overview",
    "title": "Citing",
    "category": "section",
    "text": "Currently, there is no publication about PolyChaos. Meanwhile, in case you find PolyChaos useful, feel free to get in touch, or simply participate in Github\'s gamification. ;)"
},

{
    "location": "type_hierarchy/#",
    "page": "Type Hierarchy",
    "title": "Type Hierarchy",
    "category": "page",
    "text": ""
},

{
    "location": "type_hierarchy/#Type-Hierarchy-1",
    "page": "Type Hierarchy",
    "title": "Type Hierarchy",
    "category": "section",
    "text": "Let\'s look at the types PolyChaos provides. The high-level perspective looks as such: (Image: type hierarchy)The left hand side covers types related to univariate measures; the right hand side covers multivariate measures. An arrow beginning at one type and ending at another type means that the beginning type is a field of the ending type. For example, the type OrthoPoly has a field of type Measure; the type OrthoPolyQ has a field of type OrthoPoly and a field of type Quad, and so on. Let\'s begin with the univariate case.note: Note\nIf you are unfamiliar with the mathematical background of orthogonal polynomials, please consult this tutorial."
},

{
    "location": "type_hierarchy/#Measure-1",
    "page": "Type Hierarchy",
    "title": "Measure",
    "category": "section",
    "text": "It all begins with a measure, more specifically absolutely continuous measures. What are the fields of such a type measure?Field Meaning\nname::String Name of measure\nw::Function Weight function w Omega rightarrow mathbbR\ndom::Tuple{Float64,Float64} Domain $ \\Omega$\nsymmetric::Bool Is w symmetric relative to some m in Omega, hence w(m-x) = w(m+x) for all x in Omega?\npars::Dict Additional parameters (e.g. shape parameters for Beta distributionThey are a name, a weight function w Omega rightarrow mathbbR with domain Omega (dom). If the weight function is symmetric relative to some m in Omega, the field symmetric should be set to true. Symmetry relative to m means thatforall x in Omega quad w(m-x) = w(m+x)For example, the Gaussian probability densityw(x) = frac1sqrt2pi mathrme^-x^22is symmetric relative to the origin m=0. If the weight function has any parameters, then they are stored in the dictionary pars. For example, the probability density of the Beta distribution on Omega = 01 has two positive shape parameters alpha beta  0w(x) = frac1B(alphabeta) x^alpha-1 (1-x)^beta-1This tutorial shows the above in action."
},

{
    "location": "type_hierarchy/#OrthoPoly-1",
    "page": "Type Hierarchy",
    "title": "OrthoPoly",
    "category": "section",
    "text": "Given an absolutely continuous measure we are wondering what are the monic polynomials phi_i Omega rightarrow mathbbR that are orthogonal relative to this very measure? Mathematically this readslangle phi_i phi_j rangle = int_Omega phi_i(t) phi_j(t) w(t) mathrmdt =\nbegincases\n 0  i=j \n= 0  ineq j\nendcasesThey can be constructed using the type OrthoPoly, which has the fieldsName Meaning\nname::String Name\ndeg::Int64 Maximum degree\nα::Vector{Float64} Vector of recurrence coefficients α\nβ::Vector{Float64} Vector of recurrence coefficients β\nmeas::Measure Underlying measureThe purpose of name is obvious. The integer deg stands for the maxium degree of the polynomials. Rather than storing the polynomials phi_i themselves we store the recurrence coefficients α, β that characterize the system of orthogonal polynomials. These recurrence coefficients are the single most important piece of information for the orthogonal polynomials. For several common measures, there exist analytic formulae. These are built-in to PolyChaos and should be used whenever possible.This tutorial shows the above in action."
},

{
    "location": "type_hierarchy/#Quad-1",
    "page": "Type Hierarchy",
    "title": "Quad",
    "category": "section",
    "text": "Quadrature rules are intricately related to orthogonal polynomials. An n-point quadrature rule is a pair of so-called nodes t_k and weights w_k for k=1dotsn that allow to solve integrals relative to the measureint_Omega f(t) w(t) mathrmd t approx sum_k=1^n w_k f(t_k)If the integrand f is polynomial, then the specific Gauss quadrature rules possess the remarkable property that an n-point quadrature rule can integrate polynomial integrands f of degree at most 2n-1 exactly; no approximation error is made.The fields of Quad areName Meaning\nname::String Name\nNquad::Int64 Number n of quadrature points\nnodes::Vector{Float64} Nodes\nweights::Vector{Float64} Weights\nmeas::Measure Underlying measurewith obvious meanings.This tutorial shows the above in action."
},

{
    "location": "type_hierarchy/#OrthoPolyQ-1",
    "page": "Type Hierarchy",
    "title": "OrthoPolyQ",
    "category": "section",
    "text": "As you would expect from the figure at the top, the type OrthoPolyQ is an amalgamation of OrthoPoly and Quad. It has just those two fieldsName Meaning\nop::OrthoPoly Orthogonal polynomials\nquad::Quad Quadrature ruleClearly, the underlying measures have to be the same.This tutorial shows the above in action.Make sure to check out this tutorial too."
},

{
    "location": "type_hierarchy/#MultiMeasure-1",
    "page": "Type Hierarchy",
    "title": "MultiMeasure",
    "category": "section",
    "text": "So far, everything was univariate, the weight of the measure was mapping real numbers to real numbers. PolyChaos can handle product measures too. Let\'s assume the weight function is a product of two independent Gaussian densitiesw mathbbR times mathbbR rightarrow mathbbR quad w(x) = frac1sqrt2pi mathrme^x_1^22 frac1sqrt2pi mathrme^x_2^22The type MultiMeasure serves this purpose, with its fieldsName Meaning\nname::Vector{String} Name\nw::Function Weight function of product measure\nw_uni::Vector{Function} Weight functions of underlying univariate measures\ndom::Vector{Tuple{Float64,Float64}} Domain\nsymmetric::Vector{Bool} Symmetry properties\npars::Vector{Dict} Additioanl parametersAll fields from Measure appear in vectorized versions (except for the weight w, which is the weight of the product measure) The only new field is w_uni, which stacks the univariate weight functions.This tutorial shows the above in action."
},

{
    "location": "type_hierarchy/#MultiOrthoPoly-1",
    "page": "Type Hierarchy",
    "title": "MultiOrthoPoly",
    "category": "section",
    "text": "Just as we did in the univariate case, we use MultiMeasure as a building block for multivariate orthogonal polynomials. The type MultiOrthoPoly combines product measures with the respective orthogonal polynomials and their quadrature rules. Its fields areName Meaning\nname::Vector{String} Names\ndeg::Int64 Maximum degree\ndim::Int64 Dimension of basis\nind::Matrix{Int64} Multi-index\nmeas::MultiMeasure Underlying product measure\nuni::Union{Vector{OrthoPoly},Vector{OrthoPolyQ}} Underlying univariate orthogonal polynomialsThe names of the univariate bases are stored in names; the maximum degree of the basis is deg; the overall dimension of the multivariate basis is dim; the multi-index ind maps the j-th multivariate basis to the elements of the univariate bases; the product measure is stored in meas; finally, all univariate bases are collected in uni.This tutorial shows the above in action."
},

{
    "location": "type_hierarchy/#Tensor-1",
    "page": "Type Hierarchy",
    "title": "Tensor",
    "category": "section",
    "text": "The last type we need to address is Tensor. It is used to store the results of scalar products. Its fields areName Meaning\ndim::Int64 Dimension m of tensor langle phi_i_1 phi_i_2 cdots phi_i_m-1 phi_i_m rangle\nT::SparseVector{Float64,Int64} Entries of tensor\nget::Function Function to get entries from T\nop::Union{OrthoPolyQ,MultiOrthoPoly} Underlying univariate orthogonal polynomialsThe dimension m of the tensor is the number of terms that appear in the scalar product. Let\'s assume we set m = 3, hence have langle phi_i phi_j phi_k rangle, then the concrete entry is obtained as Tensor.get([i,j,k]).This tutorial shows the above in action."
},

{
    "location": "numerical_integration/#",
    "page": "Numerical Integration",
    "title": "Numerical Integration",
    "category": "page",
    "text": ""
},

{
    "location": "numerical_integration/#NumericalIntegration-1",
    "page": "Numerical Integration",
    "title": "Numerical Integration",
    "category": "section",
    "text": "using PolyChaos\nn = 5; f(t) = sin(t)\nopq = OrthoPolyQ(\"uniform01\",n-1);\nI0 = integrate(f,opq);\nm = Measure(\"uniform01\");\nq = Quad(n-1,m);\nI1 = integrate(f,q)\nop = OrthoPoly(\"uniform01\",n-1)\nq = Quad(n,op)\nI2 = integrate(f,q)The goal of this tutorial is to solve an integral using Gauss quadrature,I = int_0^1 f(t) mathrmd t approx sum_k=1^n w_k f(t_k)where we choose f(t) = sin(t), and n = 5.Make sure to check out this tutorial too."
},

{
    "location": "numerical_integration/#Variant-0-1",
    "page": "Numerical Integration",
    "title": "Variant 0",
    "category": "section",
    "text": "using PolyChaos\nn = 5;\nf(t) = sin(t);\nopq = OrthoPolyQ(\"uniform01\",n-1);\nI0 = integrate(f,opq)\nprint(\"Numerical error: $(abs(1-cos(1)-I0))\")with negligible numerical errors."
},

{
    "location": "numerical_integration/#Variant-1-1",
    "page": "Numerical Integration",
    "title": "Variant 1",
    "category": "section",
    "text": "Let us  now solve the same problem, while elaborating what is going on under the hood. At first, we load the package by callingusing PolyChaosNow we define a measure, specifically the uniform measure mathrmdlambda(t) = w(t) mathrmd t with the weight w defined as  w mathcalW = 01 rightarrow mathbbR quad w(t) = 1This measure can be defined using the composite type Measure:m = Measure(\"uniform01\");Next, we need to compute the quadrature rule relative to the uniform measure. To do this we use the composite type Quad.q1 = Quad(n-1,m);\nnw(q)This creates a quadrature rule q named \"myq\" with n-1 nodes and weights relative to the measure m. The function nw() prints the nodes and weights. To solve the integral we call integrate()I1 = integrate(f,q1)\nprint(\"Numerical error: $(abs(1-cos(1)-I1))\")"
},

{
    "location": "numerical_integration/#Variant-2-1",
    "page": "Numerical Integration",
    "title": "Variant 2",
    "category": "section",
    "text": "There is another variant to solve the integral, which computes the quadrature rule based on the recurrence coefficients of the polynomials that are orthogonal relative to the measure m. First, we compute the orthogonal polynomials using the composite type OrthoPoly.op = OrthoPoly(\"uniform01\",n-1);\ncoeffs(op)The resulting system of orthogonal polynomials is characterized by its recursion coefficients (alpha beta), which can be extracted with the function coeffs().Now, the quadrature rule can be constructed based on op, and the integral be solved.q2 = Quad(n,op)\nnw(q)\nI2 = integrate(f,q2)\nprint(\"Numerical error: $(abs(1-cos(1)-I2))\")"
},

{
    "location": "numerical_integration/#Comparison-1",
    "page": "Numerical Integration",
    "title": "Comparison",
    "category": "section",
    "text": "We see that the different variants provide slightly different results:1-cos(1) .- [I0 I1 I2]with I0 and I2 being the same and more accurate than I1. The increased accuracy is based on the fact that for I0 and I2 the quadrature rules are based on the recursion coefficients of the underlying orthogonal polynomials. The quadrature for I1 is based on an general-purpose method that can be significantly less accurate, see also the next tutorial."
},

{
    "location": "quadrature_rules/#",
    "page": "Quadrature Rules",
    "title": "Quadrature Rules",
    "category": "page",
    "text": "using PolyChaos, LinearAlgebra\nmy_f(t) = t^2\na, b = 1.23, 3.45 # shape parameters of Jacobi weight\nint_exact = 0.353897; # reference value \nN = 4\nα, β = rm_jacobi(N,a,b)\nn_gauss, w_gauss = gauss(N,α,β)\nint_gauss = dot(w_gauss,my_f.(n_gauss))\nprint(\"first point:\\t $(n_gauss[1])\\n\")\nprint(\"end point:\\t $(n_gauss[end])\\n\")\nprint(\"error Gauss:\\t $(int_gauss-int_exact)\\n\")\nn_radau, w_radau = radau(N-1,α,β,1.)\nint_radau = dot(w_radau,my_f.(n_radau))\nprint(\"first point:\\t $(n_radau[1])\\n\")\nprint(\"end point:\\t $(n_radau[end])\\n\")\nprint(\"error Radau:\\t $(int_radau-int_exact)\")\nn_lob, w_lob = lobatto(N-2,α,β,-1.,1.)\nint_lob = dot(w_lob,my_f.(n_lob))\nprint(\"first point:\\t $(n_lob[1])\\n\")\nprint(\"end point:\\t $(n_lob[end])\\n\")\nprint(\"error Lobatto:\\t $(int_lob-int_exact)\")\nn_fej, w_fej = fejer(N)\nint_fej = dot(w_fej,my_f.(n_fej).*(1 .- n_fej).^a.*(1 .+ n_fej).^b)\nprint(\"first point:\\t $(n_fej[1])\\n\")\nprint(\"end point:\\t $(n_fej[end])\\n\")\nprint(\"error Fejer:\\t $(int_fej-int_exact)\")\nn_fej2, w_fej2 = fejer2(N)\nint_fej2 = dot(w_fej2,my_f.(n_fej2).*(1 .- n_fej2).^a.*(1 .+ n_fej2).^b)\nprint(\"first point:\\t $(n_fej2[1])\\n\")\nprint(\"end point:\\t $(n_fej2[end])\\n\")\nprint(\"error Fejer2:\\t $(int_fej2-int_exact)\")\nn_cc, w_cc = clenshaw_curtis(N)\nint_cc = dot(w_cc,my_f.(n_cc).*(1 .- n_cc).^a.*(1 .+ n_cc).^b)\nprint(\"first point:\\t\\t $(n_cc[1])\\n\")\nprint(\"end point:\\t\\t $(n_cc[end])\\n\")\nprint(\"error Clenshaw-Curtis:\\t $(int_cc-int_exact)\")"
},

{
    "location": "quadrature_rules/#QuadratureRules-1",
    "page": "Quadrature Rules",
    "title": "Quadrature Rules",
    "category": "section",
    "text": "In this tutorial we investigate how recurrence coefficients of orthogonal polynomials lead to quadrature rules.We want to solve the integralI = int_-1^1 f(t) w(t) mathrmd twith the weight functionw(t) = (1-t)^a (1+t)^bfor all t in -11 and ab-1. For the function f we choosef(t) = t^2To solve the integral we do the following:Choose number of nodes N;\nGenerate recurrence coefficients;\nGenerate quadrature rule from those recurrence coefficients.We will compare Gauss quadrature to Gauss-Radau quadrature and Gauss-Lobatto quadrature.Make sure to check out this tutorial too.Let\'s begin:using PolyChaos, LinearAlgebra\nmy_f(t) = t^2\na, b = 1.23, 3.45 # shape parameters of Jacobi weight\nint_exact = 0.353897; # reference value Now we compute N recurrence coefficients.N = 4\nα, β = rm_jacobi(N,a,b)"
},

{
    "location": "quadrature_rules/#Gauss-1",
    "page": "Quadrature Rules",
    "title": "Gauss",
    "category": "section",
    "text": "The first quadrature rule is Gauss quadrature. This method goes back to Golub and Welsch.n_gauss, w_gauss = gauss(N,α,β)\nint_gauss = dot(w_gauss,my_f.(n_gauss))\nprint(\"first point:\\t $(n_gauss[1])\\n\")\nprint(\"end point:\\t $(n_gauss[end])\\n\")\nprint(\"error Gauss:\\t $(int_gauss-int_exact)\\n\")Since Gauss quadrature has a degree of exactness of 2N-1, the value of the integral is exact."
},

{
    "location": "quadrature_rules/#Gauss-Radau-1",
    "page": "Quadrature Rules",
    "title": "Gauss-Radau",
    "category": "section",
    "text": "Gauss-Radau quadrature is a variant of Gauss quadrature that allows to specify a value of a node that has to be included. We choose to include the right end point t = 10.n_radau, w_radau = radau(N-1,α,β,1.)\nint_radau = dot(w_radau,my_f.(n_radau))\nprint(\"first point:\\t $(n_radau[1])\\n\")\nprint(\"end point:\\t $(n_radau[end])\\n\")\nprint(\"error Radau:\\t $(int_radau-int_exact)\")"
},

{
    "location": "quadrature_rules/#Gauss-Lobatto-1",
    "page": "Quadrature Rules",
    "title": "Gauss-Lobatto",
    "category": "section",
    "text": "Next, we look at Gauss-Lobatto quadrature, which allows to include two points. We choose to include the left and end point of the interval, which are t in -10 10.n_lob, w_lob = lobatto(N-2,α,β,-1.,1.)\nint_lob = dot(w_lob,my_f.(n_lob))\nprint(\"first point:\\t $(n_lob[1])\\n\")\nprint(\"end point:\\t $(n_lob[end])\\n\")\nprint(\"error Lobatto:\\t $(int_lob-int_exact)\")There are other quadratures that we subsume as all-purpose quadrature rules. These include Fejér\'s first and second rule, and Clenshaw-Curtis quadrature."
},

{
    "location": "quadrature_rules/#Fejér\'s-First-Rule-1",
    "page": "Quadrature Rules",
    "title": "Fejér\'s First Rule",
    "category": "section",
    "text": "Fejér\'s first rule does not include the end points of the interval.n_fej, w_fej = fejer(N)\nint_fej = dot(w_fej,my_f.(n_fej).*(1 .- n_fej).^a.*(1 .+ n_fej).^b)\nprint(\"first point:\\t $(n_fej[1])\\n\")\nprint(\"end point:\\t $(n_fej[end])\\n\")\nprint(\"error Fejer:\\t $(int_fej-int_exact)\")"
},

{
    "location": "quadrature_rules/#Fejér\'s-Second-Rule-1",
    "page": "Quadrature Rules",
    "title": "Fejér\'s Second Rule",
    "category": "section",
    "text": "Fejér\'s second rule does include the end points of the interval.n_fej2, w_fej2 = fejer2(N)\nint_fej2 = dot(w_fej2,my_f.(n_fej2).*(1 .- n_fej2).^a.*(1 .+ n_fej2).^b)\nprint(\"first point:\\t $(n_fej2[1])\\n\")\nprint(\"end point:\\t $(n_fej2[end])\\n\")\nprint(\"error Fejer2:\\t $(int_fej2-int_exact)\")"
},

{
    "location": "quadrature_rules/#Clenshaw-Curtis-1",
    "page": "Quadrature Rules",
    "title": "Clenshaw-Curtis",
    "category": "section",
    "text": "Clenshaw-Curtis quadrature is similar to Féjer\'s second rule, as in it includes the end points of the integration interval. For the same number of nodes it is also more accurate than Féjer\'s rules, generally speaking.n_cc, w_cc = clenshaw_curtis(N)\nint_cc = dot(w_cc,my_f.(n_cc).*(1 .- n_cc).^a.*(1 .+ n_cc).^b)\nprint(\"first point:\\t\\t $(n_cc[1])\\n\")\nprint(\"end point:\\t\\t $(n_cc[end])\\n\")\nprint(\"error Clenshaw-Curtis:\\t $(int_cc-int_exact)\")As we can see, for the same number of nodes N, the quadrature rules based on the recurrence coefficients can greatly outperform the all-purpose quadratures. So, whenever possible, use quadrature rules based on recurrence coefficients of the orthogonal polynomials relative to the underlying measure. Make sure to check out this tutorial too."
},

{
    "location": "orthogonal_polynomials_canonical/#",
    "page": "Univariate Monic Orthogonal Polynomials",
    "title": "Univariate Monic Orthogonal Polynomials",
    "category": "page",
    "text": ""
},

{
    "location": "orthogonal_polynomials_canonical/#UnivariateMonicOrthogonalPolynomials-1",
    "page": "Univariate Monic Orthogonal Polynomials",
    "title": "Univariate Monic Orthogonal Polynomials",
    "category": "section",
    "text": "Univariate monic orthogonal polynomials make up the core building block of the package. These are real polynomials  pi_k _k geq 0, which are univariate pi_k mathbbR rightarrow mathbbR and orthogonal relative to a nonnegative weight function w mathbbR rightarrow mathbbR_geq 0, and which have a leading coefficient equal to one:beginaligned\npi_k(t) = t^k + a_k-1 t^k-1 + dots + a_1 t + a_0 quad forall k = 0 1 dots \nlangle pi_k pi_l rangle = int_mathbbR pi_k(t) pi_l(t) w(t) mathrmdt =\nbegincases\n0  k neq l text and kl geq 0 \n pi_k ^2  0  k = l geq 0\nendcases\nendalignedThese univariate monic orthogonal polynomials satisfy the paramount three-term recurrence relationbeginaligned\npi_k+1(t) = (t - alpha_k) pi_k(t) - beta_k pi_k-1(t) quad k= 0 1 dots \npi_o(t) = 1 \npi_-1(t) = 0\nendalignedHence, every system of n univariate monic orthogonal polynomials  pi_k _k=0^n is isomorphic to its recurrence coefficients  alpha_k beta_k _k=0^n."
},

{
    "location": "orthogonal_polynomials_canonical/#Classical-Orthogonal-Polynomials-1",
    "page": "Univariate Monic Orthogonal Polynomials",
    "title": "Classical Orthogonal Polynomials",
    "category": "section",
    "text": "The so-called classical orthogonal polynomials are polynomials named after famous mathematicians who each discovered a special family of orthogonal polynomials, for example Hermite polynomials or Jacobi polynomials. For classical orthogonal polynomials there exist closed-form expressions of–-among others–-the recurrence coefficients. Also quadrature rules for classical orthogonal polynomials are well-studied (with dedicated packages such as FastGaussQuadrature.jl. However, more often than not these classical orthogonal polynomials are neither monic nor orthogonal, hence not normalized in any sense. For example, there is a distinction between the probabilists\' Hermite polynomials and the physicists\' Hermite polynomials. The difference is in the weight function w(t) relative to which the polynomials are orthogonal:beginaligned\ntextProbabilists w(t) = frac1sqrt2 pi  exp left( - fract^22 right) \ntextPhysicists w(t) =  exp left( - t^2 right)\nendalignedTo streamline the computations, all classical orthogonal polynomials are converted to monic orthogonal polynomials (for which, of course, the closed-form expressions persist). Currently, the following weight functions (hence classical orthogonal polynomials) are supported:Name Weight w(t) Parameters Support Classical polynomial\nhermite $ \\exp \\left( - t^2 \\right)$ - (-infty infty) Hermite\ngenhermite $ \\lvert t \\rvert^{2 \\mu}\\exp \\left( - t^2 \\right)$ mu  -frac12 (-infty infty) Generalized Hermite\nlegendre 1 - (-11) Legendre\njacobi (1-t)^alpha (1+t)^beta alpha beta  -1 (-11) Jacobi\nlaguerre exp(-t) - (0infty) Laguerre\ngenlaguerre t^alphaexp(-t) alpha-1 (0infty) Generalized Laguerre\nmeixnerpollaczek frac12 pi exp((2phi-pi)t) lvertGamma(lambda + mathrmit)rvert^2 lambda  0 0phipi (-inftyinfty) Meixner-PollaczekAdditionally, the following weight functions that are equivalent to probability density functions are supported:Name Weight w(t) Parameters Support Classical polynomial\ngaussian frac1sqrt2 pi  exp left( - fract^22 right) - (-infty infty) Probabilists\' Hermite\nuniform01 1 - (01) Legendre\nbeta01 frac1B(alphabeta)  t^alpha-1 (1-t)^beta-1 alpha beta  0 (01) Jacobi\ngamma fracbeta^alphaGamma(alpha) t^alpha-1 exp(-beta t) alpha beta  0 (0infty) Laguerre\nlogistic fracexp(-t)(1+exp(-t))^2 - (-inftyinfty) -To generate the orthogonal polynomials up to maximum degree deg, simply callusing PolyChaos\ndeg = 4\nop = OrthoPoly(\"gaussian\",deg)This generates opas an OrthoPoly type with the underlying Gaussian measure op.meas. The recurrence coefficients are accessible via coeffs().coeffs(op)By default, the constructor for OrthoPoly generates deg+1 recurrence coefficients. Sometimes, some other number Nrec may be required. This is why Nrec is a keyword for the constructor OrthoPoly.N = 100\nop_ = OrthoPoly(\"logistic\",deg;Nrec=N)Let\'s check whether we truly have more coefficients:size(coeffs(op_),1)==N"
},

{
    "location": "orthogonal_polynomials_canonical/#Arbitrary-Weights-1",
    "page": "Univariate Monic Orthogonal Polynomials",
    "title": "Arbitrary Weights",
    "category": "section",
    "text": "If you are given a weight function w that does not belong to the Table above, it is still possible to generate the respective univariate monic orthogonal polynomials. First, we define the measure by specifying a name, the weight, the support, symmetry, and parameterssupp = (-1,1)\nfunction w(t)\n    supp[1]<=t<=supp[2] ? (1. + t) : error(\"$t not in support\")\nend\nmy_meas = Measure(\"my_meas\",w,supp,false,Dict())Notice: it is advisable to define the weight such that an error is thrown for arguments outside of the support.Now, we want to construct the univariate monic orthogonal polynomials up to degree deg relative to m. The constructor ismy_op = OrthoPoly(\"my_op\",deg,my_meas;Nquad=200)By default, the recurrence coefficients are computed using the Stieltjes procuedure with Clenshaw-Curtis quadrature (with Nquad nodes and weights). Hence, the choice of Nquad influences accuracy."
},

{
    "location": "orthogonal_polynomials_canonical/#MultivariateMonicOrthogonalPolynomials-1",
    "page": "Univariate Monic Orthogonal Polynomials",
    "title": "Multivariate Monic Orthogonal Polynomials",
    "category": "section",
    "text": "Suppose we have p systems of univariate monic orthogonal polynomials, pi_k^(1) _kgeq 0   pi_k^(2) _kgeq 0 dots  pi_k^(p) _kgeq 0each system being orthogonal relative to the weights w^(1) w^(2) dots w^(p) with supports mathcalW^(1) mathcalW^(2) dots mathcalW^(p). Also, let d^(i) be the maximum degree of the i-th system of univariate orthogonal polynomials. We would like to construct a p-variate monic basis  psi_k _k geq 0 with psi mathbbR^p rightarrow mathbbR of degree at most 0 leq d leq min_i=1dotsk d^(i). Further, this basis shall be orthogonal relative to the product measure w mathcalW = mathcalW^(1) otimes mathcalW^(2) mathcalW^(1) cdots otimes mathcalW^(p) rightarrow mathbbR_geq 0 given byw(t) = prod_i=1^p w^(i)(t_i)hence satisfieslangle psi_k psi_l rangle = int_mathcalW psi_k(t) psi_l(t) w(t) mathrmd t =\nbegincases\n0  k neq l text and kl geq 0 \n psi_k ^2  0  k = l geq 0\nendcasesFor this, there exists the composite struct MultiOrthoPoly. Let\'s consider an example where we mix classical orthogonal polynomials with an arbitrary weight.deg = [3,5,6,4]\nd = minimum(deg)\n\nop1 = OrthoPoly(\"gaussian\",deg[1])\nop2 = OrthoPoly(\"uniform01\",deg[2])\nop3 = OrthoPoly(\"beta01\",deg[3],Dict(:shape_a=>2,:shape_b=>1.2))\nops = [op1,op2,op3,my_op]\nmop = MultiOrthoPoly(ops,d)The total number of  basis polynomials is stored in the field dim. The univariate basis polynomials making up the multivariate basis are stored in the field uni.mop.uniThe field ind contains the multi-index, i.e. row i stores what combination of univariate polynomials makes up the i-th multivariate polynomial. For example,i = 11\nmop.ind[i+1,:]translates mathematically topsi_11(t) = pi_0^(1)(t_1) pi_1^(2)(t_2) pi_0^(3)(t_3) pi_1^(4)(t_4)Notice that there is an offset by one, because the basis counting starts at 0, but Julia is 1-indexed. The underlying measure of mop is now of type MultiMeasure, and stored in the field measmop.measThe weight w can be evaluated as expectedmop.meas.w(0.5*ones(length(ops)))"
},

{
    "location": "multiple_discretization/#",
    "page": "Multiple Discretization",
    "title": "Multiple Discretization",
    "category": "page",
    "text": "using PolyChaos, LinearAlgebra\nγ = 0.5;\nint_exact = 1+pi/2; # exact value of the integral\nfunction my_w(t::Float64,γ::Float64)\n    γ + (1-γ)*1/sqrt(1-t^2)\nend\nN = 1000;\nn,w = fejer(N);\nint_fejer = dot(w,my_w.(n,γ))\nprint(\"Fejer error:\\t$(abs(int_exact-int_fejer))\\twith $N nodes\")\nfunction quad_gaussleg(N::Int,γ::Float64)\n    a,b=rm_legendre(N)\n    n,w=golubwelsch(a,b)\nend\nn,w = quad_gaussleg(N,γ)\nint_gaussleg = dot(w,γ .+ (1-γ)/sqrt.(1 .- n.^2))\nprint(\"Gauss-Legendre error:\\t$(abs(int_exact-int_gaussleg))\\twith $N nodes\")\nfunction quad_gausscheb(N::Int64,γ::Float64)\n    a,b = rm_chebyshev1(N)\n    n,w = golubwelsch(a,b)\nend\nn,w = quad_gausscheb(N,γ)\nint_gausscheb = dot(w,γ .+ (1-γ)*sqrt.(1 .- n.^2))\nprint(\"Gauss-Chebyshev error:\\t$(abs(int_exact-int_gausscheb))\\twith $N nodes\")\nfunction quad_gaussleg_mod(N::Int,γ::Float64)\n    n,w = quad_gaussleg(N,γ)\n    return n,γ*w\nend\nfunction quad_gausscheb_mod(N::Int,γ::Float64)\n            n,w = quad_gausscheb(N,γ)\n    return n,(1-γ)*w\nend\nN = 8\na,b = mcdiscretization(N,[n->quad_gaussleg_mod(n,γ); n->quad_gausscheb_mod(n,γ)])\nn,w = golubwelsch(a,b)\nint_mc = sum(w)\nprint(\"Discretization error:\\t$(abs(int_exact-int_mc))\\twith $N nodes\")\nΓ = 0:0.1:1;\nab = [ mcdiscretization(N,[n->quad_gaussleg_mod(n,gam); n->quad_gausscheb_mod(n,gam)]) for gam in Γ ];\nbb = hcat([ ab[i][2] for i=1:length(Γ)]...);\nb_leg = rm_legendre(N)[2];\nb_cheb = rm_chebyshev1(N)[2]\nbb[:,1]-b_cheb\nbb[:,end]-b_leg\nusing Plots\ngr()\nplot(Γ,bb\',yaxis=:log10, w=3, legend=false)\nzs, os = zeros(N), ones(N)\nscatter!(zs,b_cheb,marker=:x)\nscatter!(os,b_leg,marker=:circle)\nxlabel!(\"gamma\")\nylabel!(\"beta\")\ntitle!(\"N=$N Recurrence Coefficients Interpolating from Chebyshev to Legendre\")"
},

{
    "location": "multiple_discretization/#Multiple-Discretization-1",
    "page": "Multiple Discretization",
    "title": "Multiple Discretization",
    "category": "section",
    "text": "This tutorial shows how to compute recurrence coefficients for non-trivial weight functions, and how they are being used for quadrature. The method we use is called multiple discretization, and follows W. Gautschi\'s book \"Orthogonal Polynomials: Computation and Approximation\", specifically Section 2.2.4, and Example 2.38.Suppose we have the weight functionforall t in -11 gamma in 01 quad w(tgamma) = gamma + (1-gamma) frac1sqrt1-t^2and we would like to solveint_-1^1 f(t) w(tc) mathrmdt = sum_nu=1^N f(tau_nu) w_nuby some quadrature rule. We will see that ad-hoc quadrature rules will fail to solve the integral even for the simplest choice f equiv 1. However, finding the recurrence coefficients of the underlying orthogonal polynomials, and then finding the quadrature rule will be the way to go.Let us first try to solve the integral for f equiv 1 by Fejer\'s rule.using PolyChaos, LinearAlgebra\nγ = 0.5;\nint_exact = 1+pi/2; # exact value of the integral\nfunction my_w(t::Float64,γ::Float64)\n    γ + (1-γ)*1/sqrt(1-t^2)\nend\n\nN = 1000;\nn,w = fejer(N);\nint_fejer = dot(w,my_w.(n,γ))\nprint(\"Fejer error:\\t$(abs(int_exact-int_fejer))\\twith $N nodes\")Clearly, that is not satisfying. Well, the term gamma of the weight w makes us think of Gauss-Legendre integration, so let\'s try it instead.function quad_gaussleg(N::Int,γ::Float64)\n    a,b=rm_legendre(N)\n    n,w=golubwelsch(a,b)\nend\nn,w = quad_gaussleg(N,γ)\nint_gaussleg = dot(w,γ .+ (1-γ)/sqrt.(1 .- n.^2))\nprint(\"Gauss-Legendre error:\\t$(abs(int_exact-int_gaussleg))\\twith $N nodes\")Even worse! Well, we can factor out frac1sqrt1-t^2, making the integral amenable to a Gauss-Chebyshev rule. So, let\'s give it anothery try.function quad_gausscheb(N::Int64,γ::Float64)\n    a,b = rm_chebyshev1(N)\n    n,w = golubwelsch(a,b)\nend\nn,w = quad_gausscheb(N,γ)\nint_gausscheb = dot(w,γ .+ (1-γ)*sqrt.(1 .- n.^2))\n# int=sum(xw(:,2).*(1+sqrt(1-xw(:,1).^2)))\nprint(\"Gauss-Chebyshev error:\\t$(abs(int_exact-int_gausscheb))\\twith $N nodes\")Okay, that\'s better, but it took us a lot of nodes to get this result. Is there a different way? Indeed, there is. As we have noticed, the weight w has a lot in common with Gauss-Legendre and Gauss-Chebyshev. We can decompose the integral as followsint_-1^1 f(t) w(t) mathrmdt = sum_i=1^m int_-1^1 f(t) w_i(t) mathrmd twithbeginalign*\nw_1(t) = gamma \nw_2(t) = (1-gamma) frac1sqrt1-t^2\nendalign*To the weight w_1 we can apply Gauss-Legendre quadrature, to the weight w_2 we can apply Gauss-Chebyshev quadrature (with tiny modifications). This discretization of the measure can be used in our favor. The function mcdiscretization() takes the m discretization rules as an inputfunction quad_gaussleg_mod(N::Int,γ::Float64)\n    n,w = quad_gaussleg(N,γ)\n    return n,γ*w\nend\nfunction quad_gausscheb_mod(N::Int,γ::Float64)\n            n,w = quad_gausscheb(N,γ)\n    return n,(1-γ)*w\nend\n\nN = 8\na,b = mcdiscretization(N,[n->quad_gaussleg_mod(n,γ); n->quad_gausscheb_mod(n,γ)])\nn,w = golubwelsch(a,b)\nint_mc = sum(w)\nprint(\"Discretization error:\\t$(abs(int_exact-int_mc))\\twith $N nodes\")Et voilà, no error with fewer nodes. (For this example, we\'d need in fact just a single node.)The function mcdiscretization() is able to construct the recurrence coefficients of the orthogonal polynomials relative to the weight w. Let\'s inspect the values of the recurrence coefficients a little more. For gamma = 0, we are in the world of Chebyshev polynomials, for gamma = 1, we enter the realm of Legendre polynomials. And in between? That\'s exactly where the weight w comes in: it can be thought of as an interpolatory weight, interpolating Legendre polynomials and Chebyshev polynomials. Let\'s verify this by plotting the recurrence coefficients for several values of gamma.Γ = 0:0.1:1;\nab = [ mcdiscretization(N,[n->quad_gaussleg_mod(n,gam); n->quad_gausscheb_mod(n,gam)]) for gam in Γ ];\nbb = hcat([ ab[i][2] for i=1:length(Γ)]...);\nb_leg = rm_legendre(N)[2];\nb_cheb = rm_chebyshev1(N)[2]\nbb[:,1]-b_chebbb[:,end]-b_legLet\'s plot these values to get a better feeling.using Plots\ngr()\nplot(Γ,bb\',yaxis=:log10, w=3, legend=false)\nzs, os = zeros(N), ones(N)\nscatter!(zs,b_cheb,marker=:x)\nscatter!(os,b_leg,marker=:circle)\n\nxlabel!(\"gamma\")\nylabel!(\"beta\")\ntitle!(\"N=$N Recurrence Coefficients Interpolating from Chebyshev to Legendre\")The crosses denote the values of the β recursion coefficients for Chebyshev polynomials; the circles the β recursion coefficients for Legendre polynomials. The interpolating line in between stands for the β recursion coefficients of w(tgamma)."
},

{
    "location": "scalar_products/#",
    "page": "Computation of Scalar Products",
    "title": "Computation of Scalar Products",
    "category": "page",
    "text": ""
},

{
    "location": "scalar_products/#ComputationOfScalarProducts-1",
    "page": "Computation of Scalar Products",
    "title": "Computation of Scalar Products",
    "category": "section",
    "text": "By now, we are able to construct orthogonal polynomials, and to construct quadrature rules for a given nonnegative weight function, respectively. Now we combine both ideas to solve integrals involving the orthogonal polynomialslangle phi_i_1 phi_i_2 cdots phi_i_m-1 phi_i_m rangle\n= int phi_i_1(t) phi_i_2(t) cdots phi_i_m-1(t) phi_i_m(t) w(t) mathrmd tboth for the univariate and multivariate case. The integrand is a polynomial (possibly multivariate) that can be solved exactly with the appropriate Gauss quadrature rules.note: Note\nTo simplify notation we drop the integration interval. It is clear from the context."
},

{
    "location": "scalar_products/#Univariate-Polynomials-1",
    "page": "Computation of Scalar Products",
    "title": "Univariate Polynomials",
    "category": "section",
    "text": ""
},

{
    "location": "scalar_products/#Classical-Polynomials-1",
    "page": "Computation of Scalar Products",
    "title": "Classical Polynomials",
    "category": "section",
    "text": "Let\'s begin with a univariate basis for some classical orthogonal polynomialusing PolyChaos\ndeg, n = 4, 20\ns_α, s_β = 2.1, 3.2\nop = OrthoPoly(\"beta01\",deg,Dict(:shape_a=>s_α,:shape_b=>s_β);Nrec=n)To add the corresponding quadrature rule there is the composite struct OrthoPolyQ whose simplest constructor readsopq = OrthoPolyQ(op,n)By default, an n-point Gauss quadrature rule is create relative to the underlying measure op.meas, where n is the number of recurrence coefficients stored in op.α and op.β. The type OrthoPolyQ has just two fields: an OrthoPoly, and a Quad.To compute the squared norms phi_k ^2 = langle phi_k phi_k  rangle\n= int phi_k(t) phi_k(t) w(t) mathrmd tof the basis we call computeSP2()normsq = computeSP2(opq)For the general caselangle phi_i_1 phi_i_2 cdots phi_i_m-1 phi_i_m rangle\n= int phi_i_1(t) phi_i_2(t) cdots phi_i_m-1(t) phi_i_m(t) w(t) mathrmd tthere exists a type Tensor that requires only two arguments: the dimension m geq 1, and an OrthoPolyQm = 3\nt = Tensor(3,opq)To get the desired entries, Tensorcomes with a get() function that is called for some index a in mathbbN_0^m that has the entries a = i_1 i_2 dots i_m. For examplet.get([1,2,3])Or using comprehensionT = [ t.get([i1,i2,i3]) for i1=0:dim(opq)-1,i2=0:dim(opq)-1,i3=0:dim(opq)-1]Notice that we can cross-check the results.using LinearAlgebra\n@show normsq == LinearAlgebra.diag(T[:,:,1])\n@show normsq == LinearAlgebra.diag(T[:,1,:])\n@show normsq == LinearAlgebra.diag(T[1,:,:])Also, normsq can be computed analogously in Tensor formatt2 = Tensor(2,opq)\n@show normsq == [ t2.get([i,i]) for i=0:dim(opq)-1]"
},

{
    "location": "scalar_products/#Arbitrary-Weights-1",
    "page": "Computation of Scalar Products",
    "title": "Arbitrary Weights",
    "category": "section",
    "text": "Of course, the type OrthoPolyQ can be constructed for arbitrary weights w(t). In this case we have to compute the orthogonal basis and the respective quadrature rule. Let\'s re-work the above example by hand.using SpecialFunctions\nsupp = (0,1)\nfunction w(t)\n    supp[1]<=t<=supp[2] ? (t^(s_α-1)*(1-t)^(s_β-1)/SpecialFunctions.beta(s_α,s_β)) : error(\"$t not in support\")\nend\nmy_meas = Measure(\"my_meas\",w,supp,false,Dict())\nmy_op = OrthoPoly(\"my_op\",deg,my_meas;Nrec=n)\nmy_quad = Quad(n,my_op)\nmy_opq = OrthoPolyQ(my_op,my_quad)Now we can compute the squared norms  phi_k ^2my_normsq = computeSP2(my_opq)And the tensormy_t = Tensor(m,my_opq)\nmy_T = [ my_t.get([i1,i2,i3]) for i1=0:dim(opq)-1,i2=0:dim(opq)-1,i3=0:dim(opq)-1]Let\'s compare the results:@show abs.(normsq-my_normsq)\n@show norm(T-my_T)note: Note\nThe possibility to create quadrature rules for arbitrary weights should be reserved to cases different from classical ones."
},

{
    "location": "scalar_products/#Multivariate-Polynomials-1",
    "page": "Computation of Scalar Products",
    "title": "Multivariate Polynomials",
    "category": "section",
    "text": "For multivariate polynomials the syntax for Tensor is very much alike, except that we are dealing with the type MultiOrthoPoly now.mop = MultiOrthoPoly([opq,my_opq],deg)mt2 = Tensor(2,mop)\nmt3 = Tensor(3,mop)\nmT2 = [ mt2.get([i,i]) for i=0:dim(mop)-1 ]Notice that mT2 carries the elements of the 2-dimensional tensors for the univariate bases opq and my_opq. The encoding is given by the multi-index mop.indmop.indTo cross-check the results we can distribute the multi-index back to its univariate indices with the help of findUnivariateIndices.ind_opq = findUnivariateIndices(1,mop.ind)\nind_my_opq = findUnivariateIndices(2,mop.ind)@show mT2[ind_opq] - normsq\n@show mT2[ind_my_opq] - my_normsq;"
},

{
    "location": "pce_tutorial/#",
    "page": "Basic Usage",
    "title": "Basic Usage",
    "category": "page",
    "text": ""
},

{
    "location": "pce_tutorial/#CommonRandomVariables-1",
    "page": "Basic Usage",
    "title": "Common Random Variables",
    "category": "section",
    "text": "Polynomial chaos expansion (PCE) is a Hilbert space technique for random variables with finite variance. Mathematically equivalent to Fourier series expansions for periodic signals, PCE allows to characterize a random variable in terms of its PCE coefficients (aka Fourier coefficients). That is, the PCE of a random variable mathsfx is given bymathsfx = sum_i=0^L x_i phi_iwhere x_i are the so-called PCE coefficients, and phi_i are the orthogonal polynomials that are orthogonal relative to the probability density function of mathsfx.This tutorial walks you through the PCE of common random variables, namely Gaussian (gaussian), Beta (beta01), Uniform(uniform01), Logistic (logistic), and shows how they are implemented in PolyChaos."
},

{
    "location": "pce_tutorial/#Construction-of-Basis-1",
    "page": "Basic Usage",
    "title": "Construction of Basis",
    "category": "section",
    "text": "We begin by specifying the names and, if any, parameters for the uncertainties.using Revise\nusing PolyChaos\nα, β = 1.3, 2.2\npolynames = [\"gaussian\", \"beta01\", \"uniform01\", \"logistic\"]\npars = [Dict(), Dict(:shape_a=>α,:shape_b=>β), Dict(), Dict()]The orthogonal polynomials are constructed using OrthoPoly (here of degree at most d, and stored in the dictionary myops).d = 6\nmyops = Dict()\nfor (i,name) in enumerate(polynames)\n    myops[name]=OrthoPoly(name,d,pars[i])\nendFor example, let\'s evaluate the Gaussian basis polynomials at some pointspoints, degrees = randn(10), 0:2:d\nop_gauss=myops[\"gaussian\"]\n[ evaluate(degree,points,op_gauss) for degree in degrees ]If a quadrature rule is required, this can be added by calling OrthoPolyQmyopqs = Dict()\nfor (i,name) in enumerate(polynames)\n    myopqs[name]=OrthoPolyQ(name,d,pars[i])\nend"
},

{
    "location": "pce_tutorial/#Finding-PCE-Coefficients-1",
    "page": "Basic Usage",
    "title": "Finding PCE Coefficients",
    "category": "section",
    "text": "Having constructed the orthogonal bases, the question remains how to find the PCE coefficients for the common random variables. Every random variable can be characterized exactly by two PCE coefficients. For a Gaussian random variable, this is familiar: the mean and the variance suffice to describe a Gaussian random variable entirely. The same is true for any random variable of finite variance given the right basis. The function convert2affinePCE provides the first two PCE coefficients (hence the name affine) for the common random variables."
},

{
    "location": "pce_tutorial/#Gaussian-1",
    "page": "Basic Usage",
    "title": "Gaussian",
    "category": "section",
    "text": "Given the Gaussian random variable mathsfx sim mathcalN(mu sigma^2) with sigma  0, the affine PCE coefficients are# Gaussian\nμ, σ = 2., 0.2\npce_gaussian = convert2affinePCE(\"gaussian\",μ,σ)\n# Uniform"
},

{
    "location": "pce_tutorial/#Uniform-1",
    "page": "Basic Usage",
    "title": "Uniform",
    "category": "section",
    "text": "Given the uniform random variable mathsfx sim mathcalU(a b) with finite support ab, the affine PCE coefficients area, b = -0.3, 1.2\nconvert2affinePCE(\"uniform01\",a,b)Instead, if the expected value and standard deviation are known, the affine PCE coefficients of the uniform random variable arepce_uniform = convert2affinePCE(\"uniform01\",μ,σ;kind=:μσ)\n# notice that the zero-order coefficient IS equivalent to the expected value μ"
},

{
    "location": "pce_tutorial/#Beta-1",
    "page": "Basic Usage",
    "title": "Beta",
    "category": "section",
    "text": "Given the Beta random variable mathsfx sim mathcalB(a b alpha beta) with finite support ab and shape parameters alpha beta  0, the affine PCE coefficients areconvert2affinePCE(\"beta01\",a,b,Dict(:shape_a=>α,:shape_b=>β))Instead, if the expected value and standard deviation are known, the affine PCE coefficients of the uniform random variable arepce_beta = convert2affinePCE(\"beta01\",μ,σ,Dict(:shape_a=>α,:shape_b=>β); kind=:μσ)"
},

{
    "location": "pce_tutorial/#Logistic-1",
    "page": "Basic Usage",
    "title": "Logistic",
    "category": "section",
    "text": "Given the logstic random variable mathsfx sim mathcalL(a_1a_2) where a_20 with the probability density functionrho(t) = frac14 a_2  operatornamesech^2 left(fract-a_12a_2right)the affine PCE coefficients of the uniform random variable area1, a2 = μ, sqrt(3)*σ/pi\npce_logistic = convert2affinePCE(\"logistic\",a1,a2)"
},

{
    "location": "pce_tutorial/#Moments-1",
    "page": "Basic Usage",
    "title": "Moments",
    "category": "section",
    "text": "It is a key feature of PCE to compute moments from the PCE coefficients alone; no sampling is required."
},

{
    "location": "pce_tutorial/#Gaussian-2",
    "page": "Basic Usage",
    "title": "Gaussian",
    "category": "section",
    "text": "mean(pce_gaussian,myops[\"gaussian\"]), std(pce_gaussian,myops[\"gaussian\"])"
},

{
    "location": "pce_tutorial/#Uniform-2",
    "page": "Basic Usage",
    "title": "Uniform",
    "category": "section",
    "text": "mean(pce_uniform,myops[\"uniform01\"]), std(pce_uniform,myops[\"uniform01\"])"
},

{
    "location": "pce_tutorial/#Beta-2",
    "page": "Basic Usage",
    "title": "Beta",
    "category": "section",
    "text": "mean(pce_beta,myops[\"beta01\"]), std(pce_beta,myops[\"beta01\"])"
},

{
    "location": "pce_tutorial/#Logistic-2",
    "page": "Basic Usage",
    "title": "Logistic",
    "category": "section",
    "text": "mean(pce_logistic,myops[\"logistic\"]), std(pce_logistic,myops[\"logistic\"])"
},

{
    "location": "pce_tutorial/#Sampling-1",
    "page": "Basic Usage",
    "title": "Sampling",
    "category": "section",
    "text": "Having found the PCE coefficients, it may be useful to sample the random variables. That means, find N realizations of the random variable that obey the random variable\'s probability density function. This is done in two steps:Draw N samples from the measure (sampleMeasure()), and then\nEvaluate the basis polynomials and multiply times the PCE coefficients, i.e. sum_i=0^L x_i phi_i(xi_j) where xi_j is the j-th sample from the measure (evaluatePCE()).Both steps are combined in the function samplepCE()."
},

{
    "location": "pce_tutorial/#Gaussian-3",
    "page": "Basic Usage",
    "title": "Gaussian",
    "category": "section",
    "text": "using Statistics\nN = 1000\nξ_gaussian = sampleMeasure(N,myops[\"gaussian\"])\nsamples_gaussian = evaluatePCE(pce_gaussian,ξ_gaussian,myops[\"gaussian\"])\n# samplePCE(N,pce_gaussian,myops[\"gaussian\"])"
},

{
    "location": "pce_tutorial/#Uniform-3",
    "page": "Basic Usage",
    "title": "Uniform",
    "category": "section",
    "text": "ξ_uniform = sampleMeasure(N,myops[\"uniform01\"])\nsamples_uniform = evaluatePCE(pce_uniform,ξ_uniform,myops[\"uniform01\"])\n# samples_uniform = samplePCE(N,pce_uniform,myops[\"uniform01\"])"
},

{
    "location": "pce_tutorial/#Beta-3",
    "page": "Basic Usage",
    "title": "Beta",
    "category": "section",
    "text": "ξ_beta = sampleMeasure(N,myops[\"beta01\"])\nsamples_beta = evaluatePCE(pce_beta,ξ_beta,myops[\"beta01\"])\n# samples_beta = samplePCE(N,pce_beta,myops[\"beta01\"])"
},

{
    "location": "pce_tutorial/#Logistic-3",
    "page": "Basic Usage",
    "title": "Logistic",
    "category": "section",
    "text": "ξ_logistic = sampleMeasure(N,myops[\"logistic\"])\nsamples_logistic = evaluatePCE(pce_logistic,ξ_logistic,myops[\"logistic\"])\n# samples_logistic = samplePCE(N,pce_logistic,myops[\"logistic\"])"
},

{
    "location": "chi_squared_k1/#",
    "page": "Chi Squared, One DOF",
    "title": "Chi Squared, One DOF",
    "category": "page",
    "text": "k = 1\nusing PolyChaos\ndeg, Nrec = 2, 20\nop = OrthoPoly(\"gaussian\",deg;Nrec=Nrec);\nopq = OrthoPolyQ(op) #OR: opq = OrthoPolyQ(\"gaussian\",deg;Nrec=Nrec)\nL = dim(opq)\nmu, sig = 0., 1.\nx = [ convert2affinePCE(\"gaussian\",mu,sig); zeros(Float64,L-2) ]\nt2 = Tensor(2,opq);\nt3 = Tensor(3,opq)\ny = [ sum( x[i]*x[j]*t3.get([i-1,j-1,m-1])/t2.get([m-1,m-1])  for i=1:L, j=1:L ) for m=1:L ]\nmoms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]\nfunction myskew(y)\n   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )\n   μ = y[1]\n   σ = std(y,opq)\n   (e3-3*μ*σ^2-μ^3)/(σ^3)\nend\nprint(\"Expected value:\\t\\t$(moms_analytic(k)[1]) = $(mean(y,opq))\\n\")\nprint(\"\\t\\t\\terror = $(abs(mean(y,opq)-moms_analytic(k)[1]))\\n\")\nprint(\"Standard deviation:\\t$(moms_analytic(k)[2]) = $(std(y,opq))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[2]-std(y,opq))\\n\")\nprint(\"Skewness:\\t\\t$(moms_analytic(k)[3]) = $(myskew(y))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[3]-myskew(y))\\n\")\nusing Plots\ngr()\nNsmpl = 10000\nysmpl = samplePCE(Nsmpl,y,opq)\nhistogram(ysmpl;normalize=true,xlabel=\"t\",ylabel=\"rho(t)\")\nimport SpecialFunctions: gamma\nρ(t) = 1/(sqrt(2)*gamma(0.5))*1/sqrt(t)*exp(-0.5*t)\nt = range(0.1; stop=maximum(ysmpl), length=100)\nplot!(t,ρ.(t),w=4)"
},

{
    "location": "chi_squared_k1/#Chi-squared-Distribution-(k1)-1",
    "page": "Chi Squared, One DOF",
    "title": "Chi-squared Distribution (k=1)",
    "category": "section",
    "text": ""
},

{
    "location": "chi_squared_k1/#Theory-1",
    "page": "Chi Squared, One DOF",
    "title": "Theory",
    "category": "section",
    "text": "Given a standard random variable X sim mathcalN(01) we would like to find the random variable Y = X^2. The analytic solution is known: Y follows a chi-squared distribution with k=1 degree of freedom.Using polynomial chaos expansion (PCE), the problem can be solved using Galerkin projection. Let phi_k _k=0^n be the monic orthogonal basis relative to the probability density of X, namelyf_X(x) = frac1sqrt2 pi exp left( - fracx^22 right)Then, the PCE of X is given byX = sum_k=0^n x_k phi_kwithx_0 = 0 quad x_1 = 1 quad x_i = 0 quad forall i =2dotsnTo find the PCE coefficients y_k for Y = sum_k=^n y_k phi_k, we apply Galerkin projection, which leads toy_m langle phi_m phi_m rangle = sum_i=0^n sum_j=0^n x_i x_j langle phi_i phi_j phi_m rangle quad forall m = 0 dots nHence, knowing the scalars langle phi_m phi_m rangle, and langle phi_i phi_j phi_m rangle, the PCE coefficients y_k can be obtained immediately. From the PCE coefficients we can get the moments and compare them to the closed-form expressions.Notice: A maximum degree of 2 suffices to get the exact solution with PCE. In other words, increasing the maximum degree to values greater than 2 introduces nothing but computational overhead (and numerical errors, possibly)."
},

{
    "location": "chi_squared_k1/#Practice-1",
    "page": "Chi Squared, One DOF",
    "title": "Practice",
    "category": "section",
    "text": "First, we create a orthogonal basis relative to f_X(x) of degree at most d=2 (deg below).Notice that we consider a total of Nrec recursion coefficients, and that we also add a quadrature rule by calling OrthoPolyQ().k = 1\nusing PolyChaos\ndeg, Nrec = 2, 20\nop = OrthoPoly(\"gaussian\",deg;Nrec=Nrec);\nopq = OrthoPolyQ(op) #OR: opq = OrthoPolyQ(\"gaussian\",deg;Nrec=Nrec)Next, we define the PCE for X.L = dim(opq)\nmu, sig = 0., 1.\nx = [ convert2affinePCE(\"gaussian\",mu,sig); zeros(Float64,L-2) ]With the orthogonal basis and the quadrature at hand, we can compute the tensors t2 and t3 that store the entries langle phi_m phi_m rangle, and langle phi_i phi_j phi_m rangle, respectively.t2 = Tensor(2,opq);\nt3 = Tensor(3,opq)With the tensors at hand, we can compute the Galerkin projection.y = [ sum( x[i]*x[j]*t3.get([i-1,j-1,m-1])/t2.get([m-1,m-1])  for i=1:L, j=1:L ) for m=1:L ]Let\'s compare the moments via PCE to the closed-form expressions.moms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]\nfunction myskew(y)\n   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )\n   μ = y[1]\n   σ = std(y,opq)\n   (e3-3*μ*σ^2-μ^3)/(σ^3)\nend\n\nprint(\"Expected value:\\t\\t$(moms_analytic(k)[1]) = $(mean(y,opq))\\n\")\nprint(\"\\t\\t\\terror = $(abs(mean(y,opq)-moms_analytic(k)[1]))\\n\")\nprint(\"Standard deviation:\\t$(moms_analytic(k)[2]) = $(std(y,opq))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[2]-std(y,opq))\\n\")\nprint(\"Skewness:\\t\\t$(moms_analytic(k)[3]) = $(myskew(y))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[3]-myskew(y))\\n\")\nLet\'s plot the probability density function to compare results. We first draw samples from the measure with the help of sampleMeasure(), and then evaluate the basis at these samples and multiply times the PCE coefficients. The latter stop is done using evaluatePCE(). Finally, we compare the result agains the analytical PDF rho(t) = fracmathrme^-05tsqrt2 t  Gamma(05) of the chi-squared distribution with one degree of freedom.using Plots\ngr()\nNsmpl = 10000\n#ξ = sampleMeasure(Nsmpl,opq)\n#ysmpl = evaluatePCE(y,ξ,opq)\nysmpl = samplePCE(Nsmpl,y,opq)\nhistogram(ysmpl;normalize=true,xlabel=\"t\",ylabel=\"rho(t)\")\n\n\nimport SpecialFunctions: gamma\nρ(t) = 1/(sqrt(2)*gamma(0.5))*1/sqrt(t)*exp(-0.5*t)\nt = range(0.1; stop=maximum(ysmpl), length=100)\nplot!(t,ρ.(t),w=4)"
},

{
    "location": "chi_squared_k_greater1/#",
    "page": "Chi Squared, Several DOFs",
    "title": "Chi Squared, Several DOFs",
    "category": "page",
    "text": "k = 12\nusing PolyChaos\ndegree, Nrec = 2, 20\nop = OrthoPoly(\"gaussian\",degree;Nrec=Nrec);\nopq = OrthoPolyQ(op) #OR: opq = OrthoPolyQ(\"gaussian\",deg;Nrec=Nrec)\nmop = MultiOrthoPoly([opq for i=1:k],degree)\nL = dim(mop)\nmu, sig = 0., 1.\nx = [ assign2multi(convert2affinePCE(\"gaussian\",mu,sig),i,mop.ind) for i=1:k ]\nt2 = Tensor(2,mop)\nt3 = Tensor(3,mop)\ny = [ sum( x[i][j1]*x[i][j2]*t3.get([j1-1,j2-1,m-1])/t2.get([m-1,m-1])  for i=1:k, j1=1:L, j2=1:L ) for m=1:L ]\nmoms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]\nfunction myskew(y)\n   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )\n   μ = y[1]\n   σ = std(y,mop)\n   (e3-3*μ*σ^2-μ^3)/(σ^3)\nend\nprint(\"Expected value:\\t\\t$(moms_analytic(k)[1]) = $(mean(y,mop))\\n\")\nprint(\"\\t\\t\\terror = $(abs(mean(y,mop)-moms_analytic(k)[1]))\\n\")\nprint(\"Standard deviation:\\t$(moms_analytic(k)[2]) = $(std(y,mop))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[2]-std(y,mop))\\n\")\nprint(\"Skewness:\\t\\t$(moms_analytic(k)[3]) = $(myskew(y))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[3]-myskew(y))\\n\")\nusing Plots\ngr()\nNsmpl = 10000\nysmpl = samplePCE(Nsmpl,y,mop)\nhistogram(ysmpl;normalize=true,xlabel=\"t\",ylabel=\"rho(t)\")\nimport SpecialFunctions: gamma\nρ(t) = 1/(2^(0.5*k)*gamma(0.5*k))*t^(0.5*k-1)*exp(-0.5*t)\nt = range(0.1; stop=maximum(ysmpl), length=100)\nplot!(t,ρ.(t),w=4)"
},

{
    "location": "chi_squared_k_greater1/#Chi-squared-Distribution-(k1)-1",
    "page": "Chi Squared, Several DOFs",
    "title": "Chi-squared Distribution (k1)",
    "category": "section",
    "text": ""
},

{
    "location": "chi_squared_k_greater1/#Theory-1",
    "page": "Chi Squared, Several DOFs",
    "title": "Theory",
    "category": "section",
    "text": "Given k standard random variables X_i sim mathcalN(01) for i=1dotsk we would like to find the random variable Y = sum_i=1^k X_i^2. The analytic solution is known: Y follows a chi-squared distribution with k degrees of freedom.Using polynomial chaos expansion (PCE), the problem can be solved using Galerkin projection. Let phi_k _k=0^n be the monic orthogonal basis relative to the probability density of X = X_1 dots X_k, namelyf_X(x) =  prod_i=1^k frac1sqrt2 pi  exp left( - fracx_i^22 right)Then, the PCE of X_i is given byX_i = sum_k=0^n x_ik phi_kwithx_i0 = 0 quad x_ii+1 = 1 quad x_il = 0 quad forall l in 1dotsn setminus i+1To find the PCE coefficients y_k for Y = sum_k=^n y_k phi_k, we apply Galerkin projection, which leads toy_m langle phi_m phi_m rangle = sum_i=1^k sum_j_1=0^n sum_j_2=0^n x_ij_1 x_ij_2 langle phi_j_1 phi_j_2 phi_m rangle quad forall m = 0 dots nHence, knowing the scalars langle phi_m phi_m rangle, and langle phi_j_1 phi_j_2 phi_m rangle, the PCE coefficients y_k can be obtained immediately. From the PCE coefficients we can get the moments and compare them to the closed-form expressions.Notice: A maximum degree of 2 suffices to get the exact solution with PCE. In other words, increasing the maximum degree to values greater than 2 introduces nothing but computational overhead (and numerical errors, possibly)."
},

{
    "location": "chi_squared_k_greater1/#Practice-1",
    "page": "Chi Squared, Several DOFs",
    "title": "Practice",
    "category": "section",
    "text": "First, we create a orthogonal basis relative to f_X(x) of degree at most d=2 (degree below).Notice that we consider a total of Nrec recursion coefficients, and that we also add a quadrature rule by calling OrthoPolyQ().k = 12\nusing PolyChaos\ndegree, Nrec = 2, 20\nop = OrthoPoly(\"gaussian\",degree;Nrec=Nrec);\nopq = OrthoPolyQ(op) #OR: opq = OrthoPolyQ(\"gaussian\",deg;Nrec=Nrec)Now let\'s define a multivariate basismop = MultiOrthoPoly([opq for i=1:k],degree)Next, we define the PCE for all X_i with i = 1 dots k.L = dim(mop)\nmu, sig = 0., 1.\nx = [ assign2multi(convert2affinePCE(\"gaussian\",mu,sig),i,mop.ind) for i=1:k ]With the orthogonal basis and the quadrature at hand, we can compute the tensors t2 and t3 that store the entries langle phi_m phi_m rangle, and langle phi_j_1 phi_j_2 phi_m rangle, respectively.t2 = Tensor(2,mop)\nt3 = Tensor(3,mop)With the tensors at hand, we can compute the Galerkin projection.Notice: there are more efficient ways to do this, but let\'s keep it simple.y = [ sum( x[i][j1]*x[i][j2]*t3.get([j1-1,j2-1,m-1])/t2.get([m-1,m-1])  for i=1:k, j1=1:L, j2=1:L ) for m=1:L ]Let\'s compare the moments via PCE to the closed-form expressions.moms_analytic(k) = [k, sqrt(2k), sqrt(8/k)]\nfunction myskew(y)\n   e3 = sum( y[i]*y[j]*y[k]*t3.get([i-1,j-1,k-1]) for i=1:L,j=1:L,k=1:L )\n   μ = y[1]\n   σ = std(y,mop)\n   (e3-3*μ*σ^2-μ^3)/(σ^3)\nend\n\nprint(\"Expected value:\\t\\t$(moms_analytic(k)[1]) = $(mean(y,mop))\\n\")\nprint(\"\\t\\t\\terror = $(abs(mean(y,mop)-moms_analytic(k)[1]))\\n\")\nprint(\"Standard deviation:\\t$(moms_analytic(k)[2]) = $(std(y,mop))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[2]-std(y,mop))\\n\")\nprint(\"Skewness:\\t\\t$(moms_analytic(k)[3]) = $(myskew(y))\\n\")\nprint(\"\\t\\t\\terror = $(moms_analytic(k)[3]-myskew(y))\\n\")\nLet\'s plot the probability density function to compare results. We first draw samples from the measure with the help of sampleMeasure(), and then evaluate the basis at these samples and multiply times the PCE coefficients. The latter stop is done using evaluatePCE(). Both steps are combined in the function samplePCE(). Finally, we compare the result agains the analytical PDF rho(t) = fract^t2-1mathrme^-t22^k2  Gamma(k2) of the chi-squared distribution with one degree of freedom.using Plots\ngr()\nNsmpl = 10000\n# ξ = sampleMeasure(Nsmpl,mop)\n# ysmpl = evaluatePCE(y,ξ,mop)\nysmpl = samplePCE(Nsmpl,y,mop)\nhistogram(ysmpl;normalize=true,xlabel=\"t\",ylabel=\"rho(t)\")\n\nimport SpecialFunctions: gamma\nρ(t) = 1/(2^(0.5*k)*gamma(0.5*k))*t^(0.5*k-1)*exp(-0.5*t)\nt = range(0.1; stop=maximum(ysmpl), length=100)\nplot!(t,ρ.(t),w=4)"
},

{
    "location": "random_ode/#",
    "page": "Random ODE",
    "title": "Random ODE",
    "category": "page",
    "text": ""
},

{
    "location": "random_ode/#Galerkin-based-Solution-of-Random-Differential-Equation-1",
    "page": "Random ODE",
    "title": "Galerkin-based Solution of Random Differential Equation",
    "category": "section",
    "text": "This tutorial demonstrates how random differential equations can be solved using polynomial chaos expansions (PCE)."
},

{
    "location": "random_ode/#Theory-1",
    "page": "Random ODE",
    "title": "Theory",
    "category": "section",
    "text": "A random differential equation is an ordinary differential equation that has random parameters, hence its solution is itself a (time-varying) random variable. Perhaps the simplest non-trivial example is the following scalar, linear ordinary differential equationdotx(t) = a x(t) quad x(0) = x_0where a is the realization of a Gaussian random variable mathsfa sim mathcalN(mu sigma^2) with mean mu and variance sigma^2. Arguably, for every realization a we can solve the differential equation and obtainx(t) = x_0 mathrme^a tfrom which we find thatln (x(t)) = ln (x_0) + at sim mathcalN(ln(x_0) + mu t (sigma t)^2)In other words, the logarithm of the solution is normally distributed (so-called log-normal distribution).We\'d like to obtain this result numerically with the help of PCE. The first step is to define the (truncated) PCE for the random variable mathsfamathsfa = sum_i=0^L a_i phi_iwhere a_i are the so-called PCE coefficients, and phi_i are the orthogonal basis polynomials. As the solution to the random differential equation is itself a random variable, we treat x(t) as the realization of the random variable mathsfx(t), and define its PCEmathsfx(t) = sum_i=0^L x_i(t) phi_iThe question is how to obtain the unknown PCE coefficients x_i(t) from the known PCE coefficients a_i relative to the orthogonal basis polynomials phi_i. This can be done using Galerkin projection, which is nothing else than projecting onto the orthogonal basis. Think of a three-dimensional space, in which you have placed some three-dimensional object. If you know project the silhouett of the object onto every axis of the three-dimensional space, then you are doing a Galerkin projection. With PCE the concept is equivalent, but the imagination has a harder time. The first step for Galerkin projection is to insert the PCEssum_i=0^L dotx_i(t) phi_i = sum_j=0^L a_j phi_j sum_k=0^L x_k(t) phi_kthe second step is to project onto every basis polynomial phi_m for m = 0 1 dots L, and to exploit orthogonality of the basis. This givesdotx_m(t) langle phi_m phi_m rangle = sum_j=0^L sum_k=0^L a_j x_k(t) langle phi_l phi_k phi_m rangle quad m = 0 1 dots LOf course, the initial condition must not be forgotten:x_0(0) = x_0 quad x_m(0) = 0 quad m = 1 dots LIf we can solve this enlarged system of ordinary random differential equations, we can reconstruct the analytic solution."
},

{
    "location": "random_ode/#Practice-1",
    "page": "Random ODE",
    "title": "Practice",
    "category": "section",
    "text": "We begin by defining the random differential equationx0 = 2.0\nμ, σ = -0.5, 0.05\ntend, Δt = 3.0, 0.01Next, we define an orthogonal basis (and its quadrature rule) relative to the Gaussian measure using PolyChaos. We choose a maximum degree of L.using PolyChaos\nL, Nrec = 6, 40\nopq = OrthoPolyQ(\"gaussian\",L;Nrec=Nrec)Now we can define the PCE for mathsfa and solve the Galerkin-projected ordinary differential equation using DifferentialEquations.jl.using DifferentialEquations\n\na = [ convert2affinePCE(\"gaussian\",μ,σ); zeros(Float64,L-1) ] # PCE coefficients of a\nxinit = [ x0; zeros(Float64,L) ] # PCE coefficients of initial condition\n\nt2 = Tensor(2,opq); # \\langle \\phi_i, \\phi_j \\rangle\nt3 = Tensor(3,opq); # \\langle \\phi_i \\phi_j, \\phi_k \\rangle\n\n# Galerkin-projected random differential equation\nfunction ODEgalerkin(du,u,p,t)\n   du[:] = [ sum( p[j+1]*u[k+1]*t3.get([j,k,m])/t2.get([m,m]) for j=0:L for k=0:L) for m=0:L ]\nend\n\nprobgalerkin = ODEProblem(ODEgalerkin,xinit,(0,tend),a)\nsolgalerkin = solve(probgalerkin;saveat=0:Δt:tend)\nt, x = solgalerkin.t, solgalerkin.u;For later purposes we compute the expected value and the standard deviation at all time instants using PCE.# an advantage of PCE is that moments can be computed from the PCE coefficients alone; no sampling required\nmean_pce = [ mean(x[i],opq) for i=1:length(x)]  \nstd_pce = [ std(x[i],opq) for i=1:length(x) ]We compare the solution from PCE to a Monte-Carlo-based solution. That means to solve the ordinary differential equation for many samples of mathsfa. We first sample from the measure using sampleMeasure, and then generate samples of mathsfa using evaluatePCE. After that we solve the ODE and store the results in xmc.using Statistics\nNsmpl = 5000\nξ = sampleMeasure(Nsmpl,opq)     # sample from Gaussian measure; effectively randn() here    \nasmpl = evaluatePCE(a,ξ,opq)     # sample random variable with PCE coefficients a; effectively μ + σ*randn() here\n# or: asmpl = samplePCE(Nsmpl,a,opq)\nxmc = [ solve(ODEProblem((u,p,t)->aa*u,x0,(0,tend));saveat=0:Δt:tend).u for aa in asmpl]\nxmc = hcat(xmc...);Now we can compare the Monte Carlo mean and standard deviation to the expression from PCE for every time instant.[ mean(xmc,dims=2)-mean_pce std(xmc,dims=2)-std_pce]Clearly, the accuracy of PCE deteriorates over time. Possible remedies are to increase the dimension of PCE, and to tweak the tolerances of the integrator.Finally, we compare whether the samples follow a log-normal distribution, and compare the result to the analytic mean and standard deviation.logx_pce = [ log.(evaluatePCE(x[i],ξ,opq)) for i=1:length(t)]\n[mean.(logx_pce)-(log(x0) .+ μ*t) std.(logx_pce)-σ*t ]"
},

{
    "location": "math/#",
    "page": "Mathematical Background",
    "title": "Mathematical Background",
    "category": "page",
    "text": ""
},

{
    "location": "math/#MathematicalBackground-1",
    "page": "Mathematical Background",
    "title": "Mathematical Background",
    "category": "section",
    "text": "This section is heavily based on the book \"Orthogonal Polynomials: Computation and Approximation\" by Walter Gautschi (Oxford University Press)"
},

{
    "location": "math/#Orthogonal-Polynomials-1",
    "page": "Mathematical Background",
    "title": "Orthogonal Polynomials",
    "category": "section",
    "text": ""
},

{
    "location": "math/#Basic-Theory-1",
    "page": "Mathematical Background",
    "title": "Basic Theory",
    "category": "section",
    "text": "We always work with absolutely continuous measures for which we write mathrmd lambda (t) = w(t) mathrmdt, where the so-called weight function wis a nonnegative integrable function on the real line mathbbR, i.e. w mathcalW subseteq mathbbR rightarrow mathbbR_geq 0\nhas finite limits in case mathcalW = mathbbR, i.e.lim_t to pm infty w(t)  inftyhas finite moments of all ordersmu_r(mathrmdlambda) = int_mathcalW t^r mathrmd lambda (t) quad r = 0 1 2 dots quad textwith mu_0  0For any pair of integrable functions u v, their scalar product relative to mathrmd lambda is defined aslangle u v rangle_mathrmd lambda = int_mathcalW u(t) v(t) mathrmd lambda(t)Let mathcalP be the set of real polynomials and mathcalP_d subset mathcalP be the set of real polynomials of degree at most d on mathcalW, respectively. Monic real polynomials are real polynomials with leading coefficient equal to one, i.e. pi_k(t) = t^k + dots for k = 0 1 dotsThe polynomials uv in mathcalP with u neq v are said to be orthogonal iflangle u v rangle_mathrmd lambda = int_mathcalW u(t) v(t) mathrmd lambda(t) = 0The norm of u is given by u _ mathrmdlambda = sqrtlangle u u rangleIf the polynomials u in mathcalP has unit length  u _ mathrmdlambda = 1, it is called orthonormal.Monic orthogonal polynomials are polynomials that are monic and orthogonal, hence satisfypi_k(t) = pi_k(t mathrmd lambda) = t^k + dots\nfor k = 0 1 dots, and\nlangle pi_k pi_l rangle_mathrmdlambda = 0\nfor k neq l and k l = 0 1 dots, and langle pi_k pi_k rangle_mathrmdlambda =  pi_k ^2_ mathrmdlambda  0 for k = 0 1 dots.note: Note\nThe support mathcalW of mathrmd lambda can be an interval (finite, half-finite, infinite), or a finite number of disjoint intervals. If the support consists of a finite or denumerably infinite number of distinct points t_k at which lambda has positive jumps w_k, the measure is called a discrete measure. For a finite number N of points, the discrete measure is denoted by mathrmdlambda_N, and it is characterized by its nodes and weights  t_k w_k _k=1^N according tomathrmd lambda_N (t) = sum_k=1^N w_k delta(t-t_k)where delta is the delta-function.The inner product associated with mathrmd lambda_N islangle u v rangle_mathrmdlambda_N = int_mathcalW u(t) v(t) mathrmd lambda_N (t) = sum_k=1^N w_k u(t_k) v(t_k)There exist only N orthogonal polynomials  pi_k(mathrmd lambda_N) _k=0^N-1 that are orthogonal relative to the discrete measure mathrmd lambda_N in the senselangle pi_k(mathrmd lambda_N) pi_l(mathrmd lambda_N) rangle_mathrmdlambda_N =  pi_k(mathrmd lambda_N) _mathrmd lambda_N delta_klwhere delta_kl is the Dirac-delta, for kl = 0 1 dots N-1."
},

{
    "location": "math/#Properties-1",
    "page": "Mathematical Background",
    "title": "Properties",
    "category": "section",
    "text": ""
},

{
    "location": "math/#Symmetry-1",
    "page": "Mathematical Background",
    "title": "Symmetry",
    "category": "section",
    "text": "An absolutely continuous measure mathrmd lambda(t) = w(t) mathrmd t is symmetric (with respect to the origin) if its support is mathcalW = -aa for some 0  a leq infty, and if w(-t) = w(t) for all t in mathcalW.Similarly, a discrete measure mathrmd lambda_N (t) = sum_k=1^N w_k delta(t-t_k) is symmetric if t_k = - t_N+1-k, and w_k = w_N+1-k for k=1 2 dots N.Theorem 1.17 states that: If mathrmd lambda is symmetric, thenpi_k(-t mathrmd lambda) = (-1)^k pi_k(t mathrmd lambda) quad k=01 dotshence the parity of k decides whether pi_k is even/odd.note: Note\nSymmetry is exploited in computeSP, where symmetry need not be relative to the origin, but some arbitrary point of the support."
},

{
    "location": "math/#Three-term-Recurrence-Relation-1",
    "page": "Mathematical Background",
    "title": "Three-term Recurrence Relation",
    "category": "section",
    "text": "The fact that orthogonal polynomials can be represented in terms of a three-term recurrence formula is at the heart of all numerical methods of the package. The importance of the three-term recurrence relation is difficult to overestimate. It providesefficient means of evaluating polynomials (and derivatives),\nzeros of orthogonal polynomials by means of a eigenvalues of a symmetric, tridiagonal matrix\nacces to quadrature rules,\nnormalization constants to create orthonormal polynomials.Theorem 1.27 states:Let pi_k(cdot) = pi_k(cdot mathrmdlambda) for k = 0 1 dots be the monic orthogonal polynomials with respect to the measure mathrmd lambda. Thenbeginaligned\npi_k+1(t) = (t - alpha_k) pi_k(t) - beta_k pi_k-1(t) quad k= 0 1 dots \npi_o(t) = 1 \npi_-1(t) = 0\nendalignedwherebeginaligned\nalpha = alpha_k(mathrmd lambda) = fraclangle t pi_k pi_k rangle_mathrmd lambdalangle pi_k pi_k rangle_mathrmd lambda  k=012 dots \nbeta = beta_k(mathrmd lambda) = fraclangle pi_k pi_k rangle_mathrmd lambdalangle pi_k-1 pi_k-1 rangle_mathrmd lambda  k=12dots\nendalignedLet tildepi_k(cdot) = tildepi_k(cdot mathrmd lambda t) denote the orthonormal polynomials, thenbeginaligned\nsqrtbeta_k+1 tildepi_k(t) = (t - alpha_k) tildepi_k(t) - sqrtbeta_k tildepi_k-1(t) quad k = 0 1 dots \ntildepi_0(t) = 1 \ntildepi_-1(t) = 0\nendalignednote: Note\nWithin the package, the coefficients (α,β) are the building block to represent (monic) orthogonal polynomials.Notice that beta_0 is arbitrary. Nevertheless, it is convenient to define it asbeta_0(mathrmdlambda) = langle pi_0 pi_0 rangle_mathrmd lambda = int_mathcalW mathrmd lambda (t)because it allows to compute the norms of the polynomials based on beta_k alone pi_n _mathrmd lambda = beta_n(mathrmd lambda) beta_n-1(mathrmd lambda) cdots beta_0(mathrmd lambda) quad n = 01 dotsLet the support be mathcalW = ab for 0  ab  infty, thenbeginaligned\n a  alpha_k(mathrmd lambda)  b  k = 012 dots \n 0  beta_k(mathrmd lambda)  max(a^2 b^2)  k = 1 2 dots\nendaligned"
},

{
    "location": "math/#Quadrature-Rules-1",
    "page": "Mathematical Background",
    "title": "Quadrature Rules",
    "category": "section",
    "text": "An n-point quadrature rule for the measure mathrmd lambda t is a formula of the formint_mathcalW f(t) mathrmd lambda(t) = sum_nu = 1^n w_nu f(tau_nu) + R_n(f)The quadrature rule  (tau_nu w_nu) _nu=1^n composed of (mutually distinct) nodes tau_nu and weights w_nu provides an approximation to the integral. The respective error is given by R_n(f). If, for polynomials p in mathcalP_d, the error R_n(p) vanishes, the respective quadrature rule is said to have a degree of exactness d. Gauss quadrature rule are special quadrature rules that have a degree of exactness d = 2n - 1. That means, taking a n =3-point quadrature rule, polynomials up to degree 5 can be integrated exactly. The nodes and weights for the Gauss quadrature rules have some remarkable properties:all Gauss nodes are mutually distinct and contained in the interior of the support of mathrmd lambda;\nthe n Gauss nodes are the zeros of pi_n, the monic orthogonal polynomial of degree n relative to the measure mathrmd lambda;\nall Gauss weights are positive.The Gauss nodes and weights can be computed using the Golub-Welsch algorithm. This means to solve an eigenvalue problem of a symmetric tridiagonal matrix."
},

{
    "location": "functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "note: Note\nThe core interface of (so we hope) all essential functions are not dependent on specialized types such as OrthoPoly/OrthoPolyQ. Having said that, for exactly those essential functions there exist overloaded functions that accept specialized types such as OrthoPoly/OrthoPolyQ as arguments.Too abstract? For example, the function evaluate() that evaluates a polynomial of degree n at points x has the core interface    evaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})where a and b are the vectors of recurrence coefficients. For simplicity, there also exists the interface    evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly)which is defined as    evaluate(n::Int64,x::Vector{Float64},op::OrthoPoly) = evaluate(n,x,op.α,op.β)So fret not upon the encounter of multiply-dispatched versions of the same thing. It\'s there to simplify your life (so we hope).The idea of this approach is to make it simpler for others to copy and paste code snippets and use them in their own work.List of all functions in PolyChaos."
},

{
    "location": "functions/#PolyChaos.r_scale-Tuple{Float64,Array{Float64,1},Array{Float64,1}}",
    "page": "Functions",
    "title": "PolyChaos.r_scale",
    "category": "method",
    "text": "r_scale(c::Float64,β::Vector{Float64},α::Vector{Float64})\n\nGiven the recursion coefficients (α,β) for a system of orthogonal polynomials that are orthogonal with respect to some positive weight m(t), this function returns the recursion coefficients (α_,β_) for the scaled measure c m(t) for some positive c.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_compute-Tuple{Function,Float64,Float64}",
    "page": "Functions",
    "title": "PolyChaos.rm_compute",
    "category": "method",
    "text": "rm_compute(weight::Function,lb::Float64,ub::Float64,Npoly::Int64=4,Nquad::Int64=10;quadrature::Function=clenshaw_curtis)\n\nGiven a positive weight function with domain (lb,ub), i.e. a function w lb ub  rightarrow mathbbR_geq 0, this function creates Npoly recursion coefficients (α,β).\n\nThe keyword quadrature specifies what quadrature rule is being used.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_logistic-Tuple{Int64}",
    "page": "Functions",
    "title": "PolyChaos.rm_logistic",
    "category": "method",
    "text": "rm_logistic(N::Int)\n\nCreates N recurrence coefficients for monic polynomials that are orthogonal on (-inftyinfty) relative to w(t) = fracmathrme^-t(1 - mathrme^-t)^2\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_hermite-Tuple{Int64,Float64}",
    "page": "Functions",
    "title": "PolyChaos.rm_hermite",
    "category": "method",
    "text": "rm_hermite(N::Int,mu::Float64)\nrm_hermite(N::Int)\n\nCreates N recurrence coefficients for monic generalized Hermite polynomials that are orthogonal on (-inftyinfty) relative to w(t) = t^2 mu mathrme^-t^2\n\nThe call rm_hermite(N) is the same as rm_hermite(N,0).\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_hermite_prob-Tuple{Int64}",
    "page": "Functions",
    "title": "PolyChaos.rm_hermite_prob",
    "category": "method",
    "text": "rm_hermite_prob(N::Int)\n\nCreates N recurrence coefficients for monic probabilists\' Hermite polynomials that are orthogonal on (-inftyinfty) relative to w(t) = mathrme^-05t^2\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_laguerre-Tuple{Int64,Float64}",
    "page": "Functions",
    "title": "PolyChaos.rm_laguerre",
    "category": "method",
    "text": "rm_laguerre(N::Int,a::Float64)\nrm_laguerre(N::Int)\n\nCreates N recurrence coefficients for monic generalized Laguerre polynomials that are orthogonal on (0infty) relative to w(t) = t^a mathrme^-t.\n\nThe call rm_laguerre(N) is the same as rm_laguerre(N,0).\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_legendre-Tuple{Int64}",
    "page": "Functions",
    "title": "PolyChaos.rm_legendre",
    "category": "method",
    "text": "rm_legendre(N::Int)\n\nCreates N recurrence coefficients for monic Legendre polynomials that are orthogonal on (-11) relative to w(t) = 1.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_legendre01-Tuple{Int64}",
    "page": "Functions",
    "title": "PolyChaos.rm_legendre01",
    "category": "method",
    "text": "rm_legendre01(N::Int)\n\nCreates N recurrence coefficients for monic Legendre polynomials that are orthogonal on (01) relative to w(t) = 1.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_jacobi-Tuple{Int64,Float64,Float64}",
    "page": "Functions",
    "title": "PolyChaos.rm_jacobi",
    "category": "method",
    "text": "rm_jacobi(N::Int,a::Float64,b::Float64)\nrm_jacobi(N::Int,a::Float64)\nrm_jacobi(N::Int)\n\nCreates N recurrence coefficients for monic Jacobi polynomials that are orthogonal on (-11) relative to w(t) = (1-t)^a (1+t)^b.\n\nThe call rm_jacobi(N,a) is the same as rm_jacobi(N,a,a) and rm_jacobi(N) the same as rm_jacobi(N,0,0).\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_jacobi01-Tuple{Int64,Float64,Float64}",
    "page": "Functions",
    "title": "PolyChaos.rm_jacobi01",
    "category": "method",
    "text": "rm_jacobi01(N::Int,a::Float64,b::Float64)\nrm_jacobi01(N::Int,a::Float64)\nrm_jacobi01(N::Int)\n\nCreates N recurrence coefficients for monic Jacobi polynomials that are orthogonal on (01) relative to w(t) = (1-t)^a t^b.\n\nThe call rm_jacobi01(N,a) is the same as rm_jacobi01(N,a,a) and rm_jacobi01(N) the same as rm_jacobi01(N,0,0).\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rm_meixner_pollaczek-Tuple{Int64,Float64,Float64}",
    "page": "Functions",
    "title": "PolyChaos.rm_meixner_pollaczek",
    "category": "method",
    "text": "rm_meixner_pollaczek(N::Int,lambda::Float64,phi::Float64)\nrm_meixner_pollaczek(N::Int,lambda::Float64)\n\nCreates N recurrence coefficients for monic  Meixner-Pollaczek polynomials with parameters λ and ϕ. These are orthogonal on  -inftyinfty relative to the weight function w(t)=(2 pi)^-1 exp(2 phi-pi)t Gamma(lambda+ i t)^2.\n\nThe call rm_meixner_pollaczek(n,lambda) is the same as rm_meixner_pollaczek(n,lambda,pi/2).\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.stieltjes",
    "page": "Functions",
    "title": "PolyChaos.stieltjes",
    "category": "function",
    "text": "stieltjes(N::Int64,nodes_::Vector{Float64},weights_::Vector{Float64};removezeroweights::Bool=true)\n\nStieltjes procedure–-Given the nodes and weights the function generates the firstN recurrence coefficients of the corresponding discrete orthogonal polynomials.\n\nSet the Boolean removezeroweights to true if zero weights should be removed.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.lanczos",
    "page": "Functions",
    "title": "PolyChaos.lanczos",
    "category": "function",
    "text": "lanczos(N::Int64,nodes::Vector{Float64},weights::Vector{Float64};removezeroweights::Bool=true)\n\nLanczos procedure–-given the nodes and weights the function generates the first N recurrence coefficients of the corresponding discrete orthogonal polynomials.\n\nSet the Boolean removezeroweights to true if zero weights should be removed.\n\nThe script is adapted from the routine RKPW in W.B. Gragg and W.J. Harrod, The numerically stable reconstruction of Jacobi matrices from spectral data, Numer. Math. 44 (1984), 317-335.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.mcdiscretization",
    "page": "Functions",
    "title": "PolyChaos.mcdiscretization",
    "category": "function",
    "text": "mcdiscretization(N::Int64,quads::Vector{},discretemeasure::Matrix{Float64}=zeros(0,2);discretization::Function=stieltjes,Nmax::Integer=300,ε::Float64=1e-8,gaussquad::Bool=false)\n\nThis routine returns N recurrence coefficients of the polynomials that are orthogonal relative to a weight function w that is decomposed as a sum of m weights w_i with domains a_ib_i for i=1dotsm,\n\nw(t) = sum_i^m w_i(t) quad textwith  operatornamedom(w_i) = a_i b_i\n\nFor each weight w_i and its domain a_i b_i the function mcdiscretization() expects a quadrature rule of the form     nodes::Vector{Float64}, weights::Vector{Float64} = myquadi(N::Int64) all of which are stacked in the parameter quad     quad = [ myquad1, ..., myquadm ] If the weight function has a discrete part (specified by discretemeasure) it is added on to the discretized continuous weight function.\n\nThe function mcdiscretization() performs a sequence of discretizations of the given weight w(t), each discretization being followed by an application of the Stieltjes or Lanczos procedure (keyword discretization in [stieltjes, lanczos]) to produce approximations to the desired recurrence coefficients. The function applies to each subinterval i an N-point quadrature rule (the ith entry of quad) to discretize the weight function w_i on that subinterval. If the procedure converges to within a prescribed accuracy ε before N reaches its maximum allowed value Nmax. If the function does not converge, the function prompts an error message.\n\nThe keyword gaussquad should be set to true if Gauss quadrature rules are available for all m weights w_i(t) with i = 1 dots m.\n\nFor further information, please see W. Gautschi \"Orthogonal Polynomials: Approximation and Computation\", Section 2.2.4.\n\n\n\n\n\n"
},

{
    "location": "functions/#Recurrence-Coefficients-for-Monic-Orthogonal-Polynomials-1",
    "page": "Functions",
    "title": "Recurrence Coefficients for Monic Orthogonal Polynomials",
    "category": "section",
    "text": "The functions below provide analytic expressions for the recurrence coefficients of common orthogonal polynomials. All of these provide monic orthogonal polynomials relative to the weights.note: Note\nThe number N of recurrence coefficients has to be positive for all functions below.r_scale(c::Float64,a::Vector{Float64},b::Vector{Float64})\nrm_compute(weight::Function,lb::Float64,ub::Float64;Npoly::Int64=4,Nquad::Int64=10,quadrature::Function=clenshaw_curtis)\nrm_logistic(N::Int)\nrm_hermite(N::Int,mu::Float64)\nrm_hermite_prob(N::Int)\nrm_laguerre(N::Int,a::Float64)\nrm_legendre(N::Int)\nrm_legendre01(N::Int)\nrm_jacobi(N::Int,a::Float64,b::Float64)\nrm_jacobi01(N::Int,a::Float64,b::Float64)\nrm_meixner_pollaczek(N::Int,lambda::Float64,phi::Float64)\nstieltjes\nlanczos\nmcdiscretization"
},

{
    "location": "functions/#PolyChaos.showpoly",
    "page": "Functions",
    "title": "PolyChaos.showpoly",
    "category": "function",
    "text": "showpoly(coeffs::Vector{Float64};sym::String,digits::Integer)\nshowpoly(coeffs::Vector{Vector{Float64}};sym::String,digits::Integer)\n\nShow the monic polynomial with coefficients coeffs in a human readable way. They keyword sym sets the name of the variable, and digits controls the number of shown digits.\n\njulia> showpoly([1.2, 2.3, 3.4456])\nx^3 + 3.45x^2 + 2.3x + 1.2\njulia> showpoly([1.2, 2.3, 3.4456], sym=\"t\", digits=2)\nt^3 + 3.45t^2 + 2.3t + 1.2\n\nTailored to types from PolyChaos.jl\n\nshowpoly(d::Integer,α::Vector{Float64},β::Vector{Float64};upto::Bool,sym::String,digits::Integer)\nshowpoly(d::Integer,op::OrthoPoly;upto::Bool=true,sym::String,digits::Integer)\nshowpoly(d::Integer,opq::OrthoPolyQ;upto::Bool=true,sym::String,digits::Integer)\n\nShow the monic orthogonal polynomials with recurrence coefficients (α,β) up to degree d. Setting the keyword upto to false prints the monic polynomial of degree equal to d.\n\njulia> op = OrthoPoly(\"gaussian\",5);\njulia> showpoly(3,op)\n1\nx\nx^2 - 1.0\nx^3 - 3.0x\n\njulia> showpoly(3,op; upto=false)\nx^3 - 3.0x\n\njulia> showpoly(3,op; sym=\"t\")\n1\nt\nt^2 - 1.0\nt^3 - 3.0t\n\nThe following function calls show all orthogonal polynomials given (α,β).\n\nshowpoly(α::Vector{Float64},β::Vector{Float64};sym::String,digits::Integer)\nshowpoly(op::Union{OrthoPoly,OrthoPolyQ};sym::String,digits::Integer)\n\njulia> showpoly(op; sym=\"y\")\n1\ny\ny^2 - 1.0\ny^3 - 3.0y\ny^4 - 6.0y^2 + 3.0\ny^5 - 10.0y^3 + 15.0y\ny^6 - 15.0y^4 + 45.0y^2 - 15.0\n\nThanks @pfitzseb for providing this functionality.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.rec2coeff",
    "page": "Functions",
    "title": "PolyChaos.rec2coeff",
    "category": "function",
    "text": "rec2coeff(deg::Int,a::Vector{Float64},b::Vector{Float64})\nrec2coeff(a,b) = rec2coeff(length(a),a,b)\n\nGet the coefficients of the orthogonal polynomial of degree up to deg specified by its recurrence coefficients (a,b). The function returns the values c_i^(k) from\n\np_k (t) = t^d + sum_i=0^k-1 c_i t^i\n\nwhere k runs from 1 to deg.\n\nThe call rec2coeff(a,b) outputs all possible recurrence coefficients given (a,b).\n\n\n\n\n\n"
},

{
    "location": "functions/#Show-Orthogonal-Polynomials-1",
    "page": "Functions",
    "title": "Show Orthogonal Polynomials",
    "category": "section",
    "text": "To get a human-readable output of the orthognoal polynomials there is the function showpoly()showpolyThis function makes excessive use ofrec2coeff"
},

{
    "location": "functions/#PolyChaos.evaluate",
    "page": "Functions",
    "title": "PolyChaos.evaluate",
    "category": "function",
    "text": "Univariate\n\nevaluate(n::Int64,x::Array{Float64},a::Vector{Float64},b::Vector{Float64})\nevaluate(n::Int64,x::Float64,a::Vector{Float64},b::Vector{Float64})\nevaluate(n::Int64,x::Vector{Float64},op::OrthoPoly)\nevaluate(n::Int64,x::Float64,op::OrthoPoly)\nevaluate(n::Int64,x::Vector{Float64},opq::OrthoPolyQ) = evaluate(n,x,opq.op)\nevaluate(n::Int64,x::Float64,opq::OrthoPolyQ) = evaluate(n,[x],opq.op)\n\nEvaluate the n-th univariate basis polynomial at point(s) x The function is multiply dispatched to facilitate its use with the composite type OrthoPoly\n\nIf several basis polynomials (stored in ns) are to be evaluated at points x, then call\n\nevaluate(ns::Vector{Int64},x::Vector{Float64},op::OrthoPoly) = evaluate(ns,x,op.α,op.β)\nevaluate(ns::Vector{Int64},x::Float64,op::OrthoPoly) = evaluate(ns,[x],op)\nevaluate(ns::Vector{Int64},x::Vector{Float64},opq::OrthoPolyQ) = evaluate(ns,x,opq.op)\nevaluate(ns::Vector{Int64},x::Float64,opq::OrthoPolyQ) = evaluate(ns,[x],opq.op)\n\nIf all basis polynomials are to be evaluated at points x, then call\n\nevaluate(x::Vector{Float64},op::OrthoPoly) = evaluate(collect(0:op.deg),x,op)\nevaluate(x::Float64,op::OrthoPoly) = evaluate([x],op)\nevaluate(x::Vector{Float64},opq::OrthoPolyQ) = evaluate(x,opq.op)\nevaluate(x::Float64,opq::OrthoPolyQ) = evaluate([x],opq)\n\nwhich returns an Array of dimensions (length(x),op.deg+1).\n\nnote: Note\nn is the degree of the univariate basis polynomial\nlength(x) = N, where N is the number of points\n(a,b) are the recursion coefficients\n\nMultivariate\n\nevaluate(n::Vector{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})\nevaluate(n::Vector{Int64},x::Vector{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})\nevaluate(n::Vector{Int64},x::Matrix{Float64},op::MultiOrthoPoly)\nevaluate(n::Vector{Int64},x::Vector{Float64},op::MultiOrthoPoly)\n\nEvaluate the n-th p-variate basis polynomial at point(s) x The function is multiply dispatched to facilitate its use with the composite type MultiOrthoPoly\n\nIf several basis polynomials are to be evaluated at points x, then call\n\nevaluate(ind::Matrix{Int64},x::Matrix{Float64},a::Vector{Vector{Float64}},b::Vector{Vector{Float64}})\nevaluate(ind::Matrix{Int64},x::Matrix{Float64},op::MultiOrthoPoly)\n\nwhere ind is a matrix of multi-indices.\n\nIf all basis polynomials are to be evaluated at points x, then call\n\nevaluate(x::Matrix{Float64},mop::MultiOrthoPoly) = evaluate(mop.ind,x,mop)\n\nwhich returns an array of dimensions (mop.dim,size(x,1)).\n\nnote: Note\nn is a multi-index\nlength(n) == p, i.e. a p-variate basis polynomial\nsize(x) = (N,p), where N is the number of points\nsize(a)==size(b)=p.\n\n\n\n\n\n"
},

{
    "location": "functions/#Evaluate-Orthogonal-Polynomials-1",
    "page": "Functions",
    "title": "Evaluate Orthogonal Polynomials",
    "category": "section",
    "text": "evaluate"
},

{
    "location": "functions/#PolyChaos.computeSP2",
    "page": "Functions",
    "title": "PolyChaos.computeSP2",
    "category": "function",
    "text": "computeSP2(n::Int64,β::Vector{Float64})\ncomputeSP2(n::Int64,op::OrthoPoly) = computeSP2(n,op.β)\ncomputeSP2(op::OrthoPoly) = computeSP2(op.deg,op.β)\ncomputeSP2(opq::OrthoPolyQ) = computeSP2(opq.op)\n\nComputes the n regular scalar products aka 2-norms of the orthogonal polynomials, namely\n\nϕ_i^2 = langle phi_iphi_irangle quad forall i in  0dotsn \n\nNotice that only the values of β of the recurrence coefficients (α,β) are required. The computation is based on equation (1.3.7) from Gautschi, W. \"Orthogonal Polynomials: Computation and Approximation\". Whenever there exists an analytic expressions for β, this function should be used.\n\nThe function is multiply dispatched to facilitate its use with the composite types OrthoPoly and OrthoPolyQ.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.computeSP",
    "page": "Functions",
    "title": "PolyChaos.computeSP",
    "category": "function",
    "text": "Univariate\n\ncomputeSP(a::Vector{Int64},α::Vector{Float64},β::Vector{Float64},nodes::Vector{Float64},weights::Vector{Float64};issymmetric::Bool=false)\ncomputeSP(a::Vector{Int64},op::OrthoPoly,q::Quad;issymmetric=issymmetric(op))\ncomputeSP(a::Vector{Int64},opq::OrthoPolyQ)\n\nMultivariate\n\ncomputeSP( a::Vector{Int64},\n           α::Vector{Vector{Float64}},β::Vector{Vector{Float64}},\n           nodes::Vector{Vector{Float64}},weights::Vector{Vector{Float64}},\n           ind::Matrix{Int64};\n           issymmetric::BitArray=falses(length(α)))\ncomputeSP(a::Vector{Int64},op::Vector{OrthoPoly},q::Vector{Quad},ind::Matrix{Int64})\ncomputeSP(a::Vector{Int64},mOP::MultiOrthoPoly)\n\nComputes the scalar product\n\nlangle phi_a_1phi_a_2cdotsphi_a_n rangle\n\nwhere n = length(a). This requires to provide the recurrence coefficients (α,β) and the quadrature rule (nodes,weights), as well as the multi-index ind. If provided via the keyword issymmetric, symmetry of the weight function is exploited. All computations of the multivariate scalar products resort back to efficient computations of the univariate scalar products. Mathematically, this follows from Fubini\'s theorem.\n\nThe function is multiply dispatched to facilitate its use with OrthoPolyQ or a suitable combination of OrthoPoly and its quadrature rule Quad.\n\nnote: Note\nZero entries of a are removed automatically to simplify computations, which follows fromlangle phi_i phi_j phi_0cdotsphi_0 rangle = langle phi_i phi_j ranglebecause \\phi_0 = 1.It is checked whether enough quadrature points are supplied to solve the integral exactly.\n\n\n\n\n\n"
},

{
    "location": "functions/#Scalar-Products-1",
    "page": "Functions",
    "title": "Scalar Products",
    "category": "section",
    "text": "computeSP2\ncomputeSP"
},

{
    "location": "functions/#Quadrature-Rules-1",
    "page": "Functions",
    "title": "Quadrature Rules",
    "category": "section",
    "text": "fejer\nfejer2\nclenshaw_curtis\nquadgp\ngauss\nradau\nradau_jacobi\nradau_laguerre\nlobatto\nlobatto_jacobi"
},

{
    "location": "functions/#Statistics.mean",
    "page": "Functions",
    "title": "Statistics.mean",
    "category": "function",
    "text": "Univariate\n\nmean(x::Vector{Float64},op::OrthoPoly)\nmean(x::Vector{Float64},opq::OrthoPolyQ)\n\nMultivariate\n\nmean(x::Vector{Float64},mop::MultiOrthoPoly)\n\ncompute mean of random variable with PCE x\n\n\n\n\n\n"
},

{
    "location": "functions/#Statistics.var",
    "page": "Functions",
    "title": "Statistics.var",
    "category": "function",
    "text": "Univariate\n\nvar(x::Vector{Float64},op::OrthoPoly)\nvar(x::Vector{Float64},opq::OrthoPolyQ)\n\nMultivariate\n\nvar(x::Vector{Float64},mop::MultiOrthoPoly)\n\ncompute variance of random variable with PCE x\n\n\n\n\n\n"
},

{
    "location": "functions/#Statistics.std",
    "page": "Functions",
    "title": "Statistics.std",
    "category": "function",
    "text": "Univariate\n\nstd(x::Vector{Float64},op::OrthoPoly)\nstd(x::Vector{Float64},opq::OrthoPolyQ)\n\nMultivariate\n\nstd(x::Vector{Float64},mop::MultiOrthoPoly)\n\ncompute standard deviation of random variable with PCE x\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.sampleMeasure",
    "page": "Functions",
    "title": "PolyChaos.sampleMeasure",
    "category": "function",
    "text": "Univariate\n\nsampleMeasure(n::Int64,name::String,w::Function,dom::Tuple{Float64,Float64},symm::Bool,d::Dict;method::String=\"adaptiverejection\")\nsampleMeasure(n::Int64,m::Measure;method::String=\"adaptiverejection\")\nsampleMeasure(n::Int64,op::OrthoPoly;method::String=\"adaptiverejection\")\nsampleMeasure(n::Int64,opq::OrthoPolyQ;method::String=\"adaptiverejection\")\n\nDraw n samples from the measure m described by its\n\nname\nweight function w,\ndomain dom,\nsymmetry property symm,\nand, if applicable, parameters stored in the dictionary d.\n\nBy default an adaptive rejection sampling method is used (from AdaptiveRejectionSampling.jl), unless it is a common random variable for which Distributions.jl is used.\n\nThe function is multiply dispatched to accept OrthoPoly or OrthoPolyQ.\n\nMultivariate\n\nsampleMeasure(n::Int64,m::MultiMeasure;method::Vector{String}=[\"adaptiverejection\" for i=1:length(m.name)])\nsampleMeasure(n::Int64,mop::MultiOrthoPoly;method::Vector{String}=[\"adaptiverejection\" for i=1:length(mop.meas.name)])\n\nMultivariate extension which provides array of samples with n rows and as many columns as the multimeasure has univariate measures.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.evaluatePCE",
    "page": "Functions",
    "title": "PolyChaos.evaluatePCE",
    "category": "function",
    "text": "evaluatePCE(x::Vector{Float64},ξ::Vector{Float64},α::Vector{Float64},β::Vector{Float64})\n\nEvaluation of polynomial chaos expansion\n\nmathsfx = sum_i=0^L x_i phi_ixi_j\n\nwhere L+1 = length(x) and x_j is the jth sample where j=1dotsm with m = length(ξ).\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.samplePCE",
    "page": "Functions",
    "title": "PolyChaos.samplePCE",
    "category": "function",
    "text": "Univariate\n\nsamplePCE(n::Int64,x::Vector{Float64},op::OrthoPoly;method::String=\"adaptiverejection\")\nsamplePCE(n::Int64,x::Vector{Float64},opq::OrthoPolyQ;method::String=\"adaptiverejection\")\n\nCombines sampleMeasure and evaluatePCE, i.e. it first draws n samples from the measure, then evaluates the PCE for those samples.\n\nMultivariate\n\nsamplePCE(n::Int64,x::Vector{Float64},mop::MultiOrthoPoly;method::Vector{String}=[\"adaptiverejection\" for i=1:length(mop.meas.name)])\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.calculateAffinePCE",
    "page": "Functions",
    "title": "PolyChaos.calculateAffinePCE",
    "category": "function",
    "text": "calculateAffinePCE(α::Vector{Float64})::Vector{Float64}\n\nComputes the affine PCE coefficients x_0 and x_1 from recurrence coefficients lpha.\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.convert2affinePCE",
    "page": "Functions",
    "title": "PolyChaos.convert2affinePCE",
    "category": "function",
    "text": "convert2affinePCE(a::Vector{Float64},α0::Float64)\nconvert2affinePCE(name::String,p1::Float64,p2::Float64,d::Dict=Dict();kind::Symbol=:lbub)\n\nComputes the affine PCE coefficients x_0 and x_1 from\n\nX = a_1 + a_2 Xi = x_0 + x_1 phi_1(Xi)\n\nwhere phi_1(t) = t-alpha_0 is the first-order monic basis polynomial.\n\nFor classical polynomials the name can be given directly. The keyword kind in [:lbub, :μσ] specifies whether p1 and p2 have the meaning of lower/upper bounds or mean/standard deviation.\n\n\n\n\n\n"
},

{
    "location": "functions/#Polynomial-Chaos-1",
    "page": "Functions",
    "title": "Polynomial Chaos",
    "category": "section",
    "text": "mean\nvar\nstd\nsampleMeasure\nevaluatePCE\nsamplePCE\ncalculateAffinePCE\nconvert2affinePCE"
},

{
    "location": "functions/#PolyChaos.nw",
    "page": "Functions",
    "title": "PolyChaos.nw",
    "category": "function",
    "text": "nw(q::Quad)\nnw(opq::OrthoPolyQ)\nnw(opq::Vector{OrthoPolyQ})\nnw(mOP::MultiOrthoPoly)\n\nreturns nodes and weights in matrix form\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.coeffs",
    "page": "Functions",
    "title": "PolyChaos.coeffs",
    "category": "function",
    "text": "coeffs(op::OrthoPoly)\ncoeffs(opq::OrthoPolyQ)\ncoeffs(op::Vector{OrthoPoly})\ncoeffs(opq::Vector{OrthoPolyQ})\ncoeffs(mop::MultiOrthoPoly)\n\nreturns recurrence coefficients of in matrix form\n\n\n\n\n\n"
},

{
    "location": "functions/#PolyChaos.integrate",
    "page": "Functions",
    "title": "PolyChaos.integrate",
    "category": "function",
    "text": "integrate(f::Function,nodes::Vector{Float64},weights::Vector{Float64})\nintegrate(f::Function,q::Quad)\nintegrate(f::Function,opq::OrthogonalPolyQ)\n\nintegrate function f using quadrature rule specified via nodes, weights. For example int_0^1 6x^5 = 1 can be solved as follows:\n\njulia> opq = OrthoPolyQ(\"uniform01\",3)\njulia> integrate(x->6x^5,opq)\n1.0000000000000002\n\nnote: Note\n\n\nfunction f is assumed to return a scalar.\ninterval of integration is \"hidden\" in nodes.\n\n\n\n\n\n"
},

{
    "location": "functions/#LinearAlgebra.issymmetric",
    "page": "Functions",
    "title": "LinearAlgebra.issymmetric",
    "category": "function",
    "text": "issymmetric(m::Measure)::Bool\nissymmetric(op::OrthoPoly)::Bool\nissymmetric(q::Quad)::Bool\nissymmetric(opq::OrthoPolyQ)::Bool\n\nIs the measure symmetric (around any point in the domain)?\n\n\n\n\n\n"
},

{
    "location": "functions/#Auxiliary-Functions-1",
    "page": "Functions",
    "title": "Auxiliary Functions",
    "category": "section",
    "text": "nw\ncoeffs\nintegrate\nPolyChaos.issymmetric"
},

]}

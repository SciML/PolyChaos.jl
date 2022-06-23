export EmptyQuad,
       Quad

struct Quad{T, V <: AbstractVector{<:T}} <: AbstractQuad{T}
    name::String
    Nquad::Int # number of qudrature points
    nodes::V
    weights::V

    function Quad(name::String, N::Int, nodes, weights)
        N <= 0 && throw(DomainError(N, "number of qudrature points has to be positive"))
        !(length(nodes) == length(weights)) &&
            throw(InconsistencyError("inconsistent numbers of nodes and weights"))
        new{promote_type(eltype(nodes), eltype(weights)),
            promote_type(typeof(nodes), typeof(weights))}(lowercase(name), N, nodes,
                                                          weights)
    end
end

struct EmptyQuad{T} <: AbstractQuad{T}
    EmptyQuad() = new{Float64}()
end

# general constructor
function Quad(N::Int, α::AbstractVector{<:Real}, β::AbstractVector{<:Real})
    !(length(α) == length(β)) &&
        throw(InconsistencyError("inconsistent numbers of recurrence coefficients"))
    !(N <= length(α) - 1) && throw(DomainError(N,
                      "requested number of quadrature points $N cannot be provided with $(length(α)) recurrence coefficients"))
    nodes, weights = gauss(N, α, β)
    Quad("golubwelsch", N, nodes, weights)
end

# quadrature rules for orthoPolys
Quad(N::Int, op::AbstractOrthoPoly) = Quad(N, op.α, op.β)
Quad(op::AbstractOrthoPoly) = Quad(op.deg, op)

#####################################################
#####################################################
#####################################################
# the constructor below is probably not relevant
#####################################################
#####################################################
#####################################################
# function Quad(N::Int, weight::Function, α::Vector{<:Real}, β::Vector{<:Real}, supp::Tuple{<:Real,<:Real}, symm::Bool, d::Dict=Dict())
#     m = Measure("fun_"*String(nameof(weight)),weight,supp,symm,d)
#     Quad(N,α,β,m)
# end
#####################################################
#####################################################
#####################################################

# all-purpose constructor (last resort!)
function Quad(N::Int, w::Function, dom::Tuple{<:Real, <:Real};
              quadrature::Function = clenshaw_curtis)
    N <= 0 && throw(DomainError(N, "number of quadrature points has to be positive"))
    nodes, weights = quadgp(w, dom[1], dom[2], N; quadrature = quadrature)
    Quad("quadgp", N, nodes, weights)
end

function Quad(N::Int, measure::AbstractMeasure; quadrature::Function = clenshaw_curtis)
    typeof(measure) != Measure &&
        @warn "For measures of type $(typeof(measure)) the quadrature rule should be based on the recurrence coefficients."
    Quad(N, measure.w, (measure.dom[1], measure.dom[2]); quadrature = quadrature)
end

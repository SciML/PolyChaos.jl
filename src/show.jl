export rec2coeff, showpoly, showbasis

title_color = :underline
name_color = :light_blue

function show(measure::AbstractMeasure)
    printstyled("Measure dλ(t)=w(t)dt\n"; color = title_color)
    for name in fieldnames(typeof(measure))
        println("$(string(name)):\t$(getfield(measure, name))")
    end
end

function show(productmeasure::ProductMeasure)
    for measure in productmeasure.measures
        show(measure)
    end
end

function show(op::AbstractOrthoPoly; showmeasure::Bool = true)
    printstyled("\nUnivariate orthogonal polynomials\n"; color = title_color)
    # print("name:"); printstyled("\t\t$(op.name)\n";color=name_color)
    print("degree:")
    printstyled("\t\t$(op.deg)\n")
    print("#coeffs:")
    printstyled("\t$(length(op.α))\n")
    length(op.α) <= 7 ? n = length(op.α) : n = 7
    print("α =")
    printstyled("\t\t$(op.α[1:n])")
    n < length(op.α) ? print("...\n") : print("\n")
    print("β =")
    printstyled("\t\t$(op.β[1:n])")
    n < length(op.α) ? print("...\n") : print("\n")
    if showmeasure
        print("\n")
        show(op.measure)
    end
end

function show(mop::MultiOrthoPoly; showmeasure::Bool = false)
    p = length(mop.name)
    printstyled("\n$p-variate orthogonal polynomials\n"; color = title_color)
    print("name:")
    [printstyled("\t\t$(mop.name[i])\n"; color = name_color) for i in 1:p]
    print("deg:")
    printstyled("\t\t$(mop.deg)\n")
    print("dim:")
    printstyled("\t\t$(mop.dim)\n")
    print("ind:")
    mop.dim <= 7 ? n = mop.dim : n = 7
    [printstyled("\t\t$(mop.ind[i,:])\n") for i in 1:n]
    n < mop.dim ? print("\t\t...\n") : print("\n")
    print("\n")
    showmeasure && show(mop.measure)
end

function show(t::Tensor)
    printstyled("\n$(t.dim)-dimensional tensor\n"; color = title_color)
    print("dim:")
    printstyled("\t\t$(t.dim)\n")
    print("nonzeros:")
    printstyled("\t$(length(t.T.nzind))\n")
    show(t.op; showmeasure = false)
end

"""
```
rec2coeff(deg::Int,a::Vector{<:Real},b::Vector{<:Real})
rec2coeff(a,b) = rec2coeff(length(a),a,b)
```

Get the coefficients of the orthogonal polynomial of degree up to `deg` specified by its
recurrence coefficients `(a,b)`. The function returns the values ``c_i^{(k)}`` from

```math
p_k (t) = t^d + \\sum_{i=0}^{k-1} c_i t^i,
```

where ``k`` runs from `1` to `deg`.

The call `rec2coeff(a,b)` outputs all possible recurrence coefficients given `(a,b)`.
"""
function rec2coeff(deg::Int, a::Vector{<:Real}, b::Vector{<:Real})
    deg <= 0 &&
        throw(DomainError(deg, "degree must be positive (you asked for degree = $deg"))
    !(length(a) == length(b) && length(a) >= deg) &&
        throw(InconsistencyError("incorrect number of recurrence coefficients"))
    c = Vector{Vector{<:Real}}(undef, deg)
    @inbounds c[1] = [-a[1]]
    deg == 1 && return c
    @inbounds c[2] = [a[1] * a[2] - b[2], -a[1] - a[2]]
    deg == 2 && return c
    for k in 2:(deg - 1)
        c[k + 1] = [-a[k + 1] * c[k][1] - b[k + 1] * c[k - 1][1]]
        for i in 1:(k - 2)
            push!(c[k + 1], c[k][i] - a[k + 1] * c[k][i + 1] - b[k + 1] * c[k - 1][i + 1])
        end
        push!(c[k + 1], c[k][k - 1] - a[k + 1] * c[k][k] - b[k + 1])
        push!(c[k + 1], c[k][k] - a[k + 1])
    end
    return c
end
rec2coeff(α, β) = rec2coeff(length(α), α, β)

"""
```
showpoly(coeffs::Vector{<:Real};sym::String,digits::Integer)
```

Show the monic polynomial with coefficients `coeffs` in a human-readable way.
The keyword `sym` sets the name of the variable, and `digits` controls the number of shown digits.

```jldoctest
julia> using PolyChaos

julia> showpoly([1.2, 2.3, 3.4456])
x^3 + 3.45x^2 + 2.3x + 1.2

julia> showpoly([1.2, 2.3, 3.4456], sym = "t", digits = 2)
t^3 + 3.45t^2 + 2.3t + 1.2
```

```
showpoly(d::Integer,α::Vector{<:Real},β::Vector{<:Real}; sym::String,digits::Integer)
showpoly(d::Range,α::Vector{<:Real},β::Vector{<:Real};sym::String,digits::Integer) where Range <: OrdinalRange
```

Show the monic polynomial of degree/range `d` that has the recurrence coefficients `α`, `β`.

```jldoctest
julia> using PolyChaos

julia> α, β = rm_hermite(10)
([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.77245, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5])

julia> showpoly(3, α, β)
x^3 - 1.5x

julia> showpoly(0:2:10, α, β)
1
x^2 - 0.5
x^4 - 3.0x^2 + 0.75
x^6 - 7.5x^4 + 11.25x^2 - 1.88
x^8 - 14.0x^6 + 52.5x^4 - 52.5x^2 + 6.56
x^10 - 22.5x^8 + 157.5x^6 - 393.75x^4 + 295.31x^2 - 29.53
```

Tailored to types from `PolyChaos.jl`

```
showpoly(d::Union{Integer,Range},op::AbstractOrthoPoly;sym::String,digits::Integer) where Range <: OrdinalRange
```

Show the monic polynomial of degree/range `d` of an `AbstractOrthoPoly`.

```jldoctest
julia> using PolyChaos

julia> op = GaussOrthoPoly(5);

julia> showpoly(3, op)
x^3 - 3.0x

julia> showpoly(0:(op.deg), op; sym = "t")
1
t
t^2 - 1.0
t^3 - 3.0t
t^4 - 6.0t^2 + 3.0
t^5 - 10.0t^3 + 15.0t
```

[Thanks @pfitzseb for providing this functionality.](https://discourse.julialang.org/t/how-to-define-verbose-output-for-a-polynomial/21317/5)
"""
function showpoly(coeffs::Vector{<:Real}; sym::String = "x", digits::Integer = 2)
    io = Base.stdout
    length(coeffs) > 1 ? print(io, sym * "^", length(coeffs)) : print(io, sym)
    for (i, c) in enumerate(reverse(coeffs))
        abs(round(c, digits = digits)) == 0.0 && continue
        ex = length(coeffs) - i
        print(io, ' ', c > 0 ? '+' : '-', ' ')
        print(io, abs(round(c, digits = digits)))
        ex > 0 && print(io, sym)
        ex > 1 && print(io, '^', ex)
    end
    print(io, '\n')
end

function showpoly(d::Integer, α::Vector{<:Real}, β::Vector{<:Real}; sym::String = "x",
                  digits::Integer = 2)
    @assert d>=0 "degree has to be non-negative."
    d == 0 && return print("1\n")
    showpoly(rec2coeff(d, α, β)[end], sym = sym, digits = digits)
end

function showpoly(d::Range, α::Vector{<:Real}, β::Vector{<:Real}; sym::String = "x",
                  digits::Integer = 2) where {Range <: OrdinalRange}
    map(c -> showpoly(c, α, β; sym = sym, digits = digits), d)
    print()
end

function showpoly(d::Union{Integer, Range}, op::AbstractOrthoPoly; sym::String = "x",
                  digits::Integer = 2) where {Range <: OrdinalRange}
    showpoly(d, op.α, op.β; sym = sym, digits = digits)
end

"""
```
showbasis(α::Vector{<:Real},β::Vector{<:Real};sym::String,digits::Integer)
```

Show all basis polynomials given the recurrence coefficients `α`, `β`.
The keyword `sym` sets the name of the variable, and `digits` controls the number of shown digits.

```jldoctest
julia> using PolyChaos

julia> α, β = rm_hermite(5);

julia> showbasis(α, β)
1
x
x^2 - 0.5
x^3 - 1.5x
x^4 - 3.0x^2 + 0.75
x^5 - 5.0x^3 + 3.75x
```

Tailored to types from `PolyChaos.jl`

```
showbasis(op::AbstractOrthoPoly;sym::String,digits::Integer)
```

Show all basis polynomials of an `AbstractOrthoPoly`.

```jldoctest
julia> using PolyChaos

julia> op = LegendreOrthoPoly(4);

julia> showbasis(op)
1
x
x^2 - 0.33
x^3 - 0.6x
x^4 - 0.86x^2 + 0.09
```
"""
function showbasis(α::Vector{<:Real}, β::Vector{<:Real}; sym::String = "x",
                   digits::Integer = 2)
    showpoly(0:length(α), α, β; sym = sym, digits = digits)
end
function showbasis(op::AbstractOrthoPoly; sym::String = "x", digits::Integer = 2)
    showbasis(op.α[1:(op.deg)], op.β[1:(op.deg)]; sym = sym, digits = digits)
end

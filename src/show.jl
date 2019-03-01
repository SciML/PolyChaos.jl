export rec2coeff, showpoly

title_color =   :underline
name_color  =   :light_blue

function show(io::IO,m::Measure)
    printstyled("Measure dλ(t)=w(t)dt\n";color=title_color)
    print("name:\t\t")
    printstyled("$(m.name)\n"; color=name_color)
    println("w(t):\t\t$(m.w)")
    println("dom:\t\t($(m.dom[1]), $(m.dom[2]))")
    print("symmetric:\t")
    m.symmetric ? col=:green : col=:red
    printstyled("$(m.symmetric)\n";color=col)
    if typeof(m.pars)!=Dict{Any,Any}
        print("pars:")
        for (i,d) in enumerate(m.pars)
            printstyled("\t\t$d\n")
        end
    end
end

function show(io::IO,m::MultiMeasure)
    p = length(m.name)
    printstyled("$p-variate measure dλ(t)= ∏_{i=1}^$p w_i(t)dt_i\n";color=title_color)
    print("name:")
    [ printstyled("\t\t$(m.name[i])\n"; color=name_color) for i = 1:p]
    println("w(t):\t\t∏_{i=1}^$p w_i(t)")
    print("dom:")
    [ printstyled("\t\t($(m.dom[i][1]), $(m.dom[i][2]))\n") for i=1:p ]
    print("symmetric:")
    for i=1:p
        m.symmetric[i] ? col=:green : col=:red
        i==1 ? printstyled("\t$(m.symmetric[i])\n";color=col) : printstyled("\t\t$(m.symmetric[i])\n";color=col)
    end
    if length(findall(x->typeof(x)!=Dict{Any,Any},m.pars))>=1 print("pars:")
        for i=1:p
            print("\t\t")
            if typeof(m.pars[i])!=Dict{Any,Any}
                for (i,d) in enumerate(m.pars[i])
                    i==1 ? print("$d") : print(", $d")
                end
                print("\n")
            else
                print("∅\n")
            end
        end
    end
end

function show(io::IO,op::OrthoPoly;showmeasure::Bool=true)
    printstyled("\nUnivariate orthogonal polynomials\n";color=title_color)
    print("name:"); printstyled("\t\t$(op.name)\n";color=name_color)
    print("degree:"); printstyled("\t\t$(op.deg)\n")
    print("#coeffs:"); printstyled("\t$(length(op.α))\n")
    length(op.α)<=7 ? n=length(op.α) : n=7
    print("α ="); printstyled("\t\t$(op.α[1:n])")
    n<length(op.α) ? print("...\n") : print("\n")
    print("β ="); printstyled("\t\t$(op.β[1:n])")
    n<length(op.α) ? print("...\n") : print("\n")
    if showmeasure
        print("\n")
        show(io,op.meas)
    end
end

function show(io::IO,q::Quad;showmeasure::Bool=true)
    printstyled("\nQuadrature rule\n";color=title_color)
    print("name:"); printstyled("\t\t$(q.name)\n";color=name_color)
    print("N:"); printstyled("\t\t$(q.Nquad)\n")
    q.Nquad<=7 ? n=q.Nquad : n=7
    print("nodes"); printstyled("\t\t$(q.nodes[1:n])")
    n<q.Nquad ? print("...\n") : print("\n")
    print("weights"); printstyled("\t\t$(q.weights[1:n])")
    n<q.Nquad ? print("...\n") : print("\n")
    if showmeasure
        print("\n")
        show(io,q.meas)
    end
end

function show(io::IO,opq::OrthoPolyQ;showmeasure::Bool=false)
    show(io,opq.op)
    show(io,opq.quad;showmeasure=showmeasure)
end

function show(io::IO,mop::MultiOrthoPoly;showmeasure::Bool=false)
    p = length(mop.name)
    printstyled("\n$p-variate orthogonal polynomials\n";color=title_color)
    print("name:");
    [ printstyled("\t\t$(mop.name[i])\n";color=name_color) for i=1:p ]
    print("deg:"); printstyled("\t\t$(mop.deg)\n")
    print("dim:"); printstyled("\t\t$(mop.dim)\n")
    print("ind:");
    mop.dim<=7 ? n=mop.dim : n=7
    # printstyled("\t$(mop.ind[1,:])\n")
    [ printstyled("\t\t$(mop.ind[i,:])\n") for i=1:n ]
    n<mop.dim ? print("\t\t...\n") : print("\n")

    print("\n")
    showmeasure ? show(io,mop.meas) : ()
end

function show(io::IO,t::Tensor)
    printstyled("\n$(t.dim)-dimensional tensor\n";color=title_color)
    print("dim:"); printstyled("\t\t$(t.dim)\n")
    print("nonzeros:"); printstyled("\t$(length(t.T.nzind))\n")
    show(io,t.op;showmeasure=false)
end

"""
```
rec2coeff(deg::Int,a::Vector{Float64},b::Vector{Float64})
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
function rec2coeff(deg::Int,a::Vector{Float64},b::Vector{Float64})
    @assert deg >= 1 "degree must be positive (you asked for degree = $deg)"
    @assert length(a) == length(b) && length(a) >= deg "incorrect number of recurrence coefficients"
    c = Vector{Vector{Float64}}(undef,deg)
    @inbounds c[1] = [ -a[1] ]
    deg == 1 && return c
    @inbounds c[2] = [ a[1]*a[2] - b[2], -a[1] - a[2] ]
    deg == 2 && return c
    for k in 2:deg-1
        c[k+1] = [ -a[k+1]*c[k][1] - b[k+1]*c[k-1][1] ]
        for i = 1:k-2
            push!(c[k+1],c[k][i] - a[k+1]*c[k][i+1] - b[k+1]*c[k-1][i+1])
        end
        push!(c[k+1],c[k][k-1] - a[k+1]*c[k][k] - b[k+1])
        push!(c[k+1],c[k][k] - a[k+1])
    end
    return c
end
rec2coeff(α,β) = rec2coeff(length(α),α,β)

"""
```
showpoly(coeffs::Vector{Float64};sym::String,digits::Integer)
showpoly(coeffs::Vector{Vector{Float64}};sym::String,digits::Integer)
```
Show the monic polynomial with coefficients `coeffs` in a human readable way.
They keyword `sym` sets the name of the variable, and `digits` controls the number of shown digits.
```@repl
julia> showpoly([1.2, 2.3, 3.4456])
x^3 + 3.45x^2 + 2.3x + 1.2
julia> showpoly([1.2, 2.3, 3.4456], sym="t", digits=2)
t^3 + 3.45t^2 + 2.3t + 1.2
```

Tailored to types from `PolyChaos.jl`
```
showpoly(d::Integer,α::Vector{Float64},β::Vector{Float64};upto::Bool,sym::String,digits::Integer)
showpoly(d::Integer,op::OrthoPoly;upto::Bool=true,sym::String,digits::Integer)
showpoly(d::Integer,opq::OrthoPolyQ;upto::Bool=true,sym::String,digits::Integer)
```
Show the monic orthogonal polynomials with recurrence coefficients `(α,β)` up to degree `d`.
Setting the keyword `upto` to false prints the monic polynomial of degree equal to `d`.
```@repl
julia> op = OrthoPoly("gaussian",5);
julia> showpoly(3,op)
1
x
x^2 - 1.0
x^3 - 3.0x

julia> showpoly(3,op; upto=false)
x^3 - 3.0x

julia> showpoly(3,op; sym="t")
1
t
t^2 - 1.0
t^3 - 3.0t
```

The following function calls show all orthogonal polynomials given `(α,β)`.
```
showpoly(α::Vector{Float64},β::Vector{Float64};sym::String,digits::Integer)
showpoly(op::Union{OrthoPoly,OrthoPolyQ};sym::String,digits::Integer)
```

```@repl
julia> showpoly(op; sym="y")
1
y
y^2 - 1.0
y^3 - 3.0y
y^4 - 6.0y^2 + 3.0
y^5 - 10.0y^3 + 15.0y
y^6 - 15.0y^4 + 45.0y^2 - 15.0
```

[Thanks @pfitzseb for providing this functionality.](https://discourse.julialang.org/t/how-to-define-verbose-output-for-a-polynomial/21317/5)

"""
function showpoly(coeffs::Vector{Float64};sym::String="x",digits::Integer=2)
    io = Base.stdout
    length(coeffs) > 1 ? print(io, sym*"^", length(coeffs)) : print(io, sym)
    for (i, c) in enumerate(reverse(coeffs))
        abs(round(c, digits=digits)) == 0. && continue
        ex = length(coeffs) - i
        print(io, ' ', c > 0 ? '+' : '-', ' ')
        print(io, abs(round(c, digits=digits)))
        ex > 0 && print(io, sym)
        ex > 1 && print(io, '^', ex)
    end
    print(io, '\n')
end

function showpoly(coeffs::Vector{Vector{Float64}};sym::String="x",digits::Integer=2)
    print("1\n")
    map(c->showpoly(c,sym=sym,digits=digits),coeffs)
    print()
end

function showpoly(d::Integer,α::Vector{Float64},β::Vector{Float64};upto::Bool=true,sym::String="x",digits::Integer=2)
    @assert d >= 0 "degree has to be non-negative."
    d == 0 && return print("1")
    cfs = upto ? rec2coeff(d,α,β) : rec2coeff(d,α,β)[end]
    showpoly(cfs,sym=sym,digits=digits)
end

showpoly(d::Integer,op::OrthoPoly;upto::Bool=true,sym::String="x",digits::Integer=2) = showpoly(d,op.α,op.β;upto=upto,sym=sym,digits=digits)
showpoly(d::Integer,opq::OrthoPolyQ;upto::Bool=true,sym::String="x",digits::Integer=2) = showpoly(d,opq.op;upto=upto,sym=sym,digits=digits)

showpoly(α::Vector{Float64},β::Vector{Float64};sym::String="x",digits::Integer=2) = showpoly(length(α),α,β;sym=sym,digits=digits)
showpoly(op::Union{OrthoPoly,OrthoPolyQ};sym::String="x",digits::Integer=2) = showpoly(size(coeffs(op),1),op;sym=sym,digits=digits)

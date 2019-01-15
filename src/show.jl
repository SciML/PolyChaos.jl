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

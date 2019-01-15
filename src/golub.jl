function JacobiGW( n::Int64, a::Float64, b::Float64 )
    # Golub-Welsh for Gauss--Jacobi quadrature. This is used when max(a,b)>5.
    ab = a + b;
    ii = 2:n-1;
    abi = 2*ii + ab;
    aa = Float64[(b - a)/(2 + ab);
          (b^2 - a^2)./((abi - 2).*abi);
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))] ::Vector{Float64}
    bb = Float64[2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ;
          2 .*sqrt.(ii.*(ii .+ a).*(ii .+ b).*(ii .+ ab)./(abi.^2 .- 1))./abi] ::Vector{Float64}
    TT = SymTridiagonal(aa, bb)  # Jacobi matrix.
    x, V = eig( TT )                       # Eigenvalue decomposition.
    # Quadrature weights:
    #w = V[1,:].^2 .*( 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(2+ab) );
    w = V[1,:].^2*2^(ab+1)*beta(a+1, b+1);
    C = 2^(a+b+1)*beta(a+1, b+1)/sum(w)
    #@show C
    #w = C.*w;
    x, vec(w)
end

function golubwelsch(α,β)
    J = SymTridiagonal( α, sqrt.(β[2:end]) )
    nodes, V = eig(J)
    weights = V[1,:].^2
    return nodes, weights
end

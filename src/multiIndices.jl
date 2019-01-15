function calculateMultiIndices(d::Int, n::Int)
    #d denotes dimension of random variables/number of sources of uncertainty,
    #n the maximum degree of multivariate basis
  #function to calculate indices of multivariate basis following the algorithm
  #from "Spectral Methods for Uncertainty Quantiﬁcation; Le Maitre, Knio;2014" p.516-517
  # No::BigInt = factorial(BigInt(d+n))/(factorial(Bigint(d))*factorial(BigInt(n)));  #number of polynomials of multivariate basis
  No = numberPolynomials(d,n)
  if d==0
    print("Input zero in indicesMulti!")
    return zeros(Int64, 1,d)
  end

  # catch case n == 0 --> No-d==0
  if n==0
    inds=zeros(Int64,1,d)
    return inds
  end

  inds = vcat(zeros(Int64,1,d),Matrix(1I,d,d),zeros(Int64, No-d-1, d));  #initiate index matrix for basis
  pi=ones(Int64,No,d);

  for k=2:No
    g = 0;
    for l=1:d
      pi[k,l]=sum(pi[k-1,:])-g;
      g=g+pi[k-1,l];
    end
  end

  P=d+1;
  for k=2:n
    L = P;
    for j=1:d, m=L-pi[k,j]+1:L
        P=P+1;
        inds[P,:]=inds[m,:];
        inds[P,j]=inds[P,j]+1;
    end
  end

  return inds
end


function calculateMultiIndices(a::Array{T,1}) where T<:Integer
  #a contains maximum degrees of univariate bases
  #function to calculate indices of multivariate basis following the algorithm
  #from "Spectral Methods for Uncertainty Quantiﬁcation; Le Maitre, Knio;2014" p.516-517
  d=length(a)
  n=maximum(a)
  No::Int64 = numberPolynomials(d,n)  #number of polynomials of multivariate basis
  if d==0
    print("Input zero in indicesMulti!")
    return zeros(Int64, 1,d)
  end

  # catch case n == 0 --> No-d==0
  if n==0
    inds=zeros(Int64,1,d)
    return inds
  end

  inds = vcat(zeros(Int64,1,d),Matrix(1I,d,d),zeros(Int64, No-d-1, d));  #initiate index matrix for basis
  pi=ones(Int64,No,d);

  for k=2:No
    g = 0;
    for l=1:d
      pi[k,l]=sum(pi[k-1,:])-g;
      g=g+pi[k-1,l];
    end
  end

  P=d+1;
  for k=2:n
    L = P;
    for j=1:d, m=L-pi[k,j]+1:L
        P=P+1;
        inds[P,:]=inds[m,:];
        inds[P,j]=inds[P,j]+1;
    end
  end

  #remove elements of high order
  #TODO more efficient way to do that?
  for j=1:length(a)
    k=size(inds,1)
    l=1
    while l<=k
      if inds[l,j]>a[j]
        inds=vcat(inds[1:l-1,:],inds[l+1:end,:])
        k-=1
        l-=1
      end
      l+=1
    end
  end

  return inds
end


"""
    computes the number of polynomials with a multivariate basis
    `(d+n)!/(d!+n!)`
"""
function numberPolynomials(a::Array{T,1}) where T<:Integer
    d = length(a)
    n = maximum(a)

    x = max(d,n)
    y = min(d,n)

    return (prod(UInt128(x+1):UInt128(d+n))/factorial(UInt128(y)))
end

function numberPolynomials(d::Int64, n::Int64)
    x = max(d,n)
    y = min(d,n)
    return UInt128(prod(UInt128(x+1):UInt128(d+n))/factorial(UInt128(y)))
end

"""
    findUnivariateIndices(i::Int64,ind::Matrix{Int64})::Vector{Int64}
Given the multi-index `ind` this function returns all entries of the multivariate basis
that correspond to the `i`th univariate basis.
"""
function findUnivariateIndices(i::Int64,ind::Matrix{Int64})::Vector{Int64}
  l,p = size(ind)
  @assert i<=p "basis is $p-variate, you requested $i-variate"
  deg = ind[end,end]
  @assert deg>=0 "invalid degree"
  col = ind[:,i]
  myind = zeros(Int64,deg)
  for deg_ = 1:deg
      myind[deg_] = findfirst(x->x==deg_,col)
  end
  pushfirst!(myind,1)
end

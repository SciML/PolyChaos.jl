using PolyChaos, LinearAlgebra

function mytest(mOP::MultiOrthoPoly)
    @time P = [computeSP([i, j, k], mOP)
               for i in 0:(mOP.dim - 1), j in 0:(mOP.dim - 1), k in 0:(mOP.dim - 1)]
    @time T = Tensor(dim, mOP)
    P_ = [T.get([i, j, k]) for i in 0:(mOP.dim - 1), j in 0:(mOP.dim - 1), k in 0:(mOP.dim - 1)]
    @show norm(P - P_)
    return P_
end

##
deg, dim = 10, 3
α, β = 1.23, 3.23
## univariate
op = OrthoPoly("gaussian", deg; Nrec = 100)
opq = OrthoPolyQ(op)
mOP = MultiOrthoPoly([opq], deg)
mytest(mOP)
## multivariate
op2 = OrthoPoly("beta01", deg, Dict(:shape_a => α, :shape_b => β); Nrec = 100)
opq2 = OrthoPolyQ(op2)
mOP2 = MultiOrthoPoly([opq; opq2], deg)
mytest(mOP2)

mOP3 = MultiOrthoPoly([opq, opq, opq, opq], deg)
@time T = Tensor(dim, mOP3)

##
b = [He(deg)]
Gal = initBasis(b)
@time tensor(dim, Gal);

b2 = [He(deg), Ja(deg, α, β)]
Gal2 = initBasis(b2)
@time tensor(dim, Gal2);

b3 = [He(deg), He(deg), He(deg), He(deg), He(deg), He(deg)];
Gal3 = initBasis(b3);
@time tensor(dim, Gal3)

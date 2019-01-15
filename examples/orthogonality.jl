using PolyChaos, QuadGK


N = 10
gauss = OrthoPoly("gaussian",N)
pipj(t,i,j) = gauss.w(t)*evaluate(i,t,gauss)*evaluate(j,t,gauss)
[ quadgk(t->pipj(t,i,i),-Inf,Inf) - factorial(i) for i=1:N ]

op = OrthoPoly("uniform01",N)
pipj(t,i,j) = op.w(t)*evaluate(i,t,op)*evaluate(j,t,op)
for i = 0:N, j=0:N
    display(quadqk(t->op.w(t)*evaluate(i,t,op)*evaluate(j,t,op),0,1)[1])
end

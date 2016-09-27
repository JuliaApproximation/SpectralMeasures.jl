using Plots, ApproxFun, SpectralMeasures, Base.Test

############
### Tests
############

# Chebyshev U
a = []; b=[]
@time μ=spectralmeasure(a,b)
@test_approx_eq μ(-.99:.01:.99) sqrt(1-(-.99:.01:.99).^2)*2/π


# Chebyshev T
a = [0.]; b=[1/sqrt(2),.5]
@time μ=spectralmeasure(a,b)
@test_approx_eq μ(-.99:.01:.99) sqrt(1-(-.99:.01:.99).^2).^(-1)/π


# Legendre
n=100;a=zeros(n); b=(1:n-1)./sqrt(4*(1:n-1).^2-1)
μ=spectralmeasure(a,b)
@test_approx_eq_eps μ(-.99:.01:.99) 0.5ones(-.99:.01:.99) 0.001

# A strange hump
n=10;a=zeros(n); b=sqrt(1:(n-1))/sqrt(n)*0.5
@time μ=spectralmeasure(a,b)
@test_approx_eq μ(0.1) 2.031454229315879 # empirical


# Jacobi polynomials
α = 1.2; β = 1.1
n=107;a = (β.^2-α.^2)./((2.*(0:n)+α+β).*(2.*(1:n+1)+α+β))
b = 2*sqrt(((1:n).*(α+(1:n)).*(β+(1:n)).*(α+β+(1:n)))./((2.*(1:n)+α+β-1).*((2.*(1:n)+α+β).^2).*(2.*(1:n)+α+β+1)))
μ=spectralmeasure(a,b)
plot(μ)
x=Fun()
ν=(1-x)^α*(x+1)^β
ν/=sum(ν)

@test_approx_eq_eps μ(-.99:.01:.99) ν(-.99:.01:.99) 0.0001

# Simplest perturbation
k = 17
  a = [k/20];b=[]
  μ=spectralmeasure(a,b)
  ApproxFun.plot(μ)
k = -19
a = [k/20];b=[]
ν=spectralmeasure(a,b)
ApproxFun.plot(ν)


# Can you make this work, Sheehan? MW
ApproxFun.plot(μ-ν)


L=DiscreteLaplacian()
K=SymTriOperator(-ones(5),zeros(5))

J=L+0.5K
x,Q=eig(J)

@which eig(J)

@which ql(J-3.0*I)
A=J-3.0*I
@which ql(A.dv,A.ev,A.a,A.b)
Q,L=ql(J-3.0*I)

norm(L[1:100,1:100]  -[L[k,j] for k=1:100,j=1:100])

@which ql(J-3.0*I)

full(Q[1:100,1:100])*full(L[1:100,1:100])  |>chopm
[(v=Q*[zeros(k-1);1.0] ;
    norm(v- Q[1:length(v),k])) for k=1:100]  |>norm
using SO

show(Q*L)

n=max(size(L.K.matrix,1),length(Q.s)+3)

bandinds(L)==(-2,0)

#if bandinds(L)==(-2,0)
     # We check if L*Q is tridiagin al
tol=1E-14*(maximum(L.T)+maximum(L.K))
istri=true

show(L)

k=n
abs(L[k,k-2]*hc(Q,k-1)+L[k,k-1]*hs(Q,k-2)*hc(Q,k)+L[k,k]*hs(Q,k-2)*hs(Q,k-1)*hc(Q,k+1))

for k=3:n
    @show k
    if abs(L[k,k-2]*hc(Q,k-1)+L[k,k-1]*hs(Q,k-2)*hc(Q,k)+L[k,k]*hs(Q,k-2)*hs(Q,k-1)*hc(Q,k+1))>tol
        istri=false
        break
    end
end

istri
if istri
    issym=true
    if !isapprox(-L[1,1]*hs(Q,1),L[2,1]*hc(Q,1)*hc(Q,2)+L[2,2]*hc(Q,1)*hc(Q,3)*hs(Q,1);atol=tol)
        issym=false
    end

    if issym
        for k=2:n+1  # kth row
            if !isapprox(-L[k+1,k-1]*hs(Q,k-1)+L[k+1,k]*hc(Q,k)*hc(Q,k+1)+L[k+1,k+1]*hc(Q,k)*hc(Q,k+2)*hs(Q,k),
            -L[k,k]*hs(Q,k);atol=tol)
                issym=false
                break
            end
        end
    end

    if issym
       # result is SymTriToeplitxz

        ev=Array(Float64,max(min(size(L.K.matrix,1),size(L.K.matrix,2)),
                             length(Q.s)))
        for k=1:length(ev)
            ev[k]=-L[k,k]*hs(Q,k)
        end

        dv=Array(Float64,max(length(Q.s)+1,size(L.K.matrix,1)))
        dv[1]=hc(Q,1)*hc(Q,2)*L[1,1]
        for k=2:length(dv)
            dv[k]=-hs(Q,k-1)*L[k,k-1]+hc(Q,k)*hc(Q,k+1)*L[k,k]
        end

        t1=-L.T[1,1]*Q.s∞
        t0=-Q.s∞*L.T[2,1]+Q.c∞^2*L.T[1,1]

        si=Q.sign?1:-1
        return SymTriToeplitz(si*dv,si*ev,si*t0,si*t1)
    end
end
#end

# default constructor
TimesOperator(L,Q)

bandinds(L)
# check some ApproxFun bugs
n=100
QM=full(Q[1:n,1:n])
@test_approx_eq pad((Q*[1.] ).coefficients,n) QM[:,1]


@test_approx_eq Q.op.ops[end][1:n,1:n]\[1.;zeros(n-1)] pad(Q.op.ops[end]\[1.],n)

b=Q.op.ops[end]\[1.]

@test_approx_eq Q.op.ops[1][1:n,1:n]\pad(b,n) pad((Q.op.ops[1]\b).coefficients,n)


@time u=Q\(exp(im*x)*(Q*[1.]))
@test_approx_eq u expm(im*full(J[1:200,1:200]))[1:length(u)]

t=10.0
    @time u=Q\(exp(im*t*x)*(Q*[1.]))  # 0.04s
    scatter([real(u) imag(u)])


t=100000.0
    @time u=Q\(exp(im*t*x)*(Q*[1.]))  # 0.04s
    scatter([real(u) imag(u)])


n=100
    J[1:n,1:n]
n=10000
    @time eig(SymTridiagonal(J,1:n,1:n))


@time exp(im*t*x)
Q*[.1]


typeof(Q)



Q.op.ops[1].op
Q.op.ops[1].mat11

Q.op.ops[1]
@time Q\[1.,2.,3.,4.,5.,6.]
@time v=Q.op.ops[1]\[1.,2.,3.,4.,5.,6.]

Qt.ops
Qt=(Q.op.ops[2]')
Qt.ops[1]*(Qt.ops[2]*(Qt.ops[3]*v.coefficients))
@which Qt*v
Profile.print()
v=exp(im*t*x)*(Q*[1.])

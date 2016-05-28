using Plots, ApproxFun, SpectralMeasures, Base.Test

############
### Tests
############

# Chebyshev U
a = []; b=[]
  @time μ=spectralmeasureU(a,b)
  @time ν=spectralmeasureT(a,b)
  @test_approx_eq μ(-.99:.01:.99) ν(-.99:.01:.99)

# Chebyshev T
a = [0.]; b=[1/sqrt(2),.5]
  @time μ=spectralmeasureU(a,b)
  @time ν=spectralmeasureT(a,b)
  @test_approx_eq μ(-.99:.01:.99) ν(-.99:.01:.99)


# Legendre
n=100;a=zeros(n); b=(1:n-1)./sqrt(4*(1:n-1).^2-1)
  μ=spectralmeasureU(a,b)
  ν=spectralmeasureT(a,b)
  @test_approx_eq μ(-.99:.01:.99) ν(-.99:.01:.99)

# A strange hump
n=10;a=zeros(n); b=sqrt(1:(n-1))/sqrt(n)*0.5
  @time μ=spectralmeasureU(a,b)
  @time ν=spectralmeasureT(a,b)
  maximum(μ(-.99:.01:.99)-ν(-.99:.01:.99))


# Jacobi polynomials
α = 1.2; β = 1.1
  n=107;a = (β.^2-α.^2)./((2.*(0:n)+α+β).*(2.*(1:n+1)+α+β))
  b = 2*sqrt(((1:n).*(α+(1:n)).*(β+(1:n)).*(α+β+(1:n)))./((2.*(1:n)+α+β-1).*((2.*(1:n)+α+β).^2).*(2.*(1:n)+α+β+1)))
  μ=spectralmeasureU(a,b)
  ν=spectralmeasureT(a,b)
  @test_approx_eq μ(-.99:.01:.99) ν(-.99:.01:.99)

# Simplest perturbation
k = 17
  a = [k/20];b=[]
  μ=spectralmeasureU(a,b)
  ApproxFun.plot(μ)
k = -19
  a = [k/20];b=[]
  ν=spectralmeasureT(a,b)
  ApproxFun.plot(ν)

ApproxFun.plot(μ)
ApproxFun.plot(ν)

sum(μ)
sum(ν)

length(μ)
length(ν)

# Can you make this work, Sheehan? MW
ApproxFun.plot(μ-ν)


L=DiscreteLaplacian()
K=SymTriOperator(-ones(5),zeros(5))

J=L+0.5K
x,Q=eig(J)

@time u=Q\(exp(im*x)*(Q*[1.]))
@test_approx_eq u expm(im*full(J[1:100,1:100]))[1:length(u)]

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

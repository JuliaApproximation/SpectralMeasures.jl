using Plots, ApproxFun, SpectralMeasures, Base.Test
    import ApproxFun: A_ldiv_B_coefficients
############
### Tests
############

# Chebyshev U
a = Float64[]; b= Float64[]
@time μ=spectralmeasure(a,b)
@test μ.(-.99:.01:.99) ≈ sqrt.(1.-(-.99:.01:.99).^2)*2/π


# Chebyshev T
a = [0.]; b=[1/sqrt(2),.5]
@time μ=spectralmeasure(a,b)
@test μ.(-.99:.01:.99) ≈ sqrt.(1.-(-.99:.01:.99).^2).^(-1)/π


# Legendre
n=100;a=zeros(n); b=(1:n-1)./sqrt.(4*(1:n-1).^2-1)
μ=spectralmeasure(a,b)
@test μ.(-.99:.01:.99) ≈ 0.5ones(-.99:.01:.99) atol=0.001

# A strange hump
n=10;a=zeros(n); b=sqrt.(1:(n-1))./(2sqrt(n))
@time μ=spectralmeasure(a,b)
@test μ(0.1) ≈ 2.031454229315879 atol=1E-9 # empirical


# Jacobi polynomials
α = 1.2; β = 1.1
n=107;a = (β.^2-α.^2)./((2.*(0:n)+α+β).*(2.*(1:n+1)+α+β))
b = 2*sqrt.(((1:n).*(α+(1:n)).*(β+(1:n)).*(α+β+(1:n)))./((2.*(1:n)+α+β-1).*((2.*(1:n)+α+β).^2).*(2.*(1:n)+α+β+1)))
μ=spectralmeasure(a,b)
# plot(μ)
x=Fun()
ν=(1-x)^α*(x+1)^β
ν/=sum(ν)

@test μ.(-.99:.01:.99) ≈ ν.(-.99:.01:.99) atol=0.001


# Simplest perturbation
k = 17
  a = [k/20];b=Float64[]
  μ=spectralmeasure(a,b)
  # plot(μ)


k = -19
a = [k/20];b=Float64[]
ν=spectralmeasure(a,b)
# plot(ν)


# Can you make this work, Sheehan? MW
#ApproxFun.plot(μ-ν)

# check robustness to extra zerods
for n=3:6,m=0:2
    L=freejacobioperator()
    K=SymTriOperator(-[ones(5);zeros(m)],zeros(n))

    J=L+0.5K
    @test isa(J,SymTriToeplitz)
    @test full(J[1:3,1:6]) == [-0.5  0.5  0   0 0 0;
                                0.5 -0.5  0.5 0 0 0;
                                0    0.5 -0.5 0.5 0 0]

    A=J-3.0*I
    Q,L=ql(A)

    @test norm((Q*L-A)[1:30,1:30]) < 100eps()

    @test isa(L*Q,SpectralMeasures.SymTriToeplitz)
end



L=freejacobioperator()
K=SymTriOperator(-ones(5),zeros(4))

J=L+0.5K

x,Q=eig(J)

r = A_ldiv_B_coefficients(Q,A_mul_B_coefficients(Q,[1.0]))
@test r ≈ [1;zeros(length(r)-1)]
@test coefficients(Q\(Q*[1.0])) ≈ [1;zeros(length(r)-1)]

@test Q\(Q*[1.0]) ≈ Fun([1;zeros(length(r)-1)],domainspace(Q))
@test Q\(x*(Q*[1.0])) ≈ J*[1.0]

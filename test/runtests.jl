using ApproxFun, SpectralMeasures, Test, LinearAlgebra
    import ApproxFun: ldiv_coefficients
############
### Tests
############

@testset "Spectral measures" begin
    @testset "Chebyshev U" begin
        a = Float64[]; b= Float64[]
        @time μ=spectralmeasure(a,b)
        @test μ.(-.99:.01:.99) ≈ sqrt.(1 .-(-.99:.01:.99).^2)*2/π
    end

    @testset "Chebyshev T" begin
        a = [0.]; b=[1/sqrt(2),.5]
        @time μ=spectralmeasure(a,b)
        @test μ.(-.99:.01:.99) ≈ sqrt.(1 .-(-.99:.01:.99).^2).^(-1)/π
    end


    @testset "Legendre" begin
        n=100; a=zeros(n); b=(1:n-1)./sqrt.(4*(1:n-1).^2 .-1)
        μ=spectralmeasure(a,b)
        @test μ.(-.99:.01:.99) ≈ 0.5one.(-.99:.01:.99) atol=0.001
    end

    @testset "A strange hump" begin
        n=10;a=zeros(n); b=sqrt.(1:(n-1))./(2sqrt(n))
        @time μ=spectralmeasure(a,b)
        @test μ(0.1) ≈ 2.031454229315879 atol=1E-9 # empirical
    end

    @testset "Jacobi polynomials" begin
        α = 1.2; β = 1.1
        n=107;a = @. (β^2-α^2)/((2(0:n)+α+β)*(2*(1:n+1)+ (α+β)))
        b = @. 2sqrt(((1:n)*(α+(1:n))*(β+(1:n))*(α+β+(1:n)))/
                ((2(1:n)+α+β-1)*((2(1:n)+α+β)^2)*(2*(1:n) + α+β+1)))
        μ=spectralmeasure(a,b)
        x=Fun()
        ν=(1-x)^α*(x+1)^β
        ν/=sum(ν)

        @test μ.(-.99:.01:.99) ≈ ν.(-.99:.01:.99) atol=0.001
    end

    @testset "Simplest perturbation" begin
        k = 17
        a = [k/20];b=Float64[]
        μ=spectralmeasure(a,b)


        k = -19
        a = [k/20];b=Float64[]
        ν=spectralmeasure(a,b)
    end

    @testset "Check robustness for extra zeros" begin
        for n=3:6,m=0:2
            L=freejacobioperator()
            K=SymTriOperator(-[ones(5);zeros(m)],zeros(n))

            J=L+0.5K
            @test isa(J,SymTriPertToeplitz)
            @test Matrix(J[1:3,1:6]) == [-0.5  0.5  0   0 0 0;
                                        0.5 -0.5  0.5 0 0 0;
                                        0    0.5 -0.5 0.5 0 0]

            A=J-3.0*I
            Q,L=ql(A)

            @test norm((Q*L-A)[1:30,1:30]) < 100eps()

            @test isa(L*Q,SpectralMeasures.SymTriPertToeplitz)
        end
    end
end


@testset "eig ops" begin
    L=freejacobioperator()
    K=SymTriOperator(-ones(5),zeros(4))

    J=L+0.5K

    D,Q=eigen(J)

    r = ldiv_coefficients(Q, mul_coefficients(Q,[1.0]))
    @test r ≈ [1;zeros(length(r)-1)]
    @test coefficients(Q\(Q*[1.0])) ≈ [1;zeros(length(r)-1)]

    @test Q\(Q*[1.0]) ≈ Fun([1;zeros(length(r)-1)],domainspace(Q))
    @test Q\(D*(Q*[1.0])) ≈ J*[1.0]
end

@testset "zero b" begin
    J=SymTriPertToeplitz(ones(0),Float64[],0.0,0.5) + 1.01I
    D,Q = eigen(J)

    @test Q\(D*(Q*[1.0])) ≈ J*[1.0]

    J = -freejacobioperator()
    D,Q = eigen(J)

    @test Q\(D*(Q*[1.0])) ≈ J*[1.0]


    J=SymTriPertToeplitz(ones(0),Float64[],0.0,-0.5) + 1.01I

    Q,L = ql(J)
    @test norm((Q*L-J)[1:100,1:100]) ≤ 100eps()


    @time D,Q = eigen(J)
    @test Q\(D*(Q*[1.0])) ≈ J*[1.0]
end


@testset "Zeros in QL (#19)" begin
    J = SymTriPertToeplitz([2.0], Float64[0.0000000001], -2, 0.5)
    Q,L = ql(J)
    @test (Q*L)[1:10,1:10] ≈ J[1:10,1:10]

    J = SymTriPertToeplitz([2.0], Float64[0.0], -2, 0.5)
    Q,L = ql(J)
    @test (Q*L)[1:10,1:10] ≈ J[1:10,1:10]
end

@testset "Diagonal operator functions" begin
    Δ = freejacobioperator()
    K=SymTriOperator(-2ones(3),zeros(4))
    J = Δ + K
    D,U = eigen(J)

    v = randn(10)
    @test pad((U\(exp(D)*U*v)).coefficients,100) ≈ exp(Matrix(J[1:100,1:100]))*pad(v,100)
end

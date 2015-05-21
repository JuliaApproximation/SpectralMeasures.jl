############
### Tests
############

sum(μ)

μ.coefficients[1:4]

n=20;a=zeros(n);b=sqrt(1:(n-1))/sqrt(n)*0.5
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)


n=5;a=1./([1:n;]);b=sqrt(1:(n-1))/sqrt(n)*0.5
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)



n=51;a=zeros(n);b=sqrt(1:(n-1))/sqrt(n)*0.5
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)


n=11;a=0.2rand(n);b=rand(n-1)
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)




n=14;a=(rand(n)-0.5);b=fill(0.5,n-1)
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)



n=20;a=0.5cos([1:n]);b=fill(0.5,n-1)
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)


α=44;β=44;a=[α/50+0.02145238,0.];b=[β/20+0.012412];
    μ=spectralmeasure(a,b)
    ApproxFun.plot(μ)

values(μ)

ApproxFun.plot(μ)
μ
T,K = tkoperators(a,b)
T[1:10,1:10]
K[1:10,1:10]

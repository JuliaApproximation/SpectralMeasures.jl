using ApproxFun
using SpectralMeasure
using Gadfly
set_default_plot_format(:svg)

############
### Tests
############

# Chebyshev U
a = []; b=[]
  μ=spectralmeasureU(a,b)
  ν=spectralmeasureT(a,b)
  maximum(evaluate(μ,[-.99:.01:.99])-evaluate(ν,[-.99:.01:.99]))

# Chebyshev T
a = [0.]; b=[1/sqrt(2),.5]
  μ=spectralmeasureU(a,b,maxlength=1000000)
  ν=spectralmeasureT(a,b)
  maximum(evaluate(μ,[-.99:.01:.99])-evaluate(ν,[-.99:.01:.99]))

# Legendre
n=100;a=zeros(n); b=[1:n-1]./sqrt(4*[1:n-1].^2-1)
  μ=spectralmeasureU(a,b)
  ν=spectralmeasureT(a,b)
  maximum(evaluate(μ,[-.99:.01:.99])-evaluate(ν,[-.99:.01:.99]))

# A strange hump
n=35;a=zeros(n); b=sqrt(1:(n-1))/sqrt(n)*0.5
  μ=spectralmeasureU(a,b)
  ν=spectralmeasureT(a,b)
  maximum(evaluate(μ,[-.99:.01:.99])-evaluate(ν,[-.99:.01:.99]))

# Jacobi polynomials
α = 1.2; β = 1.1
  n=107;a = (β.^2-α.^2)./((2.*[0:n]+α+β).*(2.*[1:n+1]+α+β))
  b = 2*sqrt(([1:n].*(α+[1:n]).*(β+[1:n]).*(α+β+[1:n]))./((2.*[1:n]+α+β-1).*((2.*[1:n]+α+β).^2).*(2.*[1:n]+α+β+1)))
  μ=spectralmeasureU(a,b)
  ν=spectralmeasureT(a,b)
  maximum(evaluate(μ,[-.99:.01:.99])-evaluate(ν,[-.99:.01:.99]))

ApproxFun.plot(μ)
ApproxFun.plot(ν)

sum(μ)
sum(ν)

length(μ)
length(ν)

# Can you make this work, Sheehan? MW
ApproxFun.plot(μ-ν)


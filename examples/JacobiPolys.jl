using ApproxFun
using SpectralMeasures
using Plots; Plots.gr(legend=false)

Δ = DiscreteLaplacian()

# Legendre polynomials
n=100;a=zeros(n); b=(1:n-1)./sqrt(4*(1:n-1).^2-1)
  μ

# Ultraspherical polynomials

# Jacobi polynomials
α = 1.2; β = 1.1
  n=107;a = (β.^2-α.^2)./((2.*(0:n)+α+β).*(2.*(1:n+1)+α+β))
  b = 2*sqrt(((1:n).*(α+(1:n)).*(β+(1:n)).*(α+β+(1:n)))./((2.*(1:n)+α+β-1).*((2.*(1:n)+α+β).^2).*(2.*(1:n)+α+β+1)))
  μ

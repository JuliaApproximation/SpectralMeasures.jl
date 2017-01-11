module SpectralMeasures
using Base, Compat, ApproxFun, Plots, ComplexPhasePortrait

import Base:+,-,*,/,.*,.-,./,.+,getindex

import ApproxFun: Operator, ToeplitzOperator, DiracSpace, plot, IdentityOperator,
            TridiagonalOperator, setdomain, resizedata!, bandinds, PointSpace,
            BandedMatrix, bzeros, TimesOperator, BlockOperator, SpaceOperator, AbstractCount, UnitCount,
            MatrixSpace, ∞, ℓ⁰, domainspace, rangespace, domain

export spectralmeasure, discreteEigs, principalResolvent, discResolvent
export connectionCoeffsOperator, applyConversion, SymTriOperator, SymTriToeplitz
export triplePlot
export FreeJacobiOperator, jacobioperator, ql

include("HessenbergUnitary.jl")
include("PertToeplitz.jl")
include("helper.jl")
include("ql.jl")
include("RatFun.jl")

function spectralmeasure(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  c = Fun(Taylor,C.T.nonnegative)
  f = Fun(C*(C'*[1]),Ultraspherical(1))

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1) ,complexroots(c))))
  if length(z) > 0
     cprime = differentiate(c)
     eigs=real(map(joukowsky,z))
     weights = (z-1./z).^2./(z.*real(cprime(z)).*real(c(1./z)))
     p = Fun(DiracSpace(eigs),weights) + Fun(JacobiWeight(.5,.5,Ultraspherical(1)),[2/pi])
     q = Fun(PointSpace(eigs),ones(length(eigs))) + f
     μ = RatFun(p,q)
  else
    μ = RatFun(Fun(JacobiWeight(.5,.5,Ultraspherical(1)),[2/pi]),f)
  end
  μ
end

function principalResolvent(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Compute the necessary polynomials
  C = connectionCoeffsOperator(a,b)
  Cmu = connectionCoeffsOperator(a[2:end],b[2:end]) # Technically not Cmu from the paper
  f = Fun((C*(C'*[1])),Ultraspherical(1))
  fmu = Fun(Ultraspherical(1),Cmu*((C'*[1]).coefficients[2:end])/b[1])

  # Return the resolvent
  x->(2*sqrt(complex(x-1)).*sqrt(complex(x+1))-2*x-extrapolate(fmu,x))./extrapolate(f,x)
end

function discResolvent(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Compute the necessary polynomials
  C = SpectralMeasures.connectionCoeffsOperator(a,b)
  Cmu = SpectralMeasures.connectionCoeffsOperator(a[2:end],b[2:end]) # Technically not Cmu from the paper
  c = Fun(Taylor,C.T.nonnegative)
  cmu = Fun(Taylor,[0;Cmu.T.nonnegative]/b[1])

  # Return the rational function
  x->-cmu(x)./c(x)
end

function discreteEigs(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that C*J = Toeplitz([0,1/2])*C
  C = connectionCoeffsOperator(a,b)
  Tfun = Fun(Taylor,C.T.nonnegative)
  sort(real(map(joukowsky,filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1),complexroots(Tfun)))))
end

#Finds C such that C'(U_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and U_k is Chebyshev U
function connectionCoeffsOperator(a,b)
  n = max(2,length(a),length(b)+1)
  N = 2*n #This is sufficient only because we go from Chebyshev U
  a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
  ToeplitzVec = zeros(N)
  K = bzeros(Float64,n,N,0,N+1)
  K[1,1] = 1
  K[1,2] = -a[1]/b[1]
  K[2,2] = .5/b[1]
  # The recurrence for the first n+1 cols depend on a and b
  for j = 3:n+1
    K[1,j] = (-a[j-1]*K[1,j-1] + .5*K[2,j-1] - b[j-2]*K[1,j-2])/b[j-1]
    for i = 2:j-2
      K[i,j] = (.5*K[i-1,j-1] -a[j-1]*K[i,j-1] + .5*K[i+1,j-1] - b[j-2]*K[i,j-2])/b[j-1]
    end
    K[j-1,j] = (.5*K[j-2,j-1] -a[j-1]*K[j-1,j-1] - b[j-2]*K[j-1,j-2])/b[j-1]
    if j<n+1
      K[j,j] = .5*K[j-1,j-1]/b[j-1]
    end
  end
  ToeplitzVec[1] = K[n,n]
  ToeplitzVec[2] = K[n,n+1]
  # The recurrence for rows n+2 to 2n do not depend on a and b
  for j = n+2:N
    K[1,j] = K[2,j-1] - K[1,j-2]
    for i = 2:N-j
      K[i,j] = K[i-1,j-1] + K[i+1,j-1] - K[i,j-2]
    end
    if j < N
      K[N+1-j,j] = K[N-j,j-1] + K[N+2-j,j-1] - K[N+1-j,j-2]
    end
    ToeplitzVec[2*(j-n)-1] = K[N+1-j,j-1]
    ToeplitzVec[2*(j-n)] = K[N+1-j,j]
  end
  T = ToeplitzOperator(Float64[],chop!(ToeplitzVec))
  for j = 1:N
    for i = 1:min(j,N+1-j)
      K[i,j]-=T[i,j]
    end
  end
  T+FiniteOperator(K)
end

function triplePlot(a,b,Z=linspace(-3, 3, 300).+linspace(3,-3,300)'*im)
  # Build the measure and resolvents
  μ = spectralmeasure(a,b)
  R = principalResolvent(a,b)
  r = discResolvent(a,b)

  # Create the subplots
  p1 = plot(μ,xlims=(-2,2),ylims=(0,1.5))
  p2 = plot(portrait(R(Z),PTstepmod),xlims=(-3,3),ylims=(-3,3),aspect_ratio=1)
  plot!([-1.,1.],[0.,0.],linewidth=3,color=:black)
  p3 = plot(portrait(r(Z),PTstepmod),xlims=(-3,3),ylims=(-3,3),aspect_ratio=1)
  plot!(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=3,color=:black)

  # Plot the subplots
  l = @layout [a{0.7w}; b c]
  plot(p1,p2,p3,layout=l)
end

end  #Module

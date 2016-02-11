module SpectralMeasures
using Base, Compat, ApproxFun

import Base: +,-,*,/

import ApproxFun:BandedOperator,ToeplitzOperator,DiracSpace, plot, IdentityOperator,
TridiagonalOperator,addentries!,setdomain, SavedBandedOperator, resizedata!

export spectralmeasure, spectralmeasureRat, spectralmeasureT, discreteEigs, connectionCoeffsOperator

export DiscreteLaplacian

include("helper.jl")
include("ql.jl")
include("PertToeplitz.jl")
include("plot.jl")

spectralmeasure(a,b) = spectralmeasureRat(a,b)

function spectralmeasureRat(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  c = Fun([C.T[1,1];C.T.negative],Taylor)
  f = Fun(C'*(C*[1]),Ultraspherical{1}())

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && isreal(z),complexroots(c))))
  if length(z) > 0
    cprime = differentiate(c)
    eigs=real(map(joukowsky,z))
    weights = (z-1./z).^2./(z.*real(cprime(z)).*real(c(1./z)))
    μ = RatFun(Fun(2/pi,JacobiWeight(.5,.5,Ultraspherical{1}())),f)
    μ1 + Fun(weights,DiracSpace(eigs))
  else
    μ = RatFun(Fun(2/pi,JacobiWeight(.5,.5,Ultraspherical{1}())),f)
  end
  μ
end

function spectralmeasureT(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  c = Fun([C.T[1,1];C.T.negative],Taylor)
  f = Fun(C'*(C*[1]),Ultraspherical{1}())

  # Compute continuous part of measure
  coeffs = (2/pi)*((1-x^2)./f).coefficients
  coeffs = [coeffs[1];sqrt(2)*coeffs[2:end]] #normalised T_k polynomials
  μ = Fun(coeffs,JacobiWeight(-.5,-.5,Ultraspherical{0}()))

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && isreal(z),complexroots(c))))
  if length(z) > 0
    cprime = differentiate(c)
    eigs=real(map(joukowsky,z))
    weights = (z-1./z).^2./(z.*real(cprime(z)).*real(c(1./z)))
    μ + Fun(weights,DiracSpace(eigs))
  else
    μ
  end
end

spectralmeasure(a...;opts...)=spectralmeasureT(a...;opts...) # default to T

function discreteEigs(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  Tfun = Fun([C.T[1,1];C.T.negative],Taylor)
  sort(real(map(joukowsky,filter!(z->abs(z)<1 && isreal(z),complexroots(Tfun)))))
end

#Finds C such that C(U_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and U_k is Chebyshev U
function connectionCoeffsOperator(a,b)
  n = max(length(a),length(b)+1)
  N = 2*n+1 #This is sufficient only because we go from Chebyshev U
  a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
  ToeCol = zeros(N)
  K = zeros(N,n+1)
  K[1,1] = 1
  K[2,1] = -a[1]/b[1]
  K[2,2] = .5/b[1]
  for i = 3:n+1
    K[i,1] = (-a[i-1]*K[i-1,1] + .5*K[i-1,2] - b[i-2]*K[i-2,1])/b[i-1]
    for j = 2:i-1
      K[i,j] = (.5*K[i-1,j-1] -a[i-1]*K[i-1,j] + .5*K[i-1,j+1] - b[i-2]*K[i-2,j])/b[i-1]
    end
    K[i,i] = .5*K[i-1,i-1]/b[i-1]
  end
  for i = n+2:N
    K[i,1] = K[i-1,2] - K[i-2,1]
    for j = 2:N+1-i
      K[i,j] = K[i-1,j-1] + K[i-1,j+1] - K[i-2,j]
    end
    if i == n+2
      K[n+2,n+1] = K[n+1,n] - K[n,n+1]
    else
      K[i,N+2-i] = K[i-1,N+1-i] + K[i-1,N+3-i] - K[i-2,N+2-i]
    end
    k = 2*(i-n-1)-1
    ToeCol[k] = K[i-1,N+2-i]
    ToeCol[k+1] = K[i,N+2-i]
  end
  T = ToeplitzOperator(chop!(ToeCol[2:N]),[ToeCol[1]])
  for i = 1:N
    for j = 1:min(i,N+2-i)
      K[i,j]-=T[i,j]
    end
  end
  T+FiniteOperator(K)
end

#Finds NxN truncation of C such that C(Q_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and Q_k has Jacobi coeffs c,d
function connectionCoeffsMatrix(a,b,c,d,N)
  @assert N >= 2*max(length(a),length(b),length(c),length(d))
  a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
  c = [c;zeros(N-length(c))]; d = [d;.5+zeros(N-length(d))]

  C = zeros(N,N)
  C[1,1] = 1
  C[2,1] = (c[1]-a[1])/b[1]
  C[2,2] = d[1]/b[1]
  for i = 3:N
    C[i,1] = ((c[1]-a[i-1])*C[i-1,1] + d[1]*C[i-1,2] - b[i-2]*C[i-2,1])/b[i-1]
    for j = 2:i-1
      C[i,j] = (d[j-1]*C[i-1,j-1] + (c[j]-a[i-1])*C[i-1,j] + d[j]*C[i-1,j+1] - b[i-2]*C[i-2,j])/b[i-1]
    end
    C[i,i] = d[i-1]*C[i-1,i-1]/b[i-1]
  end
  C
end

# This is for Chebyshev U
connectionCoeffsMatrix(a,b,N) = connectionCoeffsMatrix(a,b,[],[],N)

end  #Module

module SpectralMeasures
using Base, Compat, ApproxFun, Plots

import Base: +,-,*,/

import ApproxFun:BandedOperator,ToeplitzOperator,DiracSpace, plot, IdentityOperator,
            TridiagonalOperator,addentries!,setdomain, SavedBandedOperator, resizedata!, bandinds, PointSpace

export spectralmeasure, spectralmeasureRat, spectralmeasureU, spectralmeasureT, discreteEigs, connectionCoeffsOperator

export DiscreteLaplacian, jacobioperator,ql, RatFun

include("helper.jl")
include("ql.jl")
include("PertToeplitz.jl")
include("plot.jl")
include("RatFun.jl")

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
    #error("Can't deal with discrete spectrum until PointsSpace is fully implemented.")
     cprime = differentiate(c)
     eigs=real(map(joukowsky,z))
     weights = (z-1./z).^2./(z.*real(cprime(z)).*real(c(1./z)))
     p = Fun([2/pi],JacobiWeight(.5,.5,Ultraspherical{1}())) + Fun(weights,DiracSpace(eigs))
     q = f + Fun(ones(length(eigs)),PointSpace(eigs))
     μ = RatFun(p,q)
  else
    μ = RatFun(Fun([2/pi],JacobiWeight(.5,.5,Ultraspherical{1}())),f)
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

  # Compute continuous part of measure
  coeffs = Fun(x->(2/pi)*(1-x.^2)./abs(c(x+im*sqrt(1-x.^2))).^2,Ultraspherical{0}()).coefficients
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

function spectralmeasureU(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  c = Fun([C.T[1,1];C.T.negative],Taylor)
  f = Fun(C'*(C*[1]),Ultraspherical{1}())

  # Compute continuous part of measure
  coeffs = (1./f).coefficients
  μ = Fun((2/pi)*coeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))

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
  N = 2*n #This is sufficient only because we go from Chebyshev U
  a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
  ToeplitzVec = zeros(N)
  K = zeros(N,n)
  K[1,1] = 1
  K[2,1] = -a[1]/b[1]
  K[2,2] = .5/b[1]
  # The recurrence for the top n+1 rows depend on a and b
  for i = 3:n+1
    K[i,1] = (-a[i-1]*K[i-1,1] + .5*K[i-1,2] - b[i-2]*K[i-2,1])/b[i-1]
    for j = 2:i-2
      K[i,j] = (.5*K[i-1,j-1] -a[i-1]*K[i-1,j] + .5*K[i-1,j+1] - b[i-2]*K[i-2,j])/b[i-1]
    end
    K[i,i-1] = (.5*K[i-1,i-2] -a[i-1]*K[i-1,i-1] - b[i-2]*K[i-2,i-1])/b[i-1]
    if i<n+1
      K[i,i] = .5*K[i-1,i-1]/b[i-1]
    end
  end
  ToeplitzVec[1] = K[n,n]
  ToeplitzVec[2] = K[n+1,n]
  # The recurrence for rows n+2 to 2n do not depend on a and b
  for i = n+2:N
    K[i,1] = K[i-1,2] - K[i-2,1]
    for j = 2:N-i
      K[i,j] = K[i-1,j-1] + K[i-1,j+1] - K[i-2,j]
    end
    if i < N
      K[i,N+1-i] = K[i-1,N-i] + K[i-1,N+2-i] - K[i-2,N+1-i]
    end
    ToeplitzVec[2*(i-n)-1] = K[i-1,N+1-i]
    ToeplitzVec[2*(i-n)] = K[i,N+1-i]
  end
  T = ToeplitzOperator(chop!(ToeplitzVec[2:N]),[ToeplitzVec[1]])
  for i = 1:N
    for j = 1:min(i,N+1-i)
      K[i,j]-=T[i,j]
    end
  end
  T+FiniteOperator(K)
end

#Finds NxN truncation of C such that C(Q_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and Q_k has Jacobi coeffs c,d
function connectionCoeffsMatrix(a,b,c,d,N)
  if N>max(length(a),length(b)+1,length(c),length(d)+1)
    a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
    c = [c;zeros(N-length(c))]; d = [d;.5+zeros(N-length(d))]
  end

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

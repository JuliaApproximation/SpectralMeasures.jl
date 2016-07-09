module SpectralMeasures
using Base, Compat, ApproxFun, Plots

import Base:+,-,*,/,.*,.-,./,.+,getindex

import ApproxFun:BandedOperator, ToeplitzOperator, DiracSpace, plot, IdentityOperator,
            TridiagonalOperator,addentries!,setdomain, SavedBandedOperator, resizedata!, bandinds, PointSpace,
            BandedMatrix, bzeros, TimesOperator, BlockOperator, SpaceOperator, AnySpace, AbstractCount, UnitCount,
            SubBandedMatrix, linsolve, MatrixSpace, ∞

export spectralmeasure, spectralmeasureRat, spectralmeasureU, spectralmeasureT, discreteEigs, connectionCoeffsOperator, applyConversion,
        SymTriOperator, SymTriToeplitz

export DiscreteLaplacian, jacobioperator, ql

include("helper.jl")
include("HessenbergUnitary.jl")
include("PertToeplitz.jl")
include("ql.jl")
# include("plot.jl")
include("RatFun.jl")

spectralmeasure(a,b) = spectralmeasureRat(a,b)

function spectralmeasureRat(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  c = Fun(C.T.nonnegative,Taylor)
  f = Fun(C*(C'*[1]),Ultraspherical{1}())

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1) ,complexroots(c))))
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
  c = Fun(C.T.nonnegative,Taylor)

  # Compute continuous part of measure
  coeffs = Fun(x->(2/pi)*(1-x.^2)./abs(c(x+im*sqrt(1-x.^2))).^2,Ultraspherical{0}()).coefficients
  μ = Fun(coeffs,JacobiWeight(-.5,-.5,Ultraspherical{0}()))

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1),complexroots(c))))
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
  c = Fun(C.T.nonnegative,Taylor)
  f = Fun(C*(C'*[1]),Ultraspherical{1}())

  # Compute continuous part of measure
  finv = (1./f)
  μ = Fun((2/pi)*finv.coefficients,JacobiWeight(.5,.5,space(finv)))

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1),complexroots(c))))
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
  # Finds C such that C*J = Toeplitz([0,1/2])*C
  C = connectionCoeffsOperator(a,b)
  Tfun = Fun(C.T.nonnegative,Taylor)
  sort(real(map(joukowsky,filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1),complexroots(Tfun)))))
end

#Finds C such that C'(U_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and U_k is Chebyshev U
function connectionCoeffsOperator(a,b)
  n = max(length(a),length(b)+1)
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

# Converts coefficients a^J to coefficients a^D
function applyConversion(J::SymTriToeplitz,D::SymTriToeplitz,v::Vector)
  N = length(v)
  b = zeros(N); b1 = zeros(N); b2 = zeros(N)
  for k = N:-1:1
    # before: b = b_k+1, b1 = b_k+2, (and b2 = b_k+3 is to be forgotten)
    b2 = pad((D-J[k,k]*I)*b,N)/J[k,k+1]-b1*(J[k,k+1]/J[k+1,k+2])
    b2[1] += v[k]
    b2, b1, b = b1, b, b2
    # after: b = b_k, b1 = b_k+1, b2 = b_k+2
  end
  b
end

end  #Module

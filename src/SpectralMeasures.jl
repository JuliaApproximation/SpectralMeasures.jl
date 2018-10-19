module SpectralMeasures
using Base, LinearAlgebra, ApproxFun, RecipesBase, RatFun, BandedMatrices, BlockArrays

import Base: +,-,*,/, getindex
import LinearAlgebra: adjoint, transpose

import ApproxFun: Operator, ToeplitzOperator, DiracSpace, IdentityOperator,
            TridiagonalOperator, setdomain, resizedata!, bandinds, PointSpace,
            BandedMatrix, TimesOperator, SpaceOperator, AbstractCount, UnitCount,
            MatrixSpace, ∞, ℓ⁰, domainspace, rangespace, domain, mul_coefficients,
            ldiv_coefficients, InterlaceOperator


import BlockArrays: nblocks


export spectralmeasure, discreteeigs, principalresolvent, discresolvent, validatedspectrum, spectrum
export connectioncoeffsoperator, applyconversion, SymTriOperator, SymTriPertToeplitz
export connectioncoeffsmatrix
export tripleplot
export freejacobioperator, jacobioperator, ql

include("HessenbergUnitary.jl")
include("PertToeplitz.jl")
include("helper.jl")
include("ql.jl")

function spectralmeasure(a,b)
    TT = promote_type(eltype(a),eltype(b))
  # Chop the a and b down
  a = chop!(a); b = 0.5+chop!(b-0.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(TT,n-length(a))]; b = [b;0.5+zeros(TT,n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = SpectralMeasures.connectioncoeffsoperator(a,b)
  c = Fun(Taylor,C.T.nonnegative)
  f = Fun(C*(C'*[1]),Ultraspherical(1))

  # Check for discrete eigenvalues
  z = sort(real(filter!(z->abs(z)<1 && abs(imag(z)) ≤ 10eps(TT)  && !isapprox(abs(z),1),
                        complexroots(c))))
  if length(z) > 0
     Cmu = SpectralMeasures.connectioncoeffsoperator(a[2:end],b[2:end]) # Technically not Cmu from the paper
     cmu = Fun(Taylor,[0;Cmu.T.nonnegative]/b[1]) # this is cmu from the paper
     cprime = differentiate(c)
     eigs=real(map(joukowsky,z))
     weights = (1 .- 1 ./ z.^2).*(real.(cmu.(z))./(2 .* real.(cprime.(z))))
     p = Fun(DiracSpace(eigs),weights) + Fun(JacobiWeight(0.5,0.5,Ultraspherical(1)),[2/TT(pi)])
     q = Fun(PointSpace(eigs),ones(TT,length(eigs))) + f
     μ = RationalFun(p,q)
  else
    μ = RationalFun(Fun(JacobiWeight(0.5,0.5,Ultraspherical(1)),[2/TT(pi)]),f)
  end
  μ
end

function principalresolvent(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Compute the necessary polynomials
  C = connectioncoeffsoperator(a,b)
  Cmu = connectioncoeffsoperator(a[2:end],b[2:end]) # Technically not Cmu from the paper
  f = Fun((C*(C'*[1])),Ultraspherical(1))
  fmu = Fun(Ultraspherical(1),coefficients(Cmu*((C'*[1]).coefficients[2:end])/b[1]))

  # Return the resolvent
  λ->(2*sqrt(complex(λ-1)).*sqrt(complex(λ+1))-2*λ-extrapolate(fmu,λ))./extrapolate(f,λ)
end

function discresolvent(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Compute the necessary polynomials
  C = SpectralMeasures.connectioncoeffsoperator(a,b)
  Cmu = SpectralMeasures.connectioncoeffsoperator(a[2:end],b[2:end]) # Technically not Cmu from the paper
  c = Fun(Taylor,C.T.nonnegative)
  cmu = Fun(Taylor,[0;Cmu.T.nonnegative]/b[1]) # this is the cmu from the paper

  # Return the rational function
  z->-cmu(z)./c(z)
end

function discreteeigs(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that C*J = Toeplitz([0,1/2])*C
  C = connectioncoeffsoperator(a,b)
  Tfun = Fun(Taylor,C.T.nonnegative)
  sort(real(map(joukowsky,filter!(z->abs(z)<1 && isreal(z) && !isapprox(abs(z),1),complexroots(Tfun)))))
end

#Finds C such that C'(U_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and U_k is Chebyshev U
function connectioncoeffsoperator(a,b)
  n = max(2,length(a),length(b)+1)
  N = 2*n #This is sufficient only because we go from Chebyshev U
  a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]

  elType = eltype(a)
  ToeplitzVec = zeros(elType,N)
  K = BandedMatrix(Zeros{elType}(n,N), (0,N+1))
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
  T = ToeplitzOperator(elType[],chop!(ToeplitzVec))
  for j = 1:N
    for i = 1:min(j,N+1-j)
      K[i,j]-=T[i,j]
    end
  end
  T+FiniteOperator(K,ℓ⁰,ℓ⁰)
end
end  #Module

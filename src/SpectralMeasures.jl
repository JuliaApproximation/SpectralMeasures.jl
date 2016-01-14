module SpectralMeasures
    using Base, Compat, ApproxFun, ComplexPhasePortrait, ImageView, ImageMagick

import Base: +,-,*,/

import ApproxFun:BandedOperator,ToeplitzOperator,tridql!,bandinds,DiracSpace, plot, IdentityOperator,
                    TridiagonalOperator,addentries!,setdomain, SavedBandedOperator, resizedata!

export spectralmeasure, spectralmeasureU, spectralmeasureT, ql,SymTriOperator, discreteEigs

include("helper.jl")

joukowsky(z)=.5*(z+1./z)

function discreteEigs(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  Tfun = Fun([C.T[1,1];C.T.negative],Taylor)
  eigs=sort(real(map(joukowsky,filter!(z->abs(z)<1 && isreal(z),complexroots(Tfun)))))
  eigs,C
end

function plotdiscresolvent(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  C = SpectralMeasures.connectionCoeffsOperator(a,b)
  Cmu = SpectralMeasures.connectionCoeffsOperator(a[2:end],b[2:end])
  c = Fun([C.T[1,1];C.T.negative],Taylor)
  cmu = Fun([0;Cmu.T[1,1];Cmu.T.negative]/b[1],Taylor)
  nx = 1000
  x = linspace(-2, 2, nx)
  Z = x' .+ flipdim(x, 1)*im
  portrait(-cmu(Z)./c(Z),PTstepmod)
end

function plotsplitplaneresolvent(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  C = SpectralMeasures.connectionCoeffsOperator(a,b)
  Cmu = SpectralMeasures.connectionCoeffsOperator(a[2:end],b[2:end])
  f = Fun(C'*(C*[1]),Ultraspherical{1}())
  fmu = Fun(Cmu'*((C*[1])[2:end])/b[1],Ultraspherical{1}())
  nx = 1000
  x = linspace(-2, 2, nx)
  Z = x' .+ flipdim(x, 1)*im
  portrait((2*sqrt(Z-1).*sqrt(Z+1)-2*Z-fmu(Z))./f(Z),PTstepmod)
end


function spectralmeasure(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = connectionCoeffsOperator(a,b)
  c = Fun([C.T[1,1];C.T.negative],Taylor)
  f = Fun(C'*(C*[1]),Ultraspherical{1}())
  z = sort(real(filter!(z->abs(z)<1 && isreal(z),complexroots(c))))
  μ1 = Fun((2/pi)*(1./f).coefficients, JacobiWeight(.5,.5,Ultraspherical{1}()))
  if length(z) > 0
    cprime = differentiate(c)
    eigs=real(map(joukowsky,z))
    weights = (z-1./z).^2./(z.*real(cprime(z)).*real(c(1./z)))
    μ1 + Fun(weights,DiracSpace(eigs))
  else
    μ1
  end
end

function spectralmeasureU(a,b;maxlength::Int=10000)
  # a is the first n diagonal elements of J (0 thereafter)
  # b is the first n-1 off-diagonal elements of J (.5 thereafter)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that J*C = C*Toeplitz([0,1/2])
  eigs,C=discreteEigs(a,b)
  # If there is no discrete spectrum
  if isempty(eigs)
    coeffs = forwardSubChebyshevU(C,2*n,[1],1e-15,maxlength)
    #coeffs = linsolve(C,[1];maxlength=maxlength)
    Fun((2/π)*coeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))
  # If there are discrete eigenvalues then we must deflate using QL iteration
  else
    numeigs = length(eigs)
    # Note that output Q is an array of orthogonal operators
    # Q[k] is to be interpretted as having an added kbyk identity in the top left
    eigs,a,b,Q,Qbndwdth = qlIteration(eigs,a,b,1e-15)
    # q0 is the first row of Q where Jnew = Q'*Jold*Q
    # for the spectral measure this is all we need from Q
    q0=[1.0;zeros(Qbndwdth,1)]
    for k=1:numeigs
      q0=[q0[1:k-1];Q[k]'*q0[k:end]]
    end
    # Now let us find the change of basis operator C2 for the continuous part after deflation
    C2 = connectionCoeffsOperator(a,b)
    # note that a and b have been changed (within this function call) after deflation
    n = max(length(a),length(b)+1)

    ctsfactor1 = Fun(C2'*q0[numeigs+1:end],Ultraspherical{1}())
    ctsfactor2 = Fun(forwardSubChebyshevU(C2,2*n,q0[numeigs+1:end],1e-15,maxlength),Ultraspherical{1}())
    #ctsfactor2 = Fun(linsolve(C2,q0[numeigs+1:end];maxlength=maxlength),Ultraspherical{1}())
    ctscoeffs = (ctsfactor1*ctsfactor2).coefficients

    Fun((2/pi)*ctscoeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))+
                    Fun(q0[1:numeigs].^2,DiracSpace(eigs))
  end
end

function spectralmeasureT(a,b;maxlength::Int=10000)
  # a is the first n diagonal elements of J (0 thereafter)
  # b is the first n-1 off-diagonal elements of J (.5 thereafter)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  # Finds C such that J*C = C*Toeplitz([0,1/2])
  eigs,C=discreteEigs(a,b)
  # If there is no discrete spectrum
  if isempty(eigs)
    coeffs = forwardSubChebyshevT(C,2*n,[1],1e-15,maxlength)
    Fun((1/π)*[coeffs[1];sqrt(2)*coeffs[2:end]],JacobiWeight(-.5,-.5,Ultraspherical{0}()))
  # If there are discrete eigenvalues then we must deflate using QL iteration
  else
    numeigs = length(eigs)
    # Note that output Q is an array of orthogonal operators
    # Q[k] is to be interpretted as having an added kbyk identity in the top left
    eigs,a,b,Q,Qbndwdth = qlIteration(eigs,a,b,1e-15)

    # q0 is the first row of Q where Jnew = Q'*Jold*Q
    # for the spectral measure this is all we need from Q
    q0=[1.0;zeros(Qbndwdth,1)]
    for k=1:numeigs
      q0=[q0[1:k-1];Q[k]'*q0[k:end]]
    end

    # Now let us find the change of basis operator L2 for the continuous part after deflation
    C2 = connectionCoeffsOperator(a,b)
    # note that a and b have been changed (within this function call) after deflation
    n = max(length(a),length(b)+1)

    STransCropped = FiniteOperator((1-sqrt(2))*ones(1,1))+ToeplitzOperator([0.],sqrt(2)*onesAndZeros(length(q0)-numeigs))
    factor1coeffs = STransCropped*(C2'*q0[numeigs+1:end])
    factor1 = Fun([factor1coeffs[1];sqrt(2)*factor1coeffs[2:end]],Ultraspherical{0})
    factor2coeffs = forwardSubChebyshevT(C2,2*n,q0[numeigs+1:end],1e-15,maxlength)
    factor2 = Fun([factor2coeffs[1];sqrt(2)*factor2coeffs[2:end]],Ultraspherical{0})
    ctscoeffs = (factor1*factor2).coefficients

    Fun((1/π)*ctscoeffs,JacobiWeight(-.5,-.5,Ultraspherical{0}()))+Fun(q0[1:numeigs].^2,DiracSpace(eigs))
  end
end


spectralmeasure(a...;opts...)=spectralmeasureT(a...;opts...) # default to T


# Note that output Q is an array of orthogonal operators
# Q[k] is to be interpretted as having an added kbyk identity in the top left
function qlIteration(eigs,a,b,tol)
  numeigs = length(eigs)
  Q=Array(BandedOperator{Float64},0)
  Qbndwdth = 0
  for k=1:numeigs
    t1,t0=0.5,eigs[k]
    thisQ = IdentityOperator()
    while b[1]>tol
      Qtmp,Ltmp=ql(a-t0,b,-t0,t1)
      LQ=Ltmp*Qtmp
      # Each QL step increases the size of the compact perturbation by 1, hence the +1 below
      a=Float64[LQ[k,k] for k=1:length(a)+1]+t0
      b=Float64[LQ[k,k+1] for k=1:length(b)+1]
      Qbndwdth += 1
      thisQ = thisQ*Qtmp
    end
    push!(Q,thisQ)
    # Note down the improved value of the eigenvalue and deflate
    eigs[k]=a[1]
    a=a[2:end]
    b=b[2:end]
  end
  eigs,a,b,Q,Qbndwdth
end

function onesAndZeros(n)
  v = ones(n)
  for i = 1:div(n,2)
    v[2*i] = 0
  end
  v
end

function forwardSubChebyshevU(C,bandwidth,f,tol,maxlength)
  # I want to make this a functional. No time to work out how to.
  # f = CompactFunctional(f,)
  if 3+bandwidth > length(f)
    g = [f;zeros(3+bandwidth-length(f))]
  else
    g = f
  end
  m = length(g)

  # I also want  y to be a functional
  y = zeros(2)
  y[1] = g[1]/C[1,1]
  y[2] = (g[2]-C[2,1]*g[1])/C[2,2]

  k = 3
  cvg = false
  while cvg == false
    if k+bandwidth > m
      push!(g,0)
      m+=1
    end
    #There appears to be some issue with returning submatrices of C that sometimes gives Vectors and sometimes Matrix-es.
    tmp = 0
    for i = max(1,k-bandwidth):k-1
      tmp += C[k,i]*y[i]
    end
    yk = (g[k]-tmp)/C[k,k]
    push!(y,yk)
    # If the forward L^2 error is less than tol or we reach max length, stop
    if k > bandwidth
      if ((norm(g[k+1:k+bandwidth]-C[k+1:k+bandwidth,k-bandwidth:k]*y[k-bandwidth:k]) < tol) | (k == maxlength))
        cvg = true
      end
    end
    k = k+1
  end
  y=chop!(y)
  y
end

# Adaptively solves CSx=f by forward substitution, where C is lower triangular
# and S is the (unbounded lower triangular) change of basic matrix from Chebyshev T to Chebyshev U
function forwardSubChebyshevT(Cin,bandwidth,f,tol,maxlength)
    C=SavedBandedOperator(Cin)  # This avoids recalculating C from scrate.
    datalength=300
    resizedata!(C,300)



  # I want to make this a functional. No time to work out how to.
  # f = CompactFunctional(f,)
  g = [f;0]
  m = length(g)

  # I also want x and y to be functionals
  y = zeros(2)
  x = zeros(2)
  y[1] = g[1]/C[1,1]
  x[1] = g[1]/C[1,1]
  y[2] = (g[2]-C[2,1]*g[1])/C[2,2]
  x[2] = y[2]/sqrt(2)

  k = 3
  cvg = false
  while cvg == false
    if k > m
      push!(g,0)
    end

    if k > datalength
        datalength*=2
        resizedata!(C,datalength)
    end

    #There appears to be some issue with returning submatrices of L that sometimes gives Vectors and sometimes Matrix-es.
    tmp = 0
    for i = max(1,k-bandwidth):k-1
      tmp += C[k,i]*y[i]
    end
    yk = (g[k]-tmp)/C[k,k]
    xk = (yk-y[k-2])/sqrt(2)
    push!(y,yk)
    push!(x,xk)
    if k > 20
      if norm(x[k-20:k]) < tol
        cvg = true
      end
      if k == maxlength
        cvg = true
      end
    end
    k = k+1
  end
  x=chop!(x)
  x
end

function connectionCoeffsOperator(a,b)
  n = max(length(a),length(b)+1)
  N = 2*n+1
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

function jacobimatrix(a,b,t0,t1,N)
    J = zeros(N,N)
    a = [a;t0*ones(N-length(a))]
    b = [b;t1*ones(N-length(b))]
    J[1,1]=a[1]
    for i = 1:N-1
        J[i+1,i+1] = a[i+1]
        J[i,i+1] = b[i]
        J[i+1,i] = b[i]
    end
    J
end

jacobimatrix(a,b,N) = jacobimatrix(a,b,0,.5,N)

function jacobioperator(a,b,t0,t1)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  SymTriToeplitz(ToeplitzOperator([t1],[t0,t1]),SymTriOperator(a-t0,b-t1))
end

jacobioperator(a,b) = jacobioperator(a,b,0,.5)

immutable ToeplitzGivens <: BandedOperator{Float64}
    c::Float64
    s::Float64
end

function givenstail(t0::Real,t1::Real)
    @assert t0^2-4t1^2≥0
    s∞ = (t0 - sqrt(t0^2-4t1^2))/(2t1)
    l0 = (t0 + sqrt(t0^2-4t1^2))/2
    if s∞^2 > 1
        s∞ = (t0 + sqrt(t0^2-4t1^2))/(2t1)
        l0 = (t0 - sqrt(t0^2-4t1^2))/2
    end
    c∞ = -sqrt(1-s∞^2)
    α = t1*c∞
    β = c∞*t0 - s∞*α
    l1 = 2t1
    l2 = t1*s∞
    ToeplitzGivens(c∞,s∞),ToeplitzOperator([l1,l2],[l0]),α,β
end

#bandinds(T::ToeplitzGivens)=-ceil(Int,(-36-2log(abs(c)))/log(abs(s))),1
bandinds(T::ToeplitzGivens)=floor(Int,36/log(abs(T.s))),1

function ToeplitzOperator(T::ToeplitzGivens)
    c,s=T.c,T.s
    nonneg=[c^2,s]
    m=-bandinds(T,1)
    if m ≥ 1
        neg=Array(Float64,m)
        neg[1]=-s*nonneg[1]
        for k=2:m
            neg[k]=-s*neg[k-1]
        end
    else
        neg=[]
    end
    ToeplitzOperator(neg,nonneg)
end

addentries!(T::ToeplitzGivens,A,kr::Range,::Colon)=addentries!(ToeplitzOperator(T),A,kr,:)

# This produces an orthogonal operator that is Toeplitz + compact (input is c and s)
function partialgivens(TG::ToeplitzGivens,m)
    T=ToeplitzOperator(TG)
    K=zeros(Float64,m-bandinds(T,1),m)
    neg,nonneg=T.negative,T.nonnegative
    for j=1:m-1
        if j > 1
            K[j-1,j]-=nonneg[2]
        end
        K[j,j]+=1-nonneg[1]
        for k=1:length(neg)
            K[k+j,j]-=neg[k]
        end
    end

    K[m-1,m]-=nonneg[2]
    c,s=TG.c,TG.s
    ret=c
    K[m,m]=c-nonneg[1]
    for k=1:length(neg)
        ret*=-s
        K[k+m,m]=ret-neg[k]
    end
    T+FiniteOperator(K)
end

function ql(a,b,t0,t1)
    @assert t0^2>=4t1^2
    # The Givens rotations coming from infinity (with parameters c∞ and s∞) leave us with the almost triangular
    # a[n-1]  b[n-1]   0    0    0
    # b[n-1]   a[n]   t1    0    0
    #   0       α      β    0    0
    #   0      l2     l1   l0    0
    #   0       0     l2   l1   l0

    TQ,TL,α,β=givenstail(t0,t1)

    # Here we construct this matrix as L
    n = max(length(a),length(b)+1)
    L = jacobimatrix(a,b,t0,t1,n+1)
    L[n,n+1] = t1
    #    L[n+1,n+2] = 0
    L[n+1,n+1]=β
    L[n+1,n]=α

    Q,L=tridql!(L)

    for k=1:size(Q,1)
        Q[k,k]-=1
    end
    for j=1:n+1
        L[j,j]-=TL.nonnegative[1]
        if j ≤ n
            L[j+1,j]-=TL.negative[1]
            if j ≤ n-1
                L[j+2,j]-=TL.negative[2]
            end
        end
    end

    partialgivens(TQ,n+1)*(I+FiniteOperator(Q)),TL+FiniteOperator(L)
end


include("PertToeplitz.jl")

end  #Module

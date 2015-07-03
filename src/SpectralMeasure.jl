module SpectralMeasure
    using Base, Compat, ApproxFun

import ApproxFun:BandedOperator,ToeplitzOperator,tridql!,bandinds,DiracSpace, plot, IdentityOperator,
                    TridiagonalOperator,addentries!,setdomain

export spectralmeasure,ql,SymTriOperator


include("helper.jl")

joukowsky(z)=.5*(z+1./z)

function spectralmeasure(a,b;maxlength::Int=100000)
    # a is the first n diagonal elements of J (0 thereafter)
    # b is the first n-1 off-diagonal elements of J (.5 thereafter)
    # pad to be correct lengths
    if length(b)+1!= length(a)
        n=max(length(b)+1,length(a))
        newb=pad(b,n-1)
        for k=length(b)+1:n-1
            newb[k]=0.5  # pad with toeplitz party
        end
        b=newb
        a=pad(a,n)
    end

    # Finds T,K such that L = T+K, where L takes LJL^{-1} = Toeplitz([0,1/2])
    # This is purely for determining discrete eigenvalues
    T,K=tkoperators(a,b)
    Tfun = Fun([T[1,1];T.negative],Taylor)
    eigs=sort!(real(map!(joukowsky,filter!(z->abs(z)<1 && isreal(z),complexroots(Tfun)))))

    if isempty(eigs)
      #####
      # Old style that just used the L matrix we already have from above
      L=T+K
      coeffs = linsolve(L,[1];maxlength=maxlength)
      Fun((2/π)*coeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))
      #####

      #THIS IS REALLY INEFFICIENT. I'M WORKING ON A FAST VERSION (MW)
      #L = conversionLmatrix([0,0,0],[1/sqrt(2),.5],a,b,maxlength)
      #coeffs = chop(L[1:maxlength,1],10eps())
      #Fun((1/π)*[coeffs[1];sqrt(2)*coeffs[2:end]],JacobiWeight(-.5,-.5,Ultraspherical{0}()))

    # If there are discrete eigenvalues then we must deflate using QL iteration
    else
        numeigs = length(eigs)
        Q=IdentityOperator()
        for k=1:numeigs
            t1,t0=0.5,eigs[k]
            while b[1]>10eps()
              Qtmp,Ltmp=ql(a-t0,b,-t0,t1)
              LQ=Ltmp*Qtmp
              # Each QL step increases the size of the compact perturbation by 1, hence the +1 below
              a=Float64[LQ[k,k] for k=1:length(a)+1]+t0
              b=Float64[LQ[k,k+1] for k=1:length(b)+1]
              Q = Q*Qtmp
            end
            # Note down the improved value of the eigenvalue and deflate
            eigs[k]=a[1]
            a=a[2:end]
            b=b[2:end]
        end

        # q0 is the first row of Q where Jnew = Q'*Jold*Q
        # for the spectral measure this is all we need from Q
        q0 = Q[1,1:bandinds(Q,2)+1]

        # Now let us find the change of basis operator L2 for the continuous part after deflation
        # note that a and b have been changed (within this function call) after deflation
        T,K = tkoperators(a,b)
        L2 = T+K

        ctsfactor1 = Fun(L2'*q0[numeigs+1:end],Ultraspherical{1}())
        ctsfactor2 = Fun(linsolve(L2,q0[numeigs+1:end];maxlength=maxlength),Ultraspherical{1}())
        ctscoeffs = (ctsfactor1*ctsfactor2).coefficients

        Fun([q0[1:numeigs].^2;(2/pi)*ctscoeffs],DiracSpace(JacobiWeight(.5,.5,Ultraspherical{1}()),eigs))
    end
end

function tkoperators(a,b)
    @assert length(a)-length(b)==1
    n = length(a)
    L = conversionLmatrix(a,b,2n)

    T=ToeplitzOperator(vec(L[2*n,2*n-1:-1:1]),[L[2*n,2*n]])
    K = zeros(2*n,2*n)
    for i = 1:2*n
        for j = 1:i
            K[i,j] = L[i,j]-T[i,j]
        end
    end
    K = CompactOperator(K)
    T,K
end

function conversionLmatrix(a,b,c,d,N)
  @assert N >= 2*max(length(a),length(b),length(c),length(d))
  a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
  c = [c;zeros(N-length(c))]; d = [d;.5+zeros(N-length(d))]

  L = zeros(N,N)
  L[1,1] = 1
  L[2,1] = (c[1]-a[1])/b[1]
  L[2,2] = d[1]/b[1]
  for i = 3:N
    L[i,1] = ((c[1]-a[i-1])*L[i-1,1] + d[1]*L[i-1,2] - b[i-2]*L[i-2,1])/b[i-1]
    for j = 2:i-1
      L[i,j] = (d[j-1]*L[i-1,j-1] + (c[j]-a[i-1])*L[i-1,j] + d[j]*L[i-1,j+1] - b[i-2]*L[i-2,j])/b[i-1]
    end
    L[i,i] = d[i-1]*L[i-1,i-1]/b[i-1]
  end
  L
end

# This is for Chebyshev U
conversionLmatrix(a,b,N) = conversionLmatrix(a,b,[],[],N)

function jacobimatrix(a,b,t0,t1,N)
    J = zeros(N,N)
    J[1,1]=a[1]
    a = [a;t0*ones(N-length(a))]
    b = [b;t1*ones(N-length(b))]
    for i = 1:N-1
        J[i+1,i+1] = a[i+1]
        J[i,i+1] = b[i]
        J[i+1,i] = b[i]
    end
    J
end

jacobimatrix(a,b,N) = jacobimatrix(a,b,0,.5,N)

function jacobioperator(a,b,t0,t1)
    Δ=ToeplitzOperator([t1],[t0,t1])
    B=jacobimatrix(a-t0,b-t1,t0,t1,max(length(a),length(b)+1))
    Δ+CompactOperator(B)
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

addentries!(T::ToeplitzGivens,A,kr::Range)=addentries!(ToeplitzOperator(T),A,kr)

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
    T+CompactOperator(K)
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
    n = length(a)
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

    partialgivens(TQ,n+1)*(I+CompactOperator(Q)),TL+CompactOperator(L)
end


include("PertToeplitz.jl")

end  #Module




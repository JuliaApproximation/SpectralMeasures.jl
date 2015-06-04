module SpectralMeasure
    using Base, Compat, ApproxFun

import ApproxFun:BandedOperator,ToeplitzOperator,tridql!,bandinds,DiracSpace, plot

export spectralmeasure

joukowsky(z)=.5*(z+1./z)

function spectralmeasure(a,b;maxlength::Int=10000)
    # a is the first n diagonal elements of J (0 thereafter)
    # b is the first n-1 off-diagonal elements of J (.5 thereafter)

    # Finds T,K such that L = T+K, where L takes LJL^{-1} = Toeplitz([0,1/2])
    # This is purely for determining discrete eigenvalues
    T,K=tkoperators(a,b)
    Tfun = Fun([T[1,1];T.negative],Taylor)
    eigs=sort!(real(map!(joukowsky,filter!(z->abs(z)<1 && isreal(z),complexroots(Tfun)))))

    if isempty(eigs)
      #####
      # Old style that just used the L matrix we already have from above
      #L=T+K
      #coeffs = L\[1]
      #Fun((2/π)*coeffs,JacobiWeight(.5,.5,Ultraspherical{1}()))
      #####

      #THIS IS REALLY INEFFICIENT. I'M WORKING ON A FAST VERSION (MW)
      L = conversionLmatrix([0,0,0],[1/sqrt(2),.5],a,b,maxlength)
      coeffs = chop(L[1:maxlength,1],10eps())
      Fun((1/π)*[coeffs[1];sqrt(2)*coeffs[2:end]],JacobiWeight(-.5,-.5,Ultraspherical{0}()))

    # If there are discrete eigenvalues then we must deflate using QL iteration
    else
        numeigs = length(eigs)
        Q=Array(BandedOperator{Float64},0)
        for k=1:numeigs
            t1,t0=0.5,eigs[k]
            # TODO replace this 1-step QL iteration with iteration until LQ[1,2]<ϵ
            Q1,L1=ql(a-t0,b,-t0,t1)
            push!(Q,Q1)
            LQ=L1*Q1

            a=Float64[LQ[k,k] for k=2:length(a)+1]+t0;
            b=Float64[LQ[k,k+1] for k=2:length(a)];
            # for testing purposes: @assert abs(LQ[1,2]) ≤ 10eps()
            if abs(LQ[1,2]) > 10eps()
                println("QL not converged after 1 step.")
            end
            eigs[k]=LQ[1,1]+t0
        end

        # note that a and b have been changed (within this function call) after deflation
        T,K = tkoperators(a,b)
        L2 = T+K

        q0=[1.0]
        for k=1:numeigs
            q0=[q0[1:k-1];Q[k]'*q0[k:end]]
        end
        coeffs=linsolve(L2,q0[numeigs+1:end];maxlength=maxlength)
        c=q0[numeigs+1]
        Fun([q0[1:numeigs].^2;(2c/π)*coeffs],DiracSpace(JacobiWeight(0.5,0.5,Ultraspherical{1}()),eigs))
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

end  #Module

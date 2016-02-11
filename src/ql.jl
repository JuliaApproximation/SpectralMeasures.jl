

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


immutable HessenbergOrthogonal{T} <: BandedOperator{Float64}
    c::Vector{T}
    s::Vector{T}
    c∞::T
    s∞::T
end




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
    if t0^2<4t1^2
        error("A QL decomposition only exists outside the continuous spectrum")
    end
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



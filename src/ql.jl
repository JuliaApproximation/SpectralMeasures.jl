

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


immutable HessenbergOrthogonal{uplo,T} <: BandedOperator{Float64}
    c::Vector{T}
    s::Vector{T}
    c∞::T
    s∞::T
    band::Int
end

function HessenbergOrthogonal(uplo::Char,c,s,c∞,s∞,band)
    @assert uplo=='L' || uplo=='U'
    HessenbergOrthogonal{uplo,promote_type(eltype(c),eltype(s),
                                           typeof(c∞),typeof(s∞))}(c,s,c∞,s∞,band)
end

function HessenbergOrthogonal(uplo::Char,c,s,c∞,s∞)
    @assert isapprox(s∞^2+c∞^2,1)
    @assert length(c)==length(s)+1

    @assert isapprox(abs(first(c)),1)

    for (cc,ss) in zip(c[2:end],s)
        @assert isapprox(cc^2+ss^2,1)
    end

    band=0
    n=length(s)

    cur=c[1]*c[2]
    tol=eps()


    # Compute the bandwidth of the matrix
    k=1
    for j=1:n+2
        while abs(cur) > tol
            cur*=k≤n?s[k]:s∞
            k+=1
            band+=1
        end
        # increase column and row by one
        if band==0
            # we don't need to divide or multiply by s
            if j≤n-1
                cur=c[j+1]*c[j+2]
            elseif j==n
                cur=c[j+1]*c∞
            else
                cur=c∞^2
            end
        else
            if j≤n-1
                cur*=(k≤n?s[k]:s∞)*c[j+2]/(c[j]*s[j])
            elseif j==n
                cur*=(k≤n?s[k]:s∞)*c∞/(c[j]*s[j])
            elseif j==n+1
                cur*=(k≤n?s[k]:s∞)*c∞/(c[j]*s∞)
            else
                cur*=s∞*c∞/(c∞*s∞)
            end
        end
        k+=1
    end


    HessenbergOrthogonal(uplo,c,s,c∞,s∞,band)
end


Base.ctranspose(Q::HessenbergOrthogonal{'L'})=HessenbergOrthogonal('U',Q.c,Q.s,Q.c∞,Q.s∞,Q.band)
Base.ctranspose(Q::HessenbergOrthogonal{'U'})=HessenbergOrthogonal('L',Q.c,Q.s,Q.c∞,Q.s∞,Q.band)

Base.transpose{uplo,T<:Real}(Q::HessenbergOrthogonal{uplo,T})=Q'


bandinds(Q::HessenbergOrthogonal{'L'})=-Q.band,1
bandinds(Q::HessenbergOrthogonal{'U'})=-1,Q.band


hc(Q::HessenbergOrthogonal,k)=k≤length(Q.c)?Q.c[k]:Q.c∞
hs(Q::HessenbergOrthogonal,k)=k≤length(Q.s)?Q.s[k]:Q.s∞

function addentries!(Q::HessenbergOrthogonal{'L'},A,kr::UnitRange,::Colon)
    sn=length(Q.s)
    cn=length(Q.c)
    b=bandinds(Q,1)
    f=first(kr)
    l=last(kr)


    for k=kr
        A[k,k+1]-=hs(Q,k)
    end
    for j=max(1,first(kr)+b):last(kr)
        col0=hc(Q,j)*hc(Q,j+1)
        for k=0:f-j-1
            col0*=hs(Q,j+k)*hc(Q,j+k+2)/hc(Q,j+k+1)
        end

        for k=max(0,f-j):min(-b,l-j)
            A[j+k,j]+=col0
            col0*=hs(Q,j+k)*hc(Q,j+k+2)/hc(Q,j+k+1)
        end
    end
    A
end

function addentries!(Q::HessenbergOrthogonal{'U'},A,kr::UnitRange,::Colon)
    sn=length(Q.s)
    cn=length(Q.c)
    b=bandinds(Q,2)

    # j is the row here
    for j=kr
        if j≥2
            A[j,j-1]-=hs(Q,j-1)
        end

        col0=hc(Q,j)*hc(Q,j+1)

        for k=0:b
            A[j,j+k]+=col0
            col0*=hs(Q,j+k)*hc(Q,j+k+2)/hc(Q,j+k+1)
        end
    end
    A
end

bandinds(Q::HessenbergOrthogonal{'L'}) = -Q.band, 1
bandinds(Q::HessenbergOrthogonal{'U'}) = -1, Q.band

function *(Q::HessenbergOrthogonal{'L'},v::Vector)
    slen = length(Q.s)
    vlen = length(v)
    if slen < vlen
        s = [Q.s;Q.s∞*ones(vlen-slen+1)]
        c = [Q.c[2:end];Q.c∞*ones(vlen-slen+1)]
        ret = pad(v,vlen+1)
    else
        s = [Q.s;s∞]
        c = [(Q.c)[2:end];c∞]
        ret = pad(v,slen+1)
    end

    N = length(ret) #also equal to length s, c
    ret = (Q.c)[1]*ret
    for i = 1:N-1
        ret[i:i+1] = [c[i] -s[i]; s[i] c[i]]*ret[i:i+1]
    end

    i = N
    # After this point, ret is monotonically decreasing to zero
    while abs(ret[i]) > eps()
        reti = ret[i]
        ret[i] = c[N]*reti
        push!(ret,s[N]*reti)
        i += 1
    end
    ret
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
    c,s,L=tridql!(L)
    Q=HessenbergOrthogonal('L',c,s,TQ.c,-TQ.s)
    for j=1:n+1
        L[j,j]-=TL.nonnegative[1]
        if j ≤ n
            L[j+1,j]-=TL.negative[1]
            if j ≤ n-1
                L[j+2,j]-=TL.negative[2]
            end
        end
    end
    Q,TL+FiniteOperator(L)
end







#### New ql


##
# INPUT: data for givens
# OUTPUT: data for Q and L

function givensrecurrence(α,β,γ0n,γ1n)
    n=length(α)
    @assert length(β)==n

    γ0=Array(Float64,n+1)
    γ1=Array(Float64,n)
    c=Array(Float64,n)
    s=Array(Float64,n)

    l0=Array(Float64,n+1)
    l1=Array(Float64,n)
    l2=Array(Float64,n-1)

    γ0[n+1],γ1[n]=γ0n,γ1n

    l0[n+1]=sqrt(γ0[n+1]^2+β[n]^2)
    c[n]=γ0[n+1]/l0[n+1]
    s[n]=-β[n]/l0[n+1]

    for k=n:-1:1
        l0[k+1]=sqrt(γ0[k+1]^2+β[k]^2)
        c[k]=γ0[k+1]/l0[k+1]
        s[k]=-β[k]/l0[k+1]


        γ0[k]=c[k]*α[k]+s[k]*γ1[k]
        if k>1
            γ1[k-1]=c[k]*β[k-1]
            l2[k-1]=-s[k]*β[k-1]
        end

        l1[k]=c[k]*γ1[k]-s[k]*α[k]
    end


    l0[1]=abs(γ0[1])

    sign(γ0[1]),c,s,l0,l1,l2
end

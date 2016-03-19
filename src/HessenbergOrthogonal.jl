immutable HessenbergOrthogonal{uplo,T} <: BandedOperator{Float64}
    sign::Bool
    c::Vector{T}
    s::Vector{T}
    c∞::T
    s∞::T
    band::Int
end

function HessenbergOrthogonal(uplo::Char,sign,c,s,c∞,s∞,band)
    @assert uplo=='L' || uplo=='U'
    HessenbergOrthogonal{uplo,promote_type(eltype(c),eltype(s),
                                           typeof(c∞),typeof(s∞))}(sign,c,s,c∞,s∞,band)
end

function HessenbergOrthogonal(uplo::Char,sign,c,s,c∞,s∞)
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


    HessenbergOrthogonal(uplo,sign,c,s,c∞,s∞,band)
end


Base.ctranspose(Q::HessenbergOrthogonal{'L'})=HessenbergOrthogonal('U',Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,Q.band)
Base.ctranspose(Q::HessenbergOrthogonal{'U'})=HessenbergOrthogonal('L',Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,Q.band)

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

    si=Q.sign?1:-1
    for k=kr
        A[k,k+1]-=si*hs(Q,k)
    end
    for j=max(1,first(kr)+b):last(kr)
        col0=hc(Q,j)*hc(Q,j+1)
        for k=0:f-j-1
            col0*=hs(Q,j+k)*hc(Q,j+k+2)/hc(Q,j+k+1)
        end

        for k=max(0,f-j):min(-b,l-j)
            A[j+k,j]+=si*col0
            col0*=hs(Q,j+k)*hc(Q,j+k+2)/hc(Q,j+k+1)
        end
    end
    A
end

function addentries!(Q::HessenbergOrthogonal{'U'},A,kr::UnitRange,::Colon)
    sn=length(Q.s)
    cn=length(Q.c)
    b=bandinds(Q,2)

    si=Q.sign?1:-1

    # j is the row here
    for j=kr
        if j≥2
            A[j,j-1]-=si*hs(Q,j-1)
        end

        col0=hc(Q,j)*hc(Q,j+1)

        for k=0:b
            A[j,j+k]+=si*col0
            col0*=hs(Q,j+k)*hc(Q,j+k+2)/hc(Q,j+k+1)
        end
    end
    A
end

function *(Q::HessenbergOrthogonal{'U'},v::Vector)
    si=Q.sign?1:-1

    ret = pad!(si*v,length(v)+1)
    # Compute each Givens rotation starting from the right

    for i = length(v):-1:1
        ret[i],ret[i+1] = hc(Q,i+1)*ret[i] + hs(Q,i)*ret[i+1],
                -hs(Q,i)*ret[i] + hc(Q,i+1)*ret[i+1]
    end
    ret[1]*=hc(Q,1)
    ret
end


function *(Q::HessenbergOrthogonal{'L'},v::Vector)
    N =  max(length(v),length(Q.s))+1
    si=Q.sign?1:-1
    ret = pad!(si*v,N)

    # This part does the computation we are certain we have to do
    ret[1] *= hc(Q,1)
    for i = 1:N-1
        ret[i],ret[i+1] = hc(Q,i+1)*ret[i] -hs(Q,i)*ret[i+1],
                           hs(Q,i)*ret[i] + hc(Q,i+1)*ret[i+1]
    end

    # After this point, ret is monotonically decreasing to zero
    i = N
    while abs(ret[i]) > eps()
        push!(ret,(Q.s∞)*ret[i])
        ret[i] *= Q.c∞
        i += 1
    end
    ret
end



-{uplo}(Q::HessenbergOrthogonal{uplo})=
    HessenbergOrthogonal{uplo,eltype(Q)}(!Q.sign,Q.c,Q.s,
                                          Q.c∞,
                                          Q.s∞,Q.band)



deflate{uplo}(Q::HessenbergOrthogonal{uplo})=HessenbergOrthogonal(uplo,Q.sign,
                                                                  [(Q.sign?1:(-1))*sign(Q.c[1]);Q.c],
                                                                  [0;Q.s],Q.c∞,Q.s∞,Q.band)

deflate(Q::HessenbergOrthogonal,k::Integer)=k==0?Q:deflate(deflate(Q),k-1)
abstract type UnitaryOperator{T} <: Operator{T} end

Base.inv(Q::UnitaryOperator) = Q'
Base.transpose(Q::UnitaryOperator{T}) where {T<:Real} = Q'


Ac_mul_B_coefficients(Q::UnitaryOperator,v::AbstractVector;opts...) =
    A_mul_B_coefficients(Q',v;opts...)


A_ldiv_B_coefficients(Q::UnitaryOperator,v::AbstractVector;opts...) =
    Ac_mul_B_coefficients(Q,v;opts...)



struct HessenbergUnitary{uplo,T} <: UnitaryOperator{T}
    sign::Bool
    c::Vector{T}
    s::Vector{T}
    c∞::T
    s∞::T
    band::Int

    function HessenbergUnitary{uplo,T}(sgn::Bool,c::Vector{T},s::Vector{T},
                                       c∞::T,s∞::T,bnd::Int) where {uplo,T}
        @assert isapprox(s∞^2+c∞^2,1)
        @assert length(c)==length(s)+1

        @assert isapprox(abs(first(c)),1)

        @inbounds for k=2:length(c)
            @assert -100eps() ≤ c[k]^2+s[k-1]^2-1 ≤ 100eps()
        end

        new(sgn,c,s,c∞,s∞,bnd)
    end
end

function HessenbergUnitary(::Type{Val{uplo}},sign,c,s,c∞,s∞,band) where {uplo}
    @assert uplo=='L' || uplo=='U'
    HessenbergUnitary{uplo,promote_type(eltype(c),eltype(s),
                                           typeof(c∞),typeof(s∞))}(sign,c,s,c∞,s∞,band)
end

function HessenbergUnitary(::Type{Val{uplo}},sign,c,s,c∞,s∞) where {uplo}
    @assert isapprox(s∞^2+c∞^2,1)
    @assert length(c)==length(s)+1

    @assert isapprox(abs(first(c)),1)

    @inbounds for k=2:length(c)
        @assert -100eps() ≤ c[k]^2+s[k-1]^2-1 ≤ 100eps()
    end

    band=0
    n = length(s)

    if n ≥ 1
        cur=c[1]*c[2]
    else
        cur=c[1]*c∞
    end

    tol=eps()


    # Compute the bandwidth of the matrix
    k=1
    @inbounds for j=1:n+2
        while abs(cur) > tol
            cur*=k≤n ? s[k] : s∞
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
                cur *= (k≤n ? s[k] : s∞)*c[j+2]/(c[j]*s[j])
            elseif j==n
                cur *= (k≤n ? s[k] : s∞)*c∞/(c[j]*s[j])
            elseif j==n+1
                cur *= (k≤n ? s[k] : s∞)*c∞/(c[j]*s∞)
            else
                cur *= s∞*c∞/(c∞*s∞)
            end
        end
        k+=1
    end


    HessenbergUnitary(Val{uplo},sign,c,s,c∞,s∞,band)
end


Base.ctranspose(Q::HessenbergUnitary{'L',T}) where {T<:Real} =
    HessenbergUnitary(Val{'U'},Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,Q.band)
Base.ctranspose(Q::HessenbergUnitary{'U',T}) where {T<:Real} =
    HessenbergUnitary(Val{'L'},Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,Q.band)



bandinds(Q::HessenbergUnitary{'L'}) = -Q.band,1
bandinds(Q::HessenbergUnitary{'U'}) = -1,Q.band

domainspace(::HessenbergUnitary) = ℓ⁰
rangespace(::HessenbergUnitary) = ℓ⁰




hc(c,c∞,k) = k≤length(c)?c[k]:c∞
hs(s,s∞,k) = k≤length(s)?s[k]:s∞


hc(Q::HessenbergUnitary,k) = hc(Q.c,Q.c∞,k)
hs(Q::HessenbergUnitary,k) = hs(Q.s,Q.s∞,k)

getindex(Q::HessenbergUnitary{'L'},k::Integer,j::Integer) =
    hessuni_getindex(Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,j,k)


getindex(Q::HessenbergUnitary{'U'},k::Integer,j::Integer) =
    hessuni_getindex(Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,k,j)

function hessuni_getindex(sgn::Bool,c::AbstractVector{T},s::AbstractVector{T},
                             c∞::T,s∞::T,
                             k::Integer,j::Integer) where T
    si=sgn?1:-1

    if k>j+1
        zero(T)
    elseif k≥2 && j ==k-1
        -si*hs(s,s∞,k-1)
    else
        col0=hc(c,c∞,k)*hc(c,c∞,k+1)
        for p=k+1:j
            col0*=hs(s,s∞,p-1)*hc(c,c∞,p+1)/hc(c,c∞,p)
        end

        si*col0
    end
end

function A_mul_B_coefficients(Q::HessenbergUnitary{'U'},v::Vector;opts...)
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


function A_mul_B_coefficients(Q::HessenbergUnitary{'L'},v::Vector;tolerance=eps())
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
    while abs(ret[i]) > tolerance
        push!(ret,(Q.s∞)*ret[i])
        ret[i] *= Q.c∞
        i += 1
    end
    ret
end



-(Q::HessenbergUnitary{uplo}) where {uplo}=
    HessenbergUnitary{uplo,eltype(Q)}(!Q.sign,Q.c,Q.s,
                                          Q.c∞,
                                          Q.s∞,Q.band)



deflate(Q::HessenbergUnitary{uplo}) where {uplo}=HessenbergUnitary(Val{uplo},Q.sign,
                                                                  [(Q.sign?1:(-1))*sign(Q.c[1]);Q.c],
                                                                  [0;Q.s],Q.c∞,Q.s∞,Q.band)

deflate(Q::HessenbergUnitary,k::Integer)=k==0?Q:deflate(deflate(Q),k-1)


struct BandedUnitary{uplo,T} <: UnitaryOperator{T}
    ops::Vector{HessenbergUnitary{uplo,T}}
end


Base.ctranspose(Q::BandedUnitary) = BandedUnitary(reverse!(map(ctranspose,Q.ops)))

getindex(Q::BandedUnitary,k::Integer,j::Integer) = TimesOperator(Q.ops)[k,j]
bandinds(Q::BandedUnitary) = bandinds(TimesOperator(Q.ops))


domainspace(::BandedUnitary) = ℓ⁰
rangespace(::BandedUnitary) = ℓ⁰


function A_mul_B_coefficients(Q::BandedUnitary,v::Vector)
    ret=v
    for k=length(Q.ops):-1:1
        ret=A_mul_B_coefficients(Q.ops[k],ret)
    end
    ret
end

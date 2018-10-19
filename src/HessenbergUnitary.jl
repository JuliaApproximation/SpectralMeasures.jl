abstract type UnitaryOperator{T} <: Operator{T} end

inv(Q::UnitaryOperator) = Q'
transpose(Q::UnitaryOperator{T}) where {T<:Real} = Q'

ldiv_coefficients(Q::UnitaryOperator,v::AbstractVector;opts...) =
    mul_coefficients(Q',v;opts...)



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

# we choose colstop based on s for simplicity
function hessenberg_colstop(s, s∞, j)
    n = length(s)
    tol = eps()

    # Compute the bandwidth of the matrix
    k=j
    cur = one(eltype(s))
    while abs(cur) > tol
        cur *= k ≤ n ? s[k] : s∞
        k += 1
    end
    k
end


function HessenbergUnitary(::Type{Val{uplo}},sign, c, s, c∞, s∞) where {uplo}
    @assert isapprox(s∞^2+c∞^2,1)
    @assert length(c)==length(s)+1

    @assert isapprox(abs(first(c)),1)

    @inbounds for k=2:length(c)
        @assert -100eps() ≤ c[k]^2+s[k-1]^2-1 ≤ 100eps()
    end

    band = mapreduce(j -> hessenberg_colstop(s, s∞, j)-j, max, 1:length(c))

    HessenbergUnitary(Val{uplo},sign,c,s,c∞,s∞,band)
end


adjoint(Q::HessenbergUnitary{'L',T}) where {T<:Real} =
    HessenbergUnitary(Val{'U'},Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,Q.band)
adjoint(Q::HessenbergUnitary{'U',T}) where {T<:Real} =
    HessenbergUnitary(Val{'L'},Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,Q.band)



bandwidths(Q::HessenbergUnitary{'L'}) = Q.band,1
bandwidths(Q::HessenbergUnitary{'U'}) = 1,Q.band

domainspace(::HessenbergUnitary) = ℓ⁰
rangespace(::HessenbergUnitary) = ℓ⁰




hc(c,c∞,k) = k≤length(c) ? c[k] : c∞
hs(s,s∞,k) = k≤length(s) ? s[k] : s∞


hc(Q::HessenbergUnitary,k) = hc(Q.c,Q.c∞,k)
hs(Q::HessenbergUnitary,k) = hs(Q.s,Q.s∞,k)

getindex(Q::HessenbergUnitary{'L'},k::Integer,j::Integer) =
    hessuni_getindex(Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,j,k)


getindex(Q::HessenbergUnitary{'U'},k::Integer,j::Integer) =
    hessuni_getindex(Q.sign,Q.c,Q.s,Q.c∞,Q.s∞,k,j)

function hessuni_getindex(sgn::Bool,c::AbstractVector{T},s::AbstractVector{T},
                             c∞::T,s∞::T,
                             k::Integer,j::Integer) where T
    si = sgn ? 1 : -1

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

function mul_coefficients(Q::HessenbergUnitary{'U'},v::Vector;opts...)
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


function mul_coefficients(Q::HessenbergUnitary{'L'},v::Vector;tolerance=eps())
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


adjoint(Q::BandedUnitary) = BandedUnitary(reverse!(map(adjoint,Q.ops)))

getindex(Q::BandedUnitary,k::Integer,j::Integer) = TimesOperator(Q.ops)[k,j]
bandwidths(Q::BandedUnitary) = bandwidths(TimesOperator(Q.ops))


domainspace(::BandedUnitary) = ℓ⁰
rangespace(::BandedUnitary) = ℓ⁰


function mul_coefficients(Q::BandedUnitary,v::Vector)
    ret=v
    for k=length(Q.ops):-1:1
        ret=mul_coefficients(Q.ops[k],ret)
    end
    ret
end

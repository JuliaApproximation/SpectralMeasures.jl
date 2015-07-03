


immutable SymTriOperator{T} <: TridiagonalOperator{T}
    dv::Vector{T}
    ev::Vector{T}
end
SymTriOperator(A::Vector,B::Vector)=SymTriOperator{promote_type(eltype(A),eltype(B))}(A,B)


function addentries!(S::SymTriOperator,A,kr::Range)
   for k=kr
        if k ≤ length(S.dv)
            A[k,k]+=S.dv[k]
        end
        if k ≤ length(S.ev)
            A[k,k+1]+=S.ev[k]
        end
        if 2 ≤ k ≤ length(S.ev)+1
            A[k,k-1]+=S.ev[k-1]
        end
    end
    A
end

# Represents a SymTriOperator + Symmetric ToeplitzOperator
immutable SymTriToeplitz{T} <: TridiagonalOperator{T}
    dv::Vector{T}
    ev::Vector{T}
    a::T
    b::T
end

SymTriToeplitz(dv::Vector,ev::Vector,a,b)=SymTriOperator{promote_type(eltype(dv),eltype(dv),typeof(a),typeof(b))}(dv,ev,a,b)

function SymTriToeplitz(T::ToeplitzOperator,K::SymTriOperator)
    @assert bandinds(T)==(-1,1) && issym(T)
    SymTriToeplitz(K.dv,K.ev,T.nonnegative...)
end


function addentries!(S::SymTriToeplitz,A,kr::Range)
    addentries!(SymTriOperator(S.dv,S.ev),A,kr)

   for k=kr
        if 2 ≤ k
            A[k,k-1]+=S.b
        end
        A[k,k+1]+=S.b
        A[k,k]+=S.a
    end
    A
end


immutable PertToeplitz{T} <: BandedOperator{T}
    T::ToeplitzOperator{T}
    K::CompactOperator{T}
end

bandinds(P::PertToeplitz)=min(bandinds(P.T,1),bandinds(P.K,1)),max(bandinds(P.T,2),bandinds(P.K,2))

function addentries!(P::PertToeplitz,A,kr::Range)
    addentries!(P.T,A,kr)
    addentries!(P.K,A,kr)
end


Base.issym(T::ToeplitzOperator)=length(T.negative)==length(T.nonnegative)-1&&T.negative==T.nonnegative[2:end]

# slcie of a PertToeplitz is a PertToeplitz
Base.slice(P::PertToeplitz,kr::FloatRange,jr::FloatRange)=slice(P.T,kr,jr)+slice(P.K,kr,jr)


+(T::ToeplitzOperator,K::CompactOperator)=PertToeplitz(T,K)
+(T::ToeplitzOperator,K::SymTriOperator)=SymTriToeplitz(T,K)
+(K::Union(CompactOperator,SymTriOperator),T::ToeplitzOperator)=T+K

for OP in (:+,:-)
    @eval begin
        $OP(A::SymTriToeplitz,c::UniformScaling)=SymTriToeplitz(A.dv,A.ev,$OP(A.a,c.λ),A.b)
    end
end

+(c::UniformScaling,A::SymTriToeplitz)=SymTriToeplitz(A.dv,A.ev,c.λ+A.a,A.b)
-(c::UniformScaling,A::SymTriToeplitz)=SymTriToeplitz(-A.dv,-A.ev,c.λ-A.a,-A.b)

*(c::Number,A::SymTriToeplitz)=SymTriToeplitz(c*A.dv,c*A.ev,c*A.a,c*A.b)
/(A::SymTriToeplitz,c::Number)=(1/c)*A






ql(A::SymTriToeplitz)=ql(A.dv+A.a,A.ev+A.b,A.a,A.b)

function Base.eigvals(A::SymTriToeplitz)
    if isapprox(A.a,0.) && isapprox(A.b,0.5)
        spectralmeasure(A.dv+A.a,A.ev+A.b)
    elseif isapprox(A.a,0.)
        c=2*A.b
        μ=eigvals(A/c)
        sp=space(μ)
        np=length(sp.points)  # number of points
        Fun([μ.coefficients[1:np];μ.coefficients[np+1:end]/c],setdomain(space(μ),c*domain(μ)))
    else
        μ=eigvals(A-A.a*I)
        setdomain(μ,domain(μ)+A.a)
    end
end

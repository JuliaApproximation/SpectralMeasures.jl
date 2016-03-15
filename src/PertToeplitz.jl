## New ToeplitzOperator routines


ToeplitzOperator(L::UniformScaling)=ToeplitzOperator(eltype(L)[],[L.λ])

Base.issym(T::ToeplitzOperator)=length(T.negative)==length(T.nonnegative)-1&&T.negative==T.nonnegative[2:end]


for OP in (:+,:-)
    @eval begin
        $OP(A::ToeplitzOperator,c::UniformScaling)=$OP(A,ToeplitzOperator(c))
        $OP(c::UniformScaling,A::ToeplitzOperator)=$OP(ToeplitzOperator(c),A)
        function $OP(A::ToeplitzOperator,B::ToeplitzOperator)
            n=max(length(A.negative),length(B.negative))
            m=max(length(A.nonnegative),length(B.nonnegative))
            ToeplitzOperator($OP(pad(A.negative,n),pad(B.negative,n)),
                             $OP(pad(A.nonnegative,m),pad(B.nonnegative,m)))
        end
    end
end



## Symetric tridiagonal finite dimensional operators





immutable SymTriOperator{T} <: TridiagonalOperator{T}
    dv::Vector{T}
    ev::Vector{T}
end
SymTriOperator(A::Vector,B::Vector)=SymTriOperator{promote_type(eltype(A),eltype(B))}(A,B)


function addentries!(S::SymTriOperator,A,kr::Range,::Colon)
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


function SymTriToeplitz(T::ToeplitzOperator)
    @assert issym(T)

    if isdiag(T)
        SymTriToeplitz(eltype(T)[],eltype(T)[],T.nonnegative[1],zero(eltype(T)))
    elseif bandinds(T)==(-1,1)
        SymTriToeplitz(eltype(T)[],eltype(T)[],T.nonnegative...)
    else
        error("Not a tridiagonal operator")
    end
end

function addentries!(S::SymTriToeplitz,A,kr::Range,::Colon)
    addentries!(SymTriOperator(S.dv,S.ev),A,kr,:)

   for k=kr
        if 2 ≤ k
            A[k,k-1]+=S.b
        end
        A[k,k+1]+=S.b
        A[k,k]+=S.a
    end
    A
end

## represents T + K where T is Toeplitz and K is finite-dimensional
immutable PertToeplitz{T} <: BandedOperator{T}
    T::ToeplitzOperator{T}
    K::FiniteOperator{BandedMatrix{T},T}
end

bandinds(P::PertToeplitz)=min(bandinds(P.T,1),bandinds(P.K,1)),max(bandinds(P.T,2),bandinds(P.K,2))

function addentries!(P::PertToeplitz,A,kr::Range,::Colon)
    addentries!(P.T,A,kr,:)
    addentries!(P.K,A,kr,:)
end





# slcie of a PertToeplitz is a PertToeplitz
Base.slice(P::PertToeplitz,kr::FloatRange,jr::FloatRange)=slice(P.T,kr,jr)+slice(P.K,kr,jr)


+(T::ToeplitzOperator,K::FiniteOperator)=PertToeplitz(T,K)
+(T::ToeplitzOperator,K::SymTriOperator)=SymTriToeplitz(T,K)
+(K::Union{FiniteOperator,SymTriOperator},T::ToeplitzOperator)=T+K

for OP in (:+,:-)
    @eval begin
        $OP(A::SymTriToeplitz,c::UniformScaling)=SymTriToeplitz(A.dv,A.ev,$OP(A.a,c.λ),A.b)
    end
end

+(c::UniformScaling,A::SymTriToeplitz)=SymTriToeplitz(A.dv,A.ev,c.λ+A.a,A.b)
-(c::UniformScaling,A::SymTriToeplitz)=SymTriToeplitz(-A.dv,-A.ev,c.λ-A.a,-A.b)

*(c::Number,A::SymTriToeplitz)=SymTriToeplitz(c*A.dv,c*A.ev,c*A.a,c*A.b)
/(A::SymTriToeplitz,c::Number)=(1/c)*A





for OP in (:ql,:(Base.eigvals),:(Base.eig))
    @eval $OP(A::ToeplitzOperator)=$OP(SymTriToeplitz(A))
end


ql(A::SymTriToeplitz)=ql(A.dv+A.a,A.ev+A.b,A.a,A.b)

function Base.eigvals(A::SymTriToeplitz)
    if isapprox(A.a,0.) && isapprox(A.b,0.5)
        spectralmeasure(A.dv+A.a,A.ev+A.b)
    elseif isapprox(A.a,0.)
        c=2*A.b
        μ=eigvals(A/c)
        sp=space(μ)

        np=isa(sp,DiracSpace)?length(sp.points):0  # number of points

        Fun([μ.coefficients[1:np];μ.coefficients[np+1:end]/c],setdomain(space(μ),c*domain(μ)))
    else
        μ=eigvals(A-A.a*I)
        setdomain(μ,domain(μ)+A.a)
    end
end

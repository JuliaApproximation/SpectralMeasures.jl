## New ToeplitzOperator routines


ToeplitzOperator(L::UniformScaling) = ToeplitzOperator(eltype(L)[],[L.λ])

Base.issymmetric(T::ToeplitzOperator) =
    length(T.negative)==length(T.nonnegative)-1&&T.negative==T.nonnegative[2:end]


for OP in (:+,:-)
    @eval begin
        $OP(A::ToeplitzOperator,c::UniformScaling) = $OP(A,ToeplitzOperator(c))
        $OP(c::UniformScaling,A::ToeplitzOperator) = $OP(ToeplitzOperator(c),A)
        function $OP(A::ToeplitzOperator,B::ToeplitzOperator)
            n=max(length(A.negative),length(B.negative))
            m=max(length(A.nonnegative),length(B.nonnegative))
            ToeplitzOperator($OP(pad(A.negative,n),pad(B.negative,n)),
                             $OP(pad(A.nonnegative,m),pad(B.nonnegative,m)))
        end
    end
end



## Symetric tridiagonal finite dimensional operators





struct SymTriOperator{T} <: TridiagonalOperator{T}
    dv::Vector{T}
    ev::Vector{T}
end
SymTriOperator(A::Vector,B::Vector) =
    SymTriOperator{promote_type(eltype(A),eltype(B))}(A,B)


for OP in (:domainspace,:rangespace)
    @eval $OP(::SymTriOperator) = ℓ⁰
end

function getindex(S::SymTriOperator,k::Integer,j::Integer)
    if k ≤ length(S.dv) && k == j
        S.dv[k]
    elseif k ≤ length(S.ev) && j==k+1
        S.ev[k]
    elseif 2 ≤ k ≤ length(S.ev)+1 && j==k-1
        S.ev[k-1]
    else
        zero(eltype(S))
    end
end


*(c::Number,A::SymTriOperator) = SymTriOperator(c*A.dv,c*A.ev)
*(A::SymTriOperator,c::Number) = SymTriOperator(c*A.dv,c*A.ev)

function Base.SymTridiagonal(S::SymTriOperator,kr::UnitRange,jr::UnitRange)
    n=last(kr)
    @assert n==last(jr)
    SymTridiagonal(pad(S.dv,n),pad(S.ev,n-1))
end



# Represents a SymTriOperator + Symmetric ToeplitzOperator
struct SymTriToeplitz{T} <: TridiagonalOperator{T}
    dv::Vector{T}
    ev::Vector{T}
    a::T
    b::T

    SymTriToeplitz{T}(dv::Vector{T},ev::Vector{T},a::T,b::T) where T = new(dv,ev,a,b)
    SymTriToeplitz{T}(dv::Vector,ev::Vector,a,b) where T = new(Vector{T}(dv),Vector{T}(ev),T(a),T(b))
end

SymTriToeplitz(dv::Vector,ev::Vector,a,b) =
    SymTriToeplitz{promote_type(eltype(dv),eltype(dv),typeof(a),typeof(b))}(dv,ev,a,b)

function SymTriToeplitz(T::ToeplitzOperator,K::SymTriOperator)
    @assert bandinds(T)==(-1,1) && issym(T)
    SymTriToeplitz(K.dv+T.nonnegative[1],K.ev+T.nonnegative[2],T.nonnegative...)
end

function SymTriToeplitz(K::SymTriOperator{TT}) where TT
    SymTriToeplitz(K.dv,K.ev,zero(TT),zero(TT))
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


for OP in (:domainspace,:rangespace)
    @eval $OP(::SymTriToeplitz) = ℓ⁰
end

function Base.getindex(S::SymTriToeplitz,kr::UnitCount{Int},jr::UnitCount{Int})
    k=first(kr)
    @assert k==first(jr)

    SymTriToeplitz(S.dv[k:end],S.ev[k:end],S.a,S.b)
end

function getindex(S::SymTriToeplitz,k::Integer,j::Integer)
        if 2 ≤ k && j ==k-1
            k≤length(S.ev)+1?S.ev[k-1]:S.b
        elseif j==k+1
            k≤length(S.ev)?S.ev[k]:S.b
        elseif j==k
            k≤length(S.dv)?S.dv[k]:S.a
        else
            zero(eltype(S))
        end
end



function Base.SymTridiagonal(S::SymTriToeplitz,kr::UnitRange,jr::UnitRange)
    n=last(kr)
    @assert n==last(jr)
    dv= n>length(S.dv) ?
        [S.dv;S.a*ones(eltype(S),n-length(S.dv))] :
        S.dv[1:n]

    ev= n-1>length(S.ev) ?
        [S.ev;S.b*ones(eltype(S),n-1-length(S.ev))] :
        S.ev[1:n-1]

    SymTridiagonal(dv,ev)
end




## represents T + K where T is Toeplitz and K is finite-dimensional
struct PertToeplitz{S} <: Operator{S}
    T::ToeplitzOperator{S}
    K::FiniteOperator{BandedMatrix{S},S}
end

for OP in (:domainspace,:rangespace)
    @eval $OP(::PertToeplitz) = ℓ⁰
end

bandinds(P::PertToeplitz) =
    min(bandinds(P.T,1),bandinds(P.K,1)),max(bandinds(P.T,2),bandinds(P.K,2))

getindex(P::PertToeplitz,k::Integer,j::Integer) =
    P.T[k,j]+P.K[k,j]

getindex(P::PertToeplitz,k::AbstractCount,j::AbstractCount) =
    P.T[k,j]+P.K[k,j]


for OP in (:+,:-)
    @eval begin
        function $OP(A::SymTriToeplitz,B::SymTriToeplitz)
            n_dv = max(length(A.dv),length(B.dv))
            n_ev = max(length(A.ev),length(B.ev))
            SymTriToeplitz($OP([A.dv;fill(A.a,n_dv-length(A.dv))],
                               [B.dv;fill(B.a,n_dv-length(B.dv))]),
                               $OP([A.ev;fill(A.b,n_ev-length(A.ev))],
                                   [B.ev;fill(B.b,n_ev-length(B.ev))]),
                            $OP(A.a,B.a),$OP(A.b,B.b))
        end
        $OP(T::ToeplitzOperator,K::FiniteOperator) = PertToeplitz(T,$OP(K))
        $OP(T::ToeplitzOperator,K::SymTriOperator) = SymTriToeplitz(T,$OP(K))
        $OP(S::SymTriToeplitz,K::SymTriOperator) = $OP(S,SymTriToeplitz(K))
        $OP(K::Union{FiniteOperator,SymTriOperator},T::ToeplitzOperator) = $OP(T)+K

        $OP(A::SymTriToeplitz,c::UniformScaling) =
            SymTriToeplitz($OP(A.dv,c.λ),A.ev,$OP(A.a,c.λ),A.b)
    end
end

-(A::SymTriToeplitz) = SymTriToeplitz(-A.dv,-A.ev,-A.a,-A.b)
+(c::UniformScaling,A::SymTriToeplitz) = SymTriToeplitz(c.λ+A.dv,A.ev,c.λ+A.a,A.b)
-(c::UniformScaling,A::SymTriToeplitz) = SymTriToeplitz(c.λ-A.dv,-A.ev,c.λ-A.a,-A.b)

*(c::Number,A::SymTriToeplitz) = SymTriToeplitz(c*A.dv,c*A.ev,c*A.a,c*A.b)
/(A::SymTriToeplitz,c::Number) = (1/c)*A





for OP in (:ql,:(Base.eigvals),:(Base.eig))
    @eval $OP(A::ToeplitzOperator)=$OP(SymTriToeplitz(A))
end


ql(A::SymTriToeplitz) = ql(A.dv,A.ev,A.a,A.b)

function spectralmeasure(J::SymTriToeplitz)
  μ = spectralmeasure(.5*(J.dv-J.a)/J.b,.5*J.ev/J.b)
  2*J.b*setdomain(μ,domain(μ) + J.a)
end

Base.eigvals(A::SymTriToeplitz)=domain(spectralmeasure(A))



function *(L::PertToeplitz,Q::HessenbergUnitary{'L'})
    n=max(size(L.K.matrix,1),length(Q.s)+3)

    if bandinds(L)==(-2,0)
         # We check if L*Q is tridiagin al
        tol=1E-14*(maximum(L.T)+maximum(L.K))
        istri=true
        for k=3:n
            if abs(L[k,k-2]*hc(Q,k-1)+L[k,k-1]*hs(Q,k-2)*hc(Q,k)+L[k,k]*hs(Q,k-2)*hs(Q,k-1)*hc(Q,k+1))>tol
                istri=false
                break
            end
        end
        if istri
            issym=true
            if !isapprox(-L[1,1]*hs(Q,1),L[2,1]*hc(Q,1)*hc(Q,2)+L[2,2]*hc(Q,1)*hc(Q,3)*hs(Q,1);atol=tol)
                issym=false
            end

            if issym
                for k=2:n+1  # kth row
                    if !isapprox(-L[k+1,k-1]*hs(Q,k-1)+L[k+1,k]*hc(Q,k)*hc(Q,k+1)+L[k+1,k+1]*hc(Q,k)*hc(Q,k+2)*hs(Q,k),
                    -L[k,k]*hs(Q,k);atol=tol)
                        issym=false
                        break
                    end
                end
            end

            if issym
               # result is SymTriToeplitxz

                ev=Array{Float64}(max(min(size(L.K.matrix,1),size(L.K.matrix,2)),
                                     length(Q.s)))
                for k=1:length(ev)
                    ev[k]=-L[k,k]*hs(Q,k)
                end

                dv=Array{Float64}(max(length(Q.s)+1,size(L.K.matrix,1)))
                dv[1]=hc(Q,1)*hc(Q,2)*L[1,1]
                for k=2:length(dv)
                    dv[k]=-hs(Q,k-1)*L[k,k-1]+hc(Q,k)*hc(Q,k+1)*L[k,k]
                end

                t1=-L.T[1,1]*Q.s∞
                t0=-Q.s∞*L.T[2,1]+Q.c∞^2*L.T[1,1]

                si=Q.sign?1:-1
                return SymTriToeplitz(si*dv,si*ev,si*t0,si*t1)
            end
        end
    end

    # default constructor
    TimesOperator(L,Q)
end

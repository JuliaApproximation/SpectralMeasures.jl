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
    SymTriToeplitz(K.dv+T.nonnegative[1],K.ev+T.nonnegative[2],T.nonnegative...)
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

function Base.getindex(S::SymTriToeplitz,kr::FloatRange,jr::FloatRange)
    k=first(kr)
    @assert k==first(jr)
    @assert step(kr)==step(jr)==1
    @assert last(kr)==last(jr)==Inf

    SymTriToeplitz(S.dv[k:end],S.ev[k:end],S.a,S.b)
end

function addentries!(S::SymTriToeplitz,A,kr::Range,::Colon)
   for k=kr
        if 2 ≤ k
            A[k,k-1]+=k≤length(S.ev)+1?S.ev[k-1]:S.b
        end
        A[k,k+1]+=k≤length(S.ev)?S.ev[k]:S.b
        A[k,k]+=k≤length(S.dv)?S.dv[k]:S.a
    end
    A
end

## represents T + K where T is Toeplitz and K is finite-dimensional
immutable PertToeplitz{S} <: BandedOperator{S}
    T::ToeplitzOperator{S}
    K::FiniteOperator{BandedMatrix{S},S}
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
        $OP(A::SymTriToeplitz,c::UniformScaling)=SymTriToeplitz($OP(A.dv,c.λ),A.ev,$OP(A.a,c.λ),A.b)
    end
end

-(A::SymTriToeplitz)=SymTriToeplitz(-A.dv,-A.ev,-A.a,-A.b)
+(c::UniformScaling,A::SymTriToeplitz)=SymTriToeplitz(c.λ+A.dv,A.ev,c.λ+A.a,A.b)
-(c::UniformScaling,A::SymTriToeplitz)=SymTriToeplitz(c.λ-A.dv,-A.ev,c.λ-A.a,-A.b)

*(c::Number,A::SymTriToeplitz)=SymTriToeplitz(c*A.dv,c*A.ev,c*A.a,c*A.b)
/(A::SymTriToeplitz,c::Number)=(1/c)*A





for OP in (:ql,:(Base.eigvals),:(Base.eig))
    @eval $OP(A::ToeplitzOperator)=$OP(SymTriToeplitz(A))
end


ql(A::SymTriToeplitz)=ql(A.dv,A.ev,A.a,A.b)

function spectralmeasure(A::SymTriToeplitz)
    if isapprox(A.a,0.) && isapprox(A.b,0.5)
        spectralmeasure(A.dv,A.ev)
    elseif isapprox(A.a,0.)
        c=2*A.b
        μ=spectralmeasure(A/c)
        sp=space(μ)

        np=isa(sp,DiracSpace)?length(sp.points):0  # number of points

        Fun([μ.coefficients[1:np];μ.coefficients[np+1:end]/c],setdomain(space(μ),c*domain(μ)))
    else
        μ=spectralmeasure(A-A.a*I)
        setdomain(μ,domain(μ)+A.a)
    end
end

Base.eigvals(A::SymTriToeplitz)=domain(spectralmeasure(A))



function *(L::PertToeplitz,Q::HessenbergOrthogonal{'L'})
    n=max(size(L.K.matrix,1),length(Q.s)+3)

    if bandinds(L)==(-2,0)
         # We check if L*Q is tridiagin al
        tol=1E-15
        istri=true
        for k=3:n
            if abs(L[k,k-2]*hc(Q,k-1)+L[k,k-1]*hs(Q,k-2)*hc(Q,k)+L[k,k]*hs(Q,k-2)*hs(Q,k-1)*hc(Q,k+1))>tol
                istri=false
                break
            end
        end
        if istri
            issym=true
            if !isapprox(-L[1,1]*hs(Q,1),L[2,1]*hc(Q,1)*hc(Q,2)+L[2,2]*hc(Q,1)*hc(Q,3)*hs(Q,1))
                issym=false
            end

            if issym
                for k=2:n+1  # kth row
                    if !isapprox(-L[k+1,k-1]*hs(Q,k-1)+L[k+1,k]*hc(Q,k)*hc(Q,k+1)+L[k+1,k+1]*hc(Q,k)*hc(Q,k+2)*hs(Q,k),
                    -L[k,k]*hs(Q,k))
                        issym=false
                        break
                    end
                end
            end

            if issym
               # result is SymTriToeplitxz

                ev=Array(Float64,max(min(size(L.K.matrix,1),size(L.K.matrix,2)),
                                     length(Q.s)))
                for k=1:length(ev)
                    ev[k]=-L[k,k]*hs(Q,k)
                end

                dv=Array(Float64,max(length(Q.s)+1,size(L.K.matrix,1)))
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

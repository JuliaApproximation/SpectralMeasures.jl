
# returns the parameters for the limiting Toeplitz
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
    c∞,s∞,ToeplitzOperator([l1,l2],[l0]),α,β
end




function ql(a,b,t0,t1)
    if t0^2<4t1^2
        error("A QL decomposition only exists outside the essential spectrum")
    end
    # The Givens rotations coming from infinity (with parameters c∞ and s∞) leave us with the almost triangular
    # a[n-1]  b[n-1]   0    0    0
    # b[n-1]   a[n]   t1    0    0
    #   0       α      β    0    0
    #   0      l2     l1   l0    0
    #   0       0     l2   l1   l0

    c∞,s∞,TL,α,β=givenstail(t0,t1)


    if TL[1,1] < 0
        # we want positive on L
        Q,L=ql(-a,-b,-t0,-t1)
        return -Q,L
    end

    # Here we construct this matrix as L
    n = max(length(a),length(b)+1)
    J = jacobimatrix(a,b,t0,t1,n+1)
    J[n,n+1] = t1
    #    L[n+1,n+2] = 0
    J[n+1,n+1]=β
    J[n+1,n]=α
    c,s,L=tridql!(J)


    Q=HessenbergUnitary('L',true,c,s,c∞,-s∞)


    for j=1:n+1
        L[j,j]-=TL.nonnegative[1]
        if j ≤ n
            L[j+1,j]-=TL.negative[1]
            if j ≤ n-1
                L[j+2,j]-=TL.negative[2]
            end
        end
    end
    Q,TL+FiniteOperator(L,ℓ⁰,ℓ⁰)
end

discreteeigs(J::SymTriToeplitz) = 2*J.b*discreteeigs(.5*(J.dv-J.a)/J.b,.5*J.ev/J.b) + J.a

connection_coeffs_operator(J::SymTriToeplitz) = connection_coeffs_operator(.5*(J.dv-J.a)/J.b,.5*J.ev/J.b)



immutable SpectralMap{CC,QQ,RS,T} <: Operator{T}
    C::CC
    Q::QQ
    rangespace::RS
end

SpectralMap{T}(C::Operator{T},Q::Operator{T},rs) =
    SpectralMap{typeof(C),typeof(Q),typeof(rs),T}(C,Q,rs)


domainspace(::SpectralMap) = SequenceSpace()
rangespace(S::SpectralMap) = S.rangespace


A_ldiv_B_coefficients(S::SpectralMap,v::AbstractVector;opts...) =
    A_ldiv_B_coefficients(S.Q,A_ldiv_B_coefficients(S.C,v);opts...)

A_mul_B_coefficients(S::SpectralMap,v::AbstractVector;opts...) =
    A_mul_B_coefficients(S.C,A_mul_B_coefficients(S.Q,v;opts...))

function getindex(S::SpectralMap,k::Integer,j::Integer)
    v = A_mul_B_coefficients(S,[zeros(j-1);1])
    k ≤ length(v) && return v[k]
    zero(eltype(S))
end


isbanded(S::SpectralMap) = true
function bandinds(S::SpectralMap)
    bi = bandinds(S.Q)
    bi[1],bi[2]+bandinds(S.C,2)
end

function Base.eig(Jin::SymTriToeplitz)
    Qret=Array(HessenbergUnitary{'U',Float64},0)
    λapprox=sort(discreteeigs(Jin))

    ctsspec = ApproxFun.Interval(Jin.a-2*abs(Jin.b),Jin.a+2*abs(Jin.b))

    J=Jin

    if length(λapprox) == 0
        C=connection_coeffs_operator(J)

        x=Fun(identity,Ultraspherical(1,ctsspec))

        U=SpaceOperator(C,ℓ⁰,space(x))
        return x,U
    end

    λ=Array(Float64,0)

    tol=1E-14
    for k=1:length(λapprox)
        μ=λapprox[k]

        Q,L=ql(J-μ*I)
        push!(Qret,deflate(Q',k-1))
        J=L*Q+μ*I

         while abs(J[1,2]) > tol
             # μ=J[1,1] DO NOT DO THIS. IF MU IS NOT ACCURATE, J[1,1] CAN BE AN INVALID SHIFT (MW)
             Q,L=ql(J-μ*I)
             J=L*Q+μ*I
             push!(Qret,deflate(Q',k-1))
         end

        push!(λ,J[1,1])
        J=J[2:end,2:end]
    end

    if length(λ) == 1
         Q=Qret[1]

         x=Fun(identity,PointSpace(λ[1])⊕Ultraspherical(1,ctsspec))
         C=SpaceOperator(InterlaceOperator(Diagonal([eye(length(λ)),connection_coeffs_operator(J)])),SequenceSpace(),space(x))

         U=SpectralMap(C,Q,space(x))
         return x,U
    else
        Q=BandedUnitary(reverse!(Qret))
        x=Fun(identity,mapreduce(PointSpace,⊕,λ)⊕Ultraspherical(1,ctsspec))
        C=SpaceOperator(InterlaceOperator(Diagonal([eye(length(λ)),connection_coeffs_operator(J)])),SequenceSpace(),space(x))

        U=SpectralMap(C,Q,space(x))
        return x,U
    end
end

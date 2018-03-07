
# returns the parameters for the limiting Toeplitz
function givenstail(t₀::Real,t₁::Real)
    @assert t₀^2-4t₁^2≥0
    s = (-t₀ + sqrt(t₀^2-4t₁^2))/(2t₁)
    l⁰ = (t₀ + sqrt(t₀^2-4t₁^2))/2
    # if s∞^2 > 1
    #     s∞ = (t₀ + sqrt(t₀^2-4t₁^2))/(2t₁)
    #     l0 = (t₀ - sqrt(t₀^2-4t₁^2))/2
    # end
    c = -sqrt(1-s^2)
    γ¹ = t₁*c
    γ⁰ = c*t₀ + s*γ¹
    l¹ = 2t₁  # = c*γ¹ - st₁
    l² = -t₁*s
    c,s,l⁰,l¹,l²,γ¹,γ⁰
end


ql(a,b,t₀,t₁) = ql!(copy(a),copy(b),t₀,t₁)

function ql!(a,b,t₀,t₁)
    if t₀^2 < 4t₁^2
        error("A QL decomposition only exists outside the essential spectrum")
    end
    # The Givens rotations coming from infinity (with parameters c∞ and s∞) leave us with the almost triangular
    # a[n-1]  b[n-1]   0    0    0
    # b[n-1]   a[n]   t₁    0    0
    #   0       γ¹      γ⁰    0    0
    #   0      l2     l1   l0    0
    #   0       0     l2   l1   l0


    if t₀ < 0
        # we want positive on the diagonal
        Q,L=ql(-a,-b,-t₀,-t₁)
        return -Q,L
    end

    # match dimensions in paper: a has n perts and b has n-1 perts
    n = max(length(a), length(b)+1)
    if n > length(a)
        m = length(a)
        a = resize!(a,n)
        a[m+1:end] .= t₀
    end

    if n > length(b)+1
        m = length(b)
        b = resize!(b,n-1)
        b[m+1:end] .= t₁
    end
    c∞,s∞,l⁰,l¹,l²,γ¹∞,γ⁰∞ = givenstail(t₀,t₁)
    # use recurrence for c. If we have a_0,…,a_N,t0,t0…, then
    # we only need c_-1,c_0,c_1,…,c_{N-1}.
    c=Array{eltype(c∞)}(n)
    s=Array{eltype(c∞)}(n-1)

    # ranges from 1 to N
    γ¹ = Array{eltype(c∞)}(n-1)

    # ranges from 0 to N
    γ⁰ = Array{eltype(c∞)}(n)

    γ⁰[n] = c∞*a[n] + s∞*γ¹∞  # k = N

    if n ≠ 1
        γ¹[n-1] = c∞*b[n-1]  # k = N


        k=n-1
        nrm = 1/sqrt(γ⁰[k+1]^2+b[k]^2)
        c[k+1] = γ⁰[k+1]*nrm  # K = N-1
        s[k] = -b[k]*nrm # K = N-1

        @inbounds for k=n-2:-1:1
            γ¹[k] = c[k+2]*b[k]  # k = N-1
            γ⁰[k+1] = c[k+2]*a[k+1] + s[k+1]*γ¹[k+1]  # k = N
            nrm = 1/sqrt(γ⁰[k+1]^2+b[k]^2)
            c[k+1] = γ⁰[k+1]*nrm  # K = N-1
            s[k] = -b[k]*nrm # K = N-1
        end

        γ⁰[1] = c[2]*a[1] + s[1]*γ¹[1]  # k = 0
    end


    c[1] = sign(γ⁰[1])  # k = -1

    Q = HessenbergUnitary(Val{'L'},true,c,s,c∞,s∞)

    L = BandedMatrix(eltype(c∞),n+1,n,2,0)

    L[1,1] = abs(γ⁰[1]) - l⁰
    @views L[band(0)][2:end] .=  (-).(b./s) .- l⁰
    @views L[band(-1)][1:end-1] .=  c[2:end].*γ¹ .- s.*a[1:end-1] .- l¹
    view(L,band(-1))[end] .= c∞*γ¹∞ - s∞*a[end] - l¹
    @views L[band(-2)][1:end-1] .= (-).(s[2:end].*b[1:end-1]) .- l²

    if n ≠ 1
        view(L,band(-2))[end] .= -s∞*b[end] - l²
    end

    Q,ToeplitzOperator([l¹,l²],[l⁰])+FiniteOperator(L,ℓ⁰,ℓ⁰)
end



discrete_eigs(J::SymTriPertToeplitz) =
    2*J.b*discrete_eigs(0.5*(J.dv-J.a)/J.b,0.5*J.ev/J.b) + J.a

connection_coeffs_operator(J::SymTriPertToeplitz) =
    connection_coeffs_operator(0.5*(J.dv-J.a)/J.b,0.5*J.ev/J.b)

# Spectral map takes SequenceSpace to the space of the diagonal Fun
struct SpectralMap{CC,QQ,RS,T} <: Operator{T}
    n::Int  # number of extra eigenvalues
    C::CC
    Q::QQ
    rangespace::RS
end

SpectralMap(n::Int,C::Operator{T},Q::Operator{T},rs) where {T} =
    SpectralMap{typeof(C),typeof(Q),typeof(rs),T}(n,C,Q,rs)


domainspace(::SpectralMap) = SequenceSpace()
rangespace(S::SpectralMap) = S.rangespace


function A_ldiv_B_coefficients(S::SpectralMap,v::AbstractVector;opts...)
    # leave first entries n alone
    r = copy(v)
    r[S.n+1:end] .= A_ldiv_B_coefficients(S.C,v[S.n+1:end])
    A_ldiv_B_coefficients(S.Q,r;opts...)
end

function A_mul_B_coefficients(S::SpectralMap,v::AbstractVector;opts...)
    r = A_mul_B_coefficients(S.Q,v;opts...)
    # leave first entries n alone
    r[S.n+1:end] .= A_mul_B_coefficients(S.C,r[S.n+1:end])
    r
end

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

# InvSpectralMap is the inverse of a Spectral map
struct InvSpectralMap{CC,QQ,DS,T} <: Operator{T}
    n::Int  # number of extra eigenvalues
    C::CC
    Q::QQ
    domainspace::DS
end

InvSpectralMap(n::Int,C::Operator{T},Q::Operator{T},ds) where {T} =
    InvSpectralMap{typeof(C),typeof(Q),typeof(ds),T}(n,C,Q,ds)


domainspace(S::InvSpectralMap) = S.domainspace
rangespace(::InvSpectralMap) = SequenceSpace()

function A_mul_B_coefficients(S::InvSpectralMap,v::AbstractVector;opts...)
    # leave first entries n alone
    r = copy(v)
    r[S.n+1:end] .= A_ldiv_B_coefficients(S.C,v[S.n+1:end])
    A_ldiv_B_coefficients(S.Q,r;opts...)
end

function A_ldiv_B_coefficients(S::InvSpectralMap,v::AbstractVector;opts...)
    r = A_mul_B_coefficients(S.Q,v;opts...)
    # leave first entries n alone
    r[S.n+1:end] .= A_mul_B_coefficients(S.C,r[S.n+1:end])
    r
end

function getindex(S::InvSpectralMap,k::Integer,j::Integer)
    v = A_mul_B_coefficients(S,[zeros(j-1);1])
    k ≤ length(v) && return v[k]
    zero(eltype(S))
end

isbanded(S::InvSpectralMap) = false
bandinds(S::InvSpectralMap) = [-∞;∞]

Base.inv(S::SpectralMap) = InvSpectralMap(S.n,S.C,S.Q,S.rangespace)
Base.inv(S::InvSpectralMap) = SpectralMap(S.n,S.C,S.Q,S.domainspace)

Base.eig(Jin::SymTriPertToeplitz) = eigfromguess(Jin,discrete_eigs(Jin))
function eigfromguess(Jin::SymTriPertToeplitz,approxeigs)
    Qret=Array{HessenbergUnitary{'U',Float64}}(0)
    λapprox=sort(approxeigs)

    ctsspec = ApproxFun.Segment(Jin.a-2Jin.b,Jin.a+2Jin.b)

    J=Jin

    if length(λapprox) == 0
        C=connection_coeffs_operator(J)

        D=Multiplication(Fun(identity,Ultraspherical(1,ctsspec)),Ultraspherical(1,ctsspec))

        U=SpaceOperator(C,ℓ⁰,domainspace(D))
        return D,U
    end

    λ=Array{Float64}(0)

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

         D=Multiplication(Fun(identity,PointSpace(λ)⊕Ultraspherical(1,ctsspec)),PointSpace(λ)⊕Ultraspherical(1,ctsspec))

         U=SpectralMap(length(λ),connection_coeffs_operator(J),Q,domainspace(D))
         return D,U
    else
        Q=BandedUnitary(reverse!(Qret))
        D=Multiplication(Fun(identity,PointSpace(λ)⊕Ultraspherical(1,ctsspec)),PointSpace(λ)⊕Ultraspherical(1,ctsspec))

        U=SpectralMap(length(λ),connection_coeffs_operator(J),Q,domainspace(D))
        return D,U
    end
end

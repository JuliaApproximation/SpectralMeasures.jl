
function ql(A::Matrix)
    Q,R=qr(A[end:-1:1,end:-1:1])
    Q[end:-1:1,end:-1:1],R[end:-1:1,end:-1:1]
end

joukowsky(z)=.5*(z+1./z)

function jacobimatrix(a,b,t0,t1,N)
    J = zeros(N,N)
    a = [a;t0*ones(N-length(a))]
    b = [b;t1*ones(N-length(b))]
    J[1,1]=a[1]
    for i = 1:N-1
        J[i+1,i+1] = a[i+1]
        J[i,i+1] = b[i]
        J[i+1,i] = b[i]
    end
    J
end

jacobimatrix(a,b,N) = jacobimatrix(a,b,0,.5,N)

function jacobioperator(a,b,t0,t1)
    n = max(length(a),length(b)+1)
    a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
    SymTriToeplitz(ToeplitzOperator([t1],[t0,t1]),SymTriOperator(a-t0,b-t1))
end

jacobioperator(a,b) = jacobioperator(a,b,0,.5)


DiscreteLaplacian()=ToeplitzOperator([0.5],[0.,0.5])


## tridiagonal ql

function tridql!(L::Matrix)
    n=size(L,1)

  # Now we do QL for the compact part in the top left
    cc=Array(eltype(L),n)
    ss=Array(eltype(L),n-1)

    for i = n:-1:2
        nrm=sqrt(L[i-1,i]^2+L[i,i]^2)
        c,s = L[i,i]/nrm, -L[i-1,i]/nrm
        if i > 2
            L[i-1:i,i-2:i] = [c s; -s c]*L[i-1:i,i-2:i]
            L[i-1,i]=0
        else
            L[i-1:i,i-1:i] = [c s; -s c]*L[i-1:i,i-1:i]
            L[i-1,i]=0
        end
        cc[i]=c
        ss[i-1]=s
    end
    cc[1]=sign(L[1,1])
    L[1,1]=abs(L[1,1])
    cc,ss,L
end


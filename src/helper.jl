
function ql(A::Matrix)
    Q,R=qr(A[end:-1:1,end:-1:1])
    Q[end:-1:1,end:-1:1],R[end:-1:1,end:-1:1]
end

joukowsky(z)=.5*(z+1./z)

function jacobimatrix(a,b,t0,t1,N)
    J = BandedMatrix(Float64,N,N,1,1)

    for i = 1:min(length(a),N)
        J[i,i] = a[i]
    end
    for i=length(a)+1:N
        J[i,i] = t0
    end

    for i = 1:min(length(b),N-1)
        J[i,i+1] = J[i+1,i] = b[i]
    end

    for i=length(b)+1:N-1
        J[i,i+1] = J[i+1,i] = t1
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


function tridql!(J::BandedMatrix)
    n=size(J,1)
    L=BandedMatrix(copy(J.data),J.m,2,0)

  # Now we do QL for the compact part in the top left
    cc=Array(eltype(J),n)
    ss=Array(eltype(J),n-1)

    for i = n:-1:2
        nrm=sqrt(J[i-1,i]^2+J[i,i]^2)
        c,s = J[i,i]/nrm, -J[i-1,i]/nrm

        for j=max(i-2,1):i
            L[i,j]=-s*J[i-1,j]+c*J[i,j]
            J[i-1,j]=c*J[i-1,j]+s*J[i,j]
        end
        cc[i]=c
        ss[i-1]=s
    end
    cc[1]=sign(J[1,1])
    L[1,1]=abs(J[1,1])
    cc,ss,L
end

#Finds NxN truncation of C such that C'(Q_k(s)) =  (P_k(s)),
# where P_k has Jacobi coeffs a,b and Q_k has Jacobi coeffs c,d
function connectionCoeffsMatrix(a,b,c,d,N)
  if N>max(length(a),length(b)+1,length(c),length(d)+1)
    a = [a;zeros(N-length(a))]; b = [b;.5+zeros(N-length(b))]
    c = [c;zeros(N-length(c))]; d = [d;.5+zeros(N-length(d))]
  end

  C = zeros(N,N)
  C[1,1] = 1
  C[1,2] = (c[1]-a[1])/b[1]
  C[2,2] = d[1]/b[1]
  for j = 3:N
    C[1,j] = ((c[1]-a[j-1])*C[1,j-1] + d[1]*C[2,j-1] - b[j-2]*C[1,j-2])/b[j-1]
    for i = 2:j-1
      C[i,j] = (d[i-1]*C[i-1,j-1] + (c[i]-a[j-1])*C[i,j-1] + d[i]*C[i+1,j-1] - b[j-2]*C[i,j-2])/b[j-1]
    end
    C[j,j] = d[j-1]*C[j-1,j-1]/b[j-1]
  end
  C
end

# This is for Chebyshev U
connectionCoeffsMatrix(a,b,N) = connectionCoeffsMatrix(a,b,[],[],N)

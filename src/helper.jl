
function ql(A::Matrix)
    Q,R=qr(A[end:-1:1,end:-1:1])
    Q[end:-1:1,end:-1:1],R[end:-1:1,end:-1:1]
end

joukowsky(z)=.5*(z+1./z)

function onesAndZeros(n)
  v = ones(n)
  for i = 1:div(n,2)
    v[2*i] = 0
  end
  v
end


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

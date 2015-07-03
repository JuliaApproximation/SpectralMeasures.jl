
function ql(A::Matrix)
    Q,R=qr(A[end:-1:1,end:-1:1])
    Q[end:-1:1,end:-1:1],R[end:-1:1,end:-1:1]
end
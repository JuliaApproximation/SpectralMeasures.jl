Δ=ToeplitzOperator([1.],[0.,1.])
@manipulate for n=2:10   
    K=SymTriOperator(rand(n),rand(n-1))
    eigvals(Δ+K)|>ApproxFun.plot
end


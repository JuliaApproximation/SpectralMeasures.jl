immutable RatFun{S1,T1,S2,T2}
  p::Fun{S1,T1}
  q::Fun{S2,T2}
end

evaluate(r::RatFun,x)=p(x)./q(x)


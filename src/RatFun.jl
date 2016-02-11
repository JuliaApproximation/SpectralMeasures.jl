immutable RatFun{S,T}
  p::Fun{S}
  q::Fun{T}
end

evaluate(r::RatFun,x)=p(x)./q(x)


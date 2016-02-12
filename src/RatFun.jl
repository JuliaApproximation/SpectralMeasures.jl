using ApproxFun

immutable RatFun{S1,T1,S2,T2}
  p::Fun{S1,T1}
  q::Fun{S2,T2}
end

function ratEval(r::RatFun,x)
  (r.p)(x)./(r.q)(x)
end

function ratPlot(r::RatFun)
  x = -1:0.001:1; y = ratEval(r,x);
  ApproxFun.plot(x, y, linewidth=2.0)
end

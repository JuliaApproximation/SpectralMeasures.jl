using ComplexPhasePortrait, ImageMagick


function plotdiscresolvent(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  C = SpectralMeasures.connectionCoeffsOperator(a,b)
  Cmu = SpectralMeasures.connectionCoeffsOperator(a[2:end],b[2:end])
  c = Fun([C.T[1,1];C.T.negative],Taylor)
  cmu = Fun([0;Cmu.T[1,1];Cmu.T.negative]/b[1],Taylor)
  nx = 1000
  x = linspace(-2, 2, nx)
  Z = x' .+ flipdim(x, 1)*im
  portrait(-cmu(Z)./c(Z),PTstepmod)
end

function plotsplitplaneresolvent(a,b)
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]
  C = SpectralMeasures.connectionCoeffsOperator(a,b)
  Cmu = SpectralMeasures.connectionCoeffsOperator(a[2:end],b[2:end])
  f = Fun(C'*(C*[1]),Ultraspherical{1}())
  fmu = Fun(Cmu'*((C*[1])[2:end])/b[1],Ultraspherical{1}())
  nx = 1000
  x = linspace(-2, 2, nx)
  Z = x' .+ flipdim(x, 1)*im
  portrait((2*sqrt(Z-1).*sqrt(Z+1)-2*Z-fmu(Z))./f(Z),PTstepmod)
end
using ApproxFun, SpectralMeasures

# randn and 1/n^2 decay gives good error bounds
n = 50
a = randn(n)./(1:n).^2
b = max(.05,.5+randn(n-1)./(1:n-1).^2)
spectrum, maxerr = validated_spectrum(a,b)
maxerr

# randn and 1/n decay gives okay error bounds
n = 50
a = randn(n)./(1:n)
b = max(.05,.5+randn(n-1)./(1:n-1))
spectrum, maxerr = validated_spectrum(a,b)
maxerr

# randn with no decay gives terrible error bounds
n = 20
a = randn(n)
b = max(.05,.5+randn(n-1))
spectrum, maxerr = validated_spectrum(a,b)
maxerr

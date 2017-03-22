using SpectralMeasures, ValidatedNumerics, ValidatedNumerics.RootFinding,
            ComplexPhasePortrait


function tripleplot(a,b,Z=linspace(-3, 3, 300).+linspace(3,-3,300)'*im)
  # Build the measure and resolvents
  μ = spectralmeasure(a,b)
  R = principal_resolvent(a,b)
  r = disc_resolvent(a,b)

  # Create the subplots
  p1 = plot(μ,xlims=(-2,2),ylims=(0,1.5))
  p2 = plot(portrait(R(Z),PTstepmod),xlims=(-3,3),ylims=(-3,3),aspect_ratio=1)
  plot!([-1.,1.],[0.,0.],linewidth=3,color=:black)
  p3 = plot(portrait(r(Z),PTstepmod),xlims=(-3,3),ylims=(-3,3),aspect_ratio=1)
  plot!(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=3,color=:black)

  # Plot the subplots
  l = @layout [a{0.7w}; b c]
  plot(p1,p2,p3,layout=l)
end

function validated_spectrum(a,b)
  # Chop the a and b down
  a = chop!(a); b = .5+chop!(b-.5)
  n = max(2,length(a),length(b)+1)
  a = [a;zeros(n-length(a))]; b = [b;.5+zeros(n-length(b))]

  # Finds C such that J*C = C*Toeplitz([0,1/2])
  C = SpectralMeasures.connection_coeffs_operator(a,b)
  c = Fun(Taylor,C.T.nonnegative)

  rts = find_roots(x->c(x),-1,1)
  if length(rts) > 0
     eigs=real(map(x->joukowsky(x.interval),rts))
     eigserrs = map(midpoint_radius,eigs)
     spectrum = ApproxFun.Interval(-1,1)
     maxerr = 0
     for (eig,err) in eigserrs
       maxerr = max(maxerr,err)
       spectrum = ApproxFun.Point(eig) ∪ spectrum
     end
  else
    spectrum = ApproxFun.Interval(-1,1)
    maxerr = 0
  end
  spectrum, maxerr
end

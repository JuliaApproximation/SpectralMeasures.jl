## These are the figures that appear in a first draft of the paper
# "Spectra of Jacobi operators via connection coefficients" by Marcus Webb and Sheehan Olver (2016)

using ApproxFun, SpectralMeasures, Plots; Plots.gr(linewidth=3,legend=false)
using ComplexPhasePortrait, PyPlot, FileIO; pygui(false) # PyPlot is good for Wegert Plots

x = linspace(-3, 3, 1000)
  Z = x' .+ flipdim(x, 1)*im

mkdir("Figs")

# Example for nondecreasingness of resolvents in intervals of continuity
MyResolvent = x -> real(1 + 2*sqrt(complex(x-1.)).*sqrt(complex(x+1.))-2*x + x./(9-x.^2))
  Plots.plot([-1,1],[0,0],color=:black)
  plot!([-3,-3,-3],[-3,0,3],line=:dash,color=:black,marker=:circle,markersize=2,linewidth=1)
  plot!([3,3,3],[-3,0,3],line=:dash,color=:black,marker=:circle,markersize=2,linewidth=1)
  plot!([-1,-1],[-3,3],line=:dash,color=:black,linewidth=1)
  plot!([1,1],[-3,3],line=:dash,color=:black,linewidth=1)
  plot!(-8:0.001:8,[MyResolvent(-8:0.001:-1);NaN*(-1.001:0.001:0.999);MyResolvent(1:0.001:8)],xlabel="\\lambda",ylabel="f\\(\\lambda\\)",xlims=(-7.5,7.5),ylims=(-2,2))
  Plots.savefig("Figs/secular.pdf")

# Pertrubation of the [1,1] element
for a in [0, 0.15, 0.35, 0.5, 0.75, 1.]
  Plots.plot(spectralmeasure([a],[.5]),xlims=(-2,2),ylims=(0,2),title="\\alpha=$(a)")
  Plots.savefig("Figs/spectralmeasurealpha=$(round(Int,100*a)).pdf")

  R = principal_resolvent([a],[.5])
  save("Figs/principalresolventalpha=$(round(Int,100*a)).png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolventalpha=$(round(Int,100*a)).png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolventalpha=$(round(Int,100*a)).png",transparent=true,bbox_inches="tight")
  clf()

  r = disc_resolvent([a],[.5])
  save("Figs/discresolventalpha=$(round(Int,100*a)).png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolventalpha=$(round(Int,100*a)).png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolventalpha=$(round(Int,100*a)).png",transparent=true,bbox_inches="tight")
  clf()
end

# Perturbation of the [1,2] element
for b in [.5, 0.707, .85, 1, 1.2, 1.5]
  Plots.plot(spectralmeasure([0.],[b]/sqrt(2)),xlims=(-2,2),ylims=(0,2),title="\\beta=$(b)")
  Plots.savefig("Figs/spectralmeasurebeta=$(round(Int,1000*b)).pdf")

  R = principal_resolvent([0.],[b]/sqrt(2))
  save("Figs/principalresolventbeta=$(round(Int,1000*b)).png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolventbeta=$(round(Int,1000*b)).png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolventbeta=$(round(Int,1000*b)).png",transparent=true,bbox_inches="tight")
  clf()

  r = disc_resolvent([0.],[b]/sqrt(2))
  save("Figs/discresolventbeta=$(round(Int,1000*b)).png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolventbeta=$(round(Int,1000*b)).png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolventbeta=$(round(Int,1000*b)).png",transparent=true,bbox_inches="tight")
  clf()
end



# 3 by 3 Perturbation
a = [3/4,-1/4,1/2]; b=[1,3/4]
  μ = spectralmeasure(a,b)
  R = principal_resolvent(a,b)
  r = disc_resolvent(a,b)

  Plots.plot(μ,xlims=(-3,3),ylims=(0,2))
  Plots.savefig("Figs/spectralmeasure3by3pert.pdf")

  save("Figs/principalresolvent3by3pert.png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolvent3by3pert.png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolvent3by3pert.png",transparent=true,bbox_inches="tight")
  clf()

  save("Figs/discresolvent3by3pert.png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolvent3by3pert.png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolvent3by3pert.png",transparent=true,bbox_inches="tight")
  clf()


# Legendre
for n in [1,2,3,10,30,100]
  bLeg = (1:n-1)./sqrt(4*(1:n-1).^2-1)
  Plots.plot(spectralmeasure(zeros(n),bLeg),title="Legendre approximation, n = $(n)",xlims=(-1.5,1.5),ylims=(0,1),yticks=[0,.25,.5,.75,1])
  Plots.savefig("Figs/spectralmeasureLegn=$(n).pdf")

  R = principal_resolvent([0.],bLeg)
  save("Figs/principalresolventLegn=$(n).png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolventLegn=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolventLegn=$(n).png",transparent=true,bbox_inches="tight")
  clf()

  r = disc_resolvent([0.],bLeg)
  save("Figs/discresolventLegn=$(n).png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolventLegn=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolventLegn=$(n).png",transparent=true,bbox_inches="tight")
  clf()
end

# Ultraspherical polynomials
γ = 0.6
  for n in [1,2,3,10,30,100]
  bUlt = .5*sqrt(((1:n).*(2γ+(0:n-1)))./((γ+(0:n-1)).*(γ+(1:n))))
  Plots.plot(spectralmeasure(zeros(n),bUlt),title="Ultraspherical\\($(γ)\\), n = $(n)",xlims=(-1.5,1.5),ylims=(0,1),yticks=[0,.25,.5,.75,1])
  Plots.savefig("Figs/spectralmeasureUltra0pt6n=$(n).pdf")

  R = principal_resolvent([0.],bUlt)
  save("Figs/principalresolventUltra0pt6n=$(n).png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolventUltra0pt6n=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolventUltra0pt6n=$(n).png",transparent=true,bbox_inches="tight")
  clf()

  r = disc_resolvent([0.],bUlt)
  save("Figs/discresolventUltra0pt6n=$(n).png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolventUltra0pt6n=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolventUltra0pt6n=$(n).png",transparent=true,bbox_inches="tight")
  clf()
end

# Jacobi polynomials
α = 0.4; β = 1.9
  for n in [1,2,3,10,30,100]
  aJac = (β.^2-α.^2)./((2.*(0:n)+α+β).*(2.*(1:n+1)+α+β))
  bJac = 2*sqrt(((1:n).*(α+(1:n)).*(β+(1:n)).*(α+β+(1:n)))./((2.*(1:n)+α+β-1).*((2.*(1:n)+α+β).^2).*(2.*(1:n)+α+β+1)))
  Plots.plot(spectralmeasure(aJac,bJac),title="Jacobi\\($(α),$(β)\\), n = $(n)",xlims=(-1.5,1.5),ylims=(0,1),yticks=[0,.25,.5,.75,1])
  Plots.savefig("Figs/spectralmeasureJacobin=$(n).pdf")

  R = principal_resolvent(aJac,bJac)
  save("Figs/principalresolventJacobin=$(n).png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolventJacobin=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolventJacobin=$(n).png",transparent=true,bbox_inches="tight")
  clf()

  r = disc_resolvent(aJac,bJac)
  save("Figs/discresolventJacobin=$(n).png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolventJacobin=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolventJacobin=$(n).png",transparent=true,bbox_inches="tight")
  clf()
end




# Random
mkdir("Figs")
for n in [1,2,3,10,50,100]
  srand(200)
  aRand = 3*(2rand(n)-1)./(1:n).^2
  Plots.plot(spectralmeasure(aRand,Float64[]),title="Random, n = $(n)",xlims=(-1.5,1.5),ylims=(0,1),yticks=[0,.25,.5,.75,1])
  Plots.savefig("Figs/spectralmeasurerandomn=$(n).pdf")

  R = principal_resolvent(aRand,Float64[])
  save("Figs/principalresolventrandomn=$(n).png", portrait(R(Z),PTstepmod))
  imshow(imread("Figs/principalresolventrandomn=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot([-1.,1.],[0.,0.],linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/principalresolventrandomn=$(n).png",transparent=true,bbox_inches="tight")
  clf()

  r = disc_resolvent(aRand,Float64[])
  save("Figs/discresolventrandomn=$(n).png", portrait(r(Z),PTstepmod))
  imshow(imread("Figs/discresolventrandomn=$(n).png"),extent=[-3,3,-3,3])
  PyPlot.plot(cos(linspace(0,2pi,100)),sin(linspace(0,2pi,100)),linewidth=2,color=:black)
  PyPlot.grid("on")
  PyPlot.savefig("Figs/discresolventrandomn=$(n).png",transparent=true,bbox_inches="tight")
  clf()
end

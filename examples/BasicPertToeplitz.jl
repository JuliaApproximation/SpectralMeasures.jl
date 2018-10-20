using ApproxFun, SpectralMeasures, Plots; Plots.gr(legend=false,linewidth=2,xlims=(-2,2),ylims=(0,2))

## 1by1 perturbation of Toeplitz
# J1 = [α/2 .5            ]
#      [ .5  0 .5         ]
#      [    .5  0 .5      ]
#      [       .5  0 .5   ]
#      [          ....... ]
#
#  We have several points of interest
#   α = 0 : Chebyshev of the 2nd Kind (semicircle)
#   α = -1 : Chebyshev of the 3rd kind
#   α = 1 : Chebyshev of the 4th kind
#  |α|≦ 1 : purely continuous spectrum
#  |α|> 1 : A single eigenvalue
#
#  In LightTable, "drag" k to see how the perturbation affects the spectral measure

k = 3
  plot(spectralmeasure([k/20],[.5]),title="\\alpha = $(k/10)")


## Basic 2by2 perturbation of Toeplitz
# J2 = [  0  β/⎷2            ]
#      [β/⎷2  0  .5          ]
#      [     .5   0  .5      ]
#      [         .5   0 .5   ]
#      [            ........ ]
#
#  We have several points of interest
#   β = 1/⎷2 : Chebyshev of the 2nd Kind
#   β = 1 : Chebyshev of the 1st kind
#  |β|≦ 1 : purely continuous spectrum
#  |β|> 1 : 2 eigenvalues
#
#  In LightTable, "drag" k to see how the perturbation affects the spectral measure

k = 5
  plot(spectralmeasure([0.],[.5+k/20]/sqrt(2)),title="\\beta=$(.5+k/20)")


####################
## Use the following for a closer look at what is going on with these operators

Δ = freejacobioperator()

a = 1.0
J1 = SymTriPertToeplitz([a],[.5],0.0,0.5)
C1 = connectioncoeffsoperator(J1)
Δ*C1-C1*J1
10.0*I
Q,L = ql(J1-10.0*I)

d,U = eigen(J1)



b = .5 + 0.1
J2 = SymTriPertToeplitz([0.],[b],0.0,0.5)

C2 = connectioncoeffsoperator(J2)
Δ*C2-C2*J2

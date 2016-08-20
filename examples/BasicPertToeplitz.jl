using ApproxFun
using SpectralMeasures
using Plots; Plots.gr(legend=false,linewidth=2,xlims=(-2,2),ylims=(0,2))

## 1by1 perturbation of Toeplitz
# J1 = [a/2 .5            ]
#      [ .5  0 .5         ]
#      [    .5  0 .5      ]
#      [       .5  0 .5   ]
#      [          ....... ]
#
#  We have several points of interest
#   a = 0 : Chebyshev of the 2nd Kind (semicircle)
#   a = -1 : Chebyshev of the 3rd kind
#   a = 1 : Chebyshev of the 4th kind
#  |a|≦ 1 : purely continuous spectrum
#  |a|> 1 : A single eigenvalue
#
#  In LightTable, "drag" k to see how the perturbation affects the spectral measure

k = 0
  plot(spectralmeasure([k/20],[.5]))


## Basic 2by2 perturbation of Toeplitz
# J2 = [  0  b/⎷2            ]
#      [b/⎷2  0  .5          ]
#      [     .5   0  .5      ]
#      [         .5   0 .5   ]
#      [            ........ ]
#
#  We have several points of interest
#   b = 1/⎷2 : Chebyshev of the 2nd Kind
#   b = 1 : Chebyshev of the 1st kind
#  |b|≦ 1 : purely continuous spectrum
#  |b|> 1 : 2 eigenvalues
#
#  In LightTable, "drag" k to see how the perturbation affects the spectral measure


k = 4
  plot(spectralmeasure([0.],[.5+k/20]/sqrt(2)))


####################
## Use the following for a closer look at what is going on with these operators

Δ = DiscreteLaplacian()

a = 1.0
J1 = SymTriToeplitz([a],[.5],0.0,0.5)
C1 = connectionCoeffsOperator(J1)
Δ*C1-C1*J1
10.0*I
Q,L = ql(J1-10.0*I)

d,U = eig(J1)



b = .5 + 0.1
J2 = SymTriToeplitz([0.],[b],0.0,0.5)

C2 = connectionCoeffsOperator(J2)
Δ*C2-C2*J2

# SpectralMeasures.jl

SpectralMeasures.jl is a Julia package for computing spectra and spectral measures of self-adjoint operators with structured, infinite-dimensional matrix representations. The package enables a canonical eigendecomposition for these operators, and a functional calculus too. A key point to make is that all of this is done in tailor-made data structures representing structured infinite-dimensional operators, with no error-inducing, finite-dimensional truncations.

## Overview

Syntactically, the familiar function for matrices, `eig`, will be implemented for as many operators as we can. It returns a multiplication operator, `D`, which acts on `Fun`s whose domain is the spectrum of the operator, and "spectral map", `U`, which takes vectors to `Fun`s in the same `Space` as the that associated with `D`. In Julia code, this means that if you have an operator `A`, and type

```
D,U = eig(A)
```

then the following identity holds (to machine precision):

```
A = U \ D * U
```

This is analogous to what happens in finite dimensional linear algebra, where `D` is a diagonal matrix and `U` is the matrix of eigenvectors (transposed). You can also try `eigvals(A)` or `spectrum(A)` to return the spectrum of `A` as an ApproxFun `Domain` type.


The operators `D` and `U` enable a functional calculus. The following is the application of the operator `exp(1.3*A)` to a random vector.
```
v = randn(100)
expAtv = U \ exp(1.3*D) *  U * v
```

Also, specific to infinite-dimensional operators, there is a function called `spectralmeasure`. This function produces either a `Fun` or a `RatFun` (a quotient of two `Fun`s), representing a measure supported on the spectrum of the operator --- specifically, a representation of the spectral measure from the Spectral Theorem for selfadjoint operators. At a discrete eigenvalue of the operator, the spectral measure has a Dirac delta.

The implementation of the package is based on a combination of connection coefficient matrices between orthogonal polynomials (See Webb-Olver 2018, "Spectra of Jacobi operators via connection coefficients matrices") and infinite dimensional QL iterations (See Olver-Webb 2018, "The infinite dimensional QL algorithm"  (in prep.)).

## Tridiagonal operators

### PertToeplitz

The type `SymTriPertToeplitz` implements an operator whose matrix is symmetric, tridiagonal, which is a finite-rank perturbation of a Toeplitz operator. The following is an example construction.
```
using SpectralMeasures
A = SymTriPertToeplitz([0.1;-0.1;1.2],[1.4;1.2],0.,0.5)
```
The data structure represents the infinite matrix,
```
A = [ 0.1  1.4        
      1.4 -0.1  1.2
           1.2  1.2  0.5
                0.5  0.0  0.5
                     0.5  0.0  ⋱
                           ⋱  ⋱ ]
```
Simple commands like `A*randn(100)` and `A[2017,2018]` work as you would expect.

#### Eigendecomposition and functional calculus

The spectrum of `A` can be computed using the `eigvals` function. `spectrum` can also be used.

```
eigvals(A)
```
returns a `Domain` type,
```
Point(-1.7082248431502554)∪Point(2.2839756035160352)∪【-1.0,1.0】
```
The spectrum of this particular operator consists of two discrete eigenvalues and continuous spectrum on the interval [-1,1].

You can use familiar syntax to compute an eigendecomposition:

```
D,U = eig(A)
```

The first output, `D`, is a multiplication operator, which multiplies `Fun`s with space `PointSpace([-1.70822…,2.28398…])∪Ultraspherical(1)`. It has matrix representation,
```
D = [ -1.70822  0.0      0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  ⋯
       0.0      2.28398  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  ⋱
       0.0      0.0      0.0  0.5  0.0  0.0  0.0  0.0  0.0  0.0  ⋱
       0.0      0.0      0.5  0.0  0.5  0.0  0.0  0.0  0.0  0.0  ⋱
       0.0      0.0      0.0  0.5  0.0  0.5  0.0  0.0  0.0  0.0  ⋱
       0.0      0.0      0.0  0.0  0.5  0.0  0.5  0.0  0.0  0.0  ⋱
       0.0      0.0      0.0  0.0  0.0  0.5  0.0  0.5  0.0  0.0  ⋱
       0.0      0.0      0.0  0.0  0.0  0.0  0.5  0.0  0.5  0.0  ⋱
       0.0      0.0      0.0  0.0  0.0  0.0  0.0  0.5  0.0  0.5  ⋱
       0.0      0.0      0.0  0.0  0.0  0.0  0.0  0.0  0.5  0.0  ⋱
        ⋮        ⋱        ⋱    ⋱   ⋱   ⋱    ⋱   ⋱   ⋱    ⋱  ⋱ ]
```
 The second output, `U` is a an operator which maps vectors to `Fun`s in that space. It has matrix representation,

```
U = [ -0.575073   0.742758   -0.324516    0.104914   -0.0339182   0.0109656   -0.00354511   0.00114611   -0.000370532   0.000119791  ⋱
       0.368164   0.57433     0.711466    0.164031    0.0378177   0.00871899   0.00201019   0.000463455   0.000106851   2.46348e-5   ⋱
      -1.07598    0.160981    1.1016     -3.23631     1.29719     0.217432    -0.0841252    9.59302e-16   4.7358e-16    2.61726e-16  ⋱
       0.235551  -0.401105   -0.206107    2.39879    -3.01888     1.21307      0.217432    -0.0841252     1.21084e-15   5.3256e-16   ⋱
                  0.0841252  -0.160117    0.0113248   2.31467    -3.01888      1.21307      0.217432     -0.0841252     1.23512e-15  ⋱
                              0.0350522  -0.244242    0.0113248   2.31467     -3.01888      1.21307       0.217432     -0.0841252    ⋱
                                          0.0350522  -0.244242    0.0113248    2.31467     -3.01888       1.21307       0.217432     ⋱
                                                      0.0350522  -0.244242     0.0113248    2.31467      -3.01888       1.21307      ⋱
                                                                  0.0350522   -0.244242     0.0113248     2.31467      -3.01888      ⋱
                                                                               0.0350522   -0.244242      0.0113248     2.31467      ⋱
                                                                                             ⋱             ⋱            ⋱          ⋱ ]
```

 The operator `A - U \ D * U` has entries which are approximately zero to machine precision.

Functions of `A` applied to a vector can be computed as follows.

```
f = x -> exp(1.4im*x)
v = randn(100)
fAv = U\ f(D) * U * v
```
The key point to note is that no finite-dimensional truncation has occurred; the entire structured operators have been used at each stage of the computation.

#### Spectral measure and resolvent

We can compute the spectral measure of `A`:
```
μ = spectralmeasure(A)
```
μ is a RationalFun from the RatFun package, and has Dirac deltas at the discrete eigenvalues of A. Let us plot this.
```
using Plots
plot(μ)
```
<img src=images/spectralmeasure.png width=400px>

We can also plot the resolvent of the operator (which is the Cauchy transform of the spectral measure):
```
r = principalresolvent(A)
using ComplexPhasePortrait
Z=linspace(-3, 3, 600)'.+linspace(3,-3,600)*im
plot(portrait(r(Z),PTstepmod),xlims=(-3,3),ylims=(-3,3),aspect_ratio=1)
```
<img src=images/principalresolvent.png width=400px>

This plot is a Wegert plot, a.k.a. a complex phase portrait. See <a href=https://github.com/JuliaHolomorphic/ComplexPhasePortrait.jl>here</a> for an explanation.

#### QL iterations

We can compute a QL decomposition in the expected way. However, note that QL decompositions do not exist if the continuous spectrum contains zero (see Olver-Webb 2018, "The infinite dimensional QL algorithm"). Therefore, we must use a shift.
```
Q,L = ql(A+1.7082248431502554)
```

The operators `Q` and `L` are stored in their entirety as structured operators of type `HessenbergUnitary` and `PertToeplitz`, respectively.
```
Q = [ -0.575073      0.818102
       0.742758      0.52211       0.41918
      -0.324516     -0.228114      0.859147      0.323294
       0.104914      0.0737479    -0.277758      0.895481     0.323294
      -0.0339182    -0.0238423     0.0897975    -0.289504     0.895481     0.323294
       0.0109656     0.00770808   -0.029031      0.093595    -0.289504     0.895481     0.323294
      -0.00354511   -0.00249198    0.00938557   -0.0302588    0.093595    -0.289504     0.895481    0.323294
       0.00114611    0.000805644  -0.0030343     0.00978249  -0.0302588    0.093595    -0.289504    0.895481   0.323294
      -0.000370532  -0.00026046    0.000980974  -0.00316262   0.00978249  -0.0302588    0.093595   -0.289504   0.895481  0.323294
       0.000119791   8.42053e-5   -0.000317143   0.00102246  -0.00316262   0.00978249  -0.0302588   0.093595  -0.289504  0.895481  ⋱
        ⋱             ⋱            ⋱             ⋱           ⋱            ⋱           ⋱           ⋱         ⋱       ⋱        ⋱ ]

L = [ 2.22045e-16
      2.21027      1.71128
      0.586852     1.70511   2.86273
                   0.387953  1.38795   1.54658
                             0.161647  1.0       1.54658
                                       0.161647  1.0       1.54658
                                                 0.161647  1.0       1.54658
                                                           0.161647  1.0       1.54658
                                                                     0.161647  1.0       1.54658
                                                                               0.161647  1.0      1.54658
                                                                                          ⋱        ⋱       ⋱ ]
```

A QL iteration can be performed by forming the operator,
```
A2 = L*Q - 1.7082248431502554I
```
```
A2 = [ -1.70822      1.81655e-16
        1.81655e-16  0.993475     0.717334
                     0.717334     1.46603   0.925505
                                  0.925505  0.125423   0.5
                                            0.5       -2.22045e-16   0.5
                                                       0.5          -2.22045e-16   0.5
                                                                     0.5          -2.22045e-16   0.5
                                                                                   0.5          -2.22045e-16   0.5
                                                                                                 0.5          -2.22045e-16   0.5
                                                                                                               0.5          -2.22045e-16  ⋱
                                                                                                                              ⋱           ⋱ ]
```
Since the shift was approximately an eigenvalue of `A`, the (1,2) entry of the QL iterate becomes approximately zero. Repeated iteration will yield convergence just like in the finite dimensional case (see Olver-Webb 2018, "The infinite dimensional QL algorithm" (in prep.))

## Future operators
Tridiagonal bi-infinite PertToeplitz, Tridiagonal PertPeriodic. General banded operators. Work on all of these is currently in progress. Work on continuous Schrödinger operators and nonsymmetric operators is something which may be done in the distant future.

# SpectralMeasures.jl

SpectralMeasures.jl is a Julia package for computing spectra and spectral measures of self-adjoint operators with a structured, infinite-dimensional matrix representation. The package enables a canonical eigendecomposition for these operators, and whence a functional calculus too. A key point to make is that all of this is done in tailor-made data structures representing structured infinite-dimensional operators, with no error-inducing truncations.

## Overview

Syntactically, the familiar function for matrices, `eig`, will be implemented for as many operators as we can. It returns a multiplication operator, `D`, which acts on `Fun`s whose domain is the spectrum of the operator, and "spectral map", `U`, which takes vectors to `Fun`s in the same `Space` as the that associated with `D`. In Julia code, this means that if you have an operator `A`, and type

```
D,U = eig(A)
```

then the following identity holds (to machine precision):

```
A = U \ D * U
```

This is analogous to what happens in finite dimensional linear algebra, where `D` is a diagonal matrix and `U` is the matrix of eigenvectors (transposed). You can also try `eigvals(A)` to return the spectrum of `A` as an ApproxFun `Domain` type.


The operators `D` and `U` enable a functional calculus. The following is the application of the operator `exp(1.3*A)` to a random vector.
```
v = randn(100)
expAtv = U \ exp(1.3*D) *  U * v
```

Also, specific to infinite-dimensional operators, there is a function called `spectralmeasure`. This function produces either a `Fun` or a `RatFun` (a quotient of two `Fun`s), representing a measure supported on the spectrum of the operator --- specifically, a representation of the spectral measure from the Spectral Theorem for selfadjoint operators. At a discrete eigenvalue of the operator, the spectral measure has a Dirac delta.

The implementation of the package is based on a combination of connection coefficient matrices between orthogonal polynomials (See Webb-Olver 2018, "Spectra of Jacobi operators via connection coefficients matrices") and infinite dimensional QL iterations (See Olver-Webb 2018, "The infinite dimensional QL algorithm").

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

```
D,U = eig(A)
```

```
eigvals(A)
```
returns
```
Point(-1.7082248431502554)∪Point(2.2839756035160352)∪【-1.0,1.0】
```
There are two discrete eigenvalues and a single interval of continuous spectrum.

(TODO: add functional calculus example)

#### Spectral measure and resolvent

We can compute the spectral measure of `A`,
```
μ = spectralmeasure(A)
```
μ is a RationalFun from the RatFun package, and has Dirac deltas at the discrete eigenvalues of A. Let us plot this.
```
using Plots
plot(μ)
```
<img src=images/spectralmeasure.png>

We can also plot the resolvent of the operator (which is the Cauchy transform of the spectral measure):
```
r = principalresolvent(A)
using ComplexPhasePortrait
Z=linspace(-3, 3, 300).+linspace(3,-3,300)'*im
plot(portrait(r(Z),PTstepmod),xlims=(-3,3),ylims=(-3,3),aspect_ratio=1)
```
<img src=images/principalresolvent.png>

This plot is a Wegert plot, a.k.a. a complex phase portrait. See <a href=https://github.com/JuliaHolomorphic/ComplexPhasePortrait.jl>here</a>.

#### QL iterations

We can compute a QL decomposition in the expected way. However, note that QL decompositions do not exist if the continuous spectrum contains zero (see Olver-Webb 2018, The infinite dimensional QL algorithm). Therefore, we must use a shift.
```
Q,L = ql(A+1.7082248431502554)
```

(TODO: QL iterations)

## Future operators
Tridiagonal bi-infinite PertToeplitz, Tridiagonal PertPeriodic. General banded operators. Work on all of these is currently in progress.

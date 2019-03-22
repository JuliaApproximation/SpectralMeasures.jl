using InfiniteArrays, FillArrays, BlockArrays, BlockBandedMatrices

Δ = mortar(Tridiagonal(Vcat([ones(2,1)],Fill([1.0 0; 0.0 1.0],∞)), Vcat([zeros(1,1)],Fill(zeros(2,2),∞)), Vcat([ones(1,2)],Fill([1.0 0; 0.0 1.0],∞))))

Vcat([ones(2,1)],Fill([1.0 0; 0.0 1.0],∞))
Vcat(zeros(1,1),Fill(zeros(2,2),∞))
Vcat(ones(1,2),Fill([1.0 0; 0.0 1.0],∞))[
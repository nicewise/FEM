"""
'''julia
import Pardiso
'''
ps = Pardiso.MKLPardisoSolver:
	Matrix type: Real nonsymmetric
	Phase: Analysis, numerical factorization, solve, iterative refinement
"""
function MKL(A, b)
    ps = Pardiso.MKLPardisoSolver() # 默认Matrix type: Real nonsymmetric
    Pardiso.set_nprocs!(ps, 6)
    x = similar(b)
    Pardiso.solve!(ps, x, A, b)
    x
end

"""
'''julia
import SymRCM
'''
"""
function reorder_and_slash(A, b)
    d = SymRCM.symrcm(A)
    x = zeros(length(b))
    x[d] = A[d, d] \ b[d]
    x
end

"""
'''julia
import LinearSolve
'''
and other packages solver needs
"""
@inline function mysolver(A, b)
    prob = LinearSolve.LinearProblem(A, b)
    sol = LinearSolve.solve(prob, LinearSolve.MKLPardisoFactorize(), nprocs=8)
    sol.u
end

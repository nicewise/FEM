"""
'''julia
import Pardiso
'''
ps = Pardiso.MKLPardisoSolver:
	Matrix type: Real nonsymmetric
	Phase: Analysis, numerical factorization, solve, iterative refinement
"""
struct MKLSolver <: AbstractSolver
    ps
end
function MKLSolver(nprocs::Int)
    ps = Pardiso.MKLPardisoSolver() # 默认Matrix type: Real nonsymmetric
    Pardiso.set_nprocs!(ps, nprocs)
    MKLSolver(ps)
end
function solveit(solver::MKLSolver, A, b)
    x = similar(b)
    Pardiso.solve!(ps, x, A, b)
    return x
end

"""
'''julia
import LinearSolve
'''
and other packages solver needs
"""
struct LinearSolveDotjl <: AbstractSolver
    solver
end
function solveit(solver::LinearSolveDotjl, A, b)
    prob = LinearSolve.LinearProblem(A, b)
    sol = LinearSolve.solve(prob, solver.solver)
    return sol.u
end

"""
'''julia
import SymRCM
'''
"""
struct reorder_and_slash <: AbstractSolver end
function solveit(solver::reorder_and_slash, A, b)
    d = SymRCM.symrcm(A)
    x = zeros(length(b))
    x[d] = A[d, d] \ b[d]
    return x
end

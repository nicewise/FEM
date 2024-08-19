function displacement_load(; K, free_dofs, load_dofs, load_value, R, factor = 1, solver = LeftDivisionSolver())
    K_ff = K[free_dofs, free_dofs]
    K_fl = K[free_dofs, load_dofs]
    b_f = R - K_fl * load_value * factor
    du = zeros(size(K, 1))
    du[load_dofs] = load_value * factor
    du[free_dofs] = solveit(solver, K_ff, b_f)
    return du
end

function displacement_load(a::AbstractAnalysis; K, factor = 1, solver = a.solver)
    displacement_load(K = K, free_dofs = a.free_dofs, load_dofs = a.load_dofs, load_value = a.load_value, R = a.R, factor = factor, solver = a.solver)
end

function displacement_converge(; K, free_dofs, R, solver = LeftDivisionSolver())
    K_ff = K[free_dofs, free_dofs]
    du = zeros(size(K, 1))
    du[free_dofs] = solveit(solver, K_ff, R)
    return du
end

function displacement_converge(a::AbstractAnalysis; K, solver = a.solver)
    displacement_load(K = K, free_dofs = a.free_dofs, R = a.R, solver = a.solver)
end

struct simple_analysis <: AbstractAnalysis
    free_dofs::Vector{Int}
    fixed_dofs::Vector{Int}
    load_dofs::Vector{Int}
    load_value::Vector{Float64}
    solver::AbstractSolver
    f_ext::Vector{Float64}
    f_int::Vector{Float64}
    R::Vector{Float64}
    tolerance::Float64
end
simple_analysis(fixed_dofs, load_dofs, load_value, dofs; solver = LeftDivisionSolver(), tolerance = 1e-4) = begin
    free_dofs = setdiff(1:dofs, fixed_dofs, load_dofs)
    simple_analysis(free_dofs,
                    fixed_dofs,
                    load_dofs,
                    load_value,
                    solver,
                    zeros(length(free_dofs)),
                    zeros(dofs),
                    zeros(length(free_dofs)),
                    tolerance)
end

#=

mutable struct pd_implicit_static <: AbstractAnalysis
    free_dofs::Vector{Int}
    fixed_dofs::Vector{Int}
    load_dofs::Vector{Int}
    load_value::Vector{Float64}
    K::SparseMatrixCSC{Float64, Int64}
    u_f::Vector{Float64}
    f_ext::Vector{Float64}
    f_int::Vector{Float64}
    R::Vector{Float64}
    load_factor::Vector{Float64}
    time::Vector{Float64}
    broken_bond::Vector{Vector{Int}}
    tolerance::Float64
    convergence::Vector{Bool}
    solver::AbstractSolver
end
pd_implicit_static(fixed_dofs, load_dofs, load_value, dofs; solver::AbstractSolver = LeftDivisionSolver(), tolerance = 1e-3) = begin
    free_dofs = setdiff(1:dofs, fixed_dofs, load_dofs)
    implicit_static_analysis(
        free_dofs,
        fixed_dofs,
        load_dofs,
        load_value,
        spzeros(Float64, 0, 0),
        zeros(length(free_dofs)),
        zeros(length(free_dofs)),
        zeros(dofs),
        zeros(length(free_dofs)),
        Float64[],
        Float64[],
        Vector{Vector{Int}}(),
        tolerance,
        Vector{Bool}(),
        solver
    )
end

function load_attempt(a::pd_implicit_static, factor::Float64; solver::AbstractSolver = a.solver)
    K_ff = a.K[a.free_dofs, a.free_dofs]
    K_fl = a.K[a.free_dofs, a.load_dofs]
    b_f = -K_fl * a.load_value * factor
    u = zeros(size(a.K, 1))
    u[a.load_dofs] = a.load_value * factor
    u[a.free_dofs] = solve(solver, K_ff, b_f)
    return u
end

function R_gen!(a::pd_implicit_static)
    a.R = a.f_ext - a.f_int[a.free_dofs]
    conv = norm(a.R, 2)^2 / (1 + norm(a.f_ext, 2)^2) <= a.tolerance
    return conv
end

function save_step_data_h5(i::Int, fn::String, u, f, d)
    fn = fn * ".h5"
    h5open(fn, "r+") do fid
        ux, uy = u[1:2:end], u[2:2:end]
        fx, fy = f[1:2:end], f[2:2:end]
        write(fid, string(i)*"/displacement/ux", ux)
        write(fid, string(i)*"/displacement/uy", uy)
        write(fid, string(i)*"/force/fx", fx)
        write(fid, string(i)*"/force/fy", fy)
        write(fid, string(i)*"/damage", d)
    end
    return nothing
end
=#

mutable struct implicit_static_analysis
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
    solver::Function
end
implicit_static_analysis(fixed_dofs, load_dofs, load_value, dofs; solver::Function = \, tolerance = 1e-3) = begin
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

function load_attempt(solve::implicit_static_analysis, factor::Float64; solver::Function = solve.solver)
    K_ff = solve.K[solve.free_dofs, solve.free_dofs]
    K_fl = solve.K[solve.free_dofs, solve.load_dofs]
    b_f = -K_fl * solve.load_value * factor
    u = zeros(size(solve.K, 1))
    u[solve.load_dofs] = solve.load_value * factor
    u[solve.free_dofs] = solver(K_ff, b_f)
    return u
end

function R_gen!(solve::implicit_static_analysis)
    solve.R = solve.f_ext - solve.f_int[solve.free_dofs]
    conv = norm(solve.R, 2)^2 / (1 + norm(solve.f_ext, 2)^2) <= solve.tolerance
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

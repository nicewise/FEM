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
    displacement_converge(K = K, free_dofs = a.free_dofs, R = a.R, solver = a.solver)
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


mutable struct testcon <: AbstractRecorder
    σ::Vector{Vector{Float64}}
    ε::Vector{Vector{Float64}}
    statev::Vector{Vector{Float64}}
    iteration::Vector{Int}
    time::Vector{Float64}

    testcon() = new([], [], [], [], [])
end
function record!(rc::testcon, p::AbstractPoint)
    push!(rc.σ, copy(p.σ))
    push!(rc.ε, copy(p.ε))
    push!(rc.statev, copy(p.statev))
    return nothing
end

function testcon!(rc::AbstractRecorder, tc::simple_analysis, p::AbstractPoint, cl::AbstractConstitutiveLaw, nincr; kwargs...)
    for incr in 1:nincr
        incr_time = @elapsed begin
            dε = displacement_load(tc, K = p.D)
            de = dε
            push!(rc.iteration, 0)
            for iter in 1:30
                p.dε .= dε
                constitutive_law_apply!(cl, p; kwargs...)
                tc.f_int .= p.σ
                tc.R .= tc.f_ext - tc.f_int[tc.free_dofs]
                if sqrt( (de'*de) / (dε'*dε) ) <= tc.tolerance
                    rc.iteration[incr] = iter
                    break
                end
                de = displacement_converge(tc, K = p.D)
                dε += de
                rc.iteration[incr] = iter
            end
            p.ε .+= dε
        end
        record!(rc, p)
        push!(rc.time, incr_time)
    end
end


mutable struct fem_recorder <: AbstractRecorder
    fn::String
    output_dir::String
    vtk_dir::String
    interpolation_method
    points
    statev_index
    convergence::Vector{Bool}
    iteration::Vector{Int}
    time::Vector{Float64}
end
function fem_recorder(fn::String, output_dir::String, ps::Vector{<:AbstractPoint};
                      interpolation_method = ScatteredInterpolation.Shepard(),
                      statev_index = 1:length(ps[1].statev))

    output_dir = abspath(output_dir)
    vtk_dir = joinpath(output_dir, fn*"_vtk")
    if !ispath(vtk_dir)
        mkpath(vtk_dir)
    end
    statev_index = isempty(statev_index) ? nothing : statev_index
    position = Vector{Vector{Float64}}(undef, length(ps))
    @threads for i in eachindex(ps)
        position[i] = ps[i].position
    end
    return fem_recorder(fn, output_dir, vtk_dir, interpolation_method, hcat(position...), statev_index, [], [], [])
end
function sort_point_data(rc::fem_recorder, ps::Vector{<:AbstractPoint})
    n1 = length(ps)
    ε = Vector{Vector{Float64}}(undef, n1)
    σ = Vector{Vector{Float64}}(undef, n1)
    statev = Vector{Vector{Float64}}(undef, n1)
    @threads for i in eachindex(ps)
        ε[i] = ps[i].ε
        σ[i] = ps[i].σ
        statev[i] = ps[i].statev
    end
    return ε, σ, statev
end
function record!(rc::fem_recorder,
                 a::AbstractAnalysis,
                 fe::AbstractFem,
                 u;
                 i::Int)
    fn = joinpath(rc.output_dir, rc.fn*".jld2")
    JLD2.jldopen(fn, "a+") do file
        file[string(i)*"/displacement"] = u
        file[string(i)*"/force"] = a.f_int
        ε, σ, statev = sort_point_data(rc, fe.quadrature_points)
        file[string(i)*"/strain"] = ε
        file[string(i)*"/stress"] = σ
        if !isnothing(rc.statev_index)
            file[string(i)*"/statev"] = statev
        end
    end
    return nothing
end
function record!(rc::fem_recorder)
    fn = joinpath(rc.output_dir, rc.fn*".jld2")
    JLD2.jldopen(fn, "a+") do file
        file["convergence"] = rc.convergence
        file["iteration"] = rc.iteration
        file["time"] = rc.time
    end
    return nothing
end

function displacement_control!(rc::AbstractRecorder,
                               a::simple_analysis,
                               fe::AbstractFem,
                               mesh::AbstractMesh,
                               cl::AbstractConstitutiveLaw;
                               nincr::Int,
                               miter::Int = 30,
                               u::Vector{Float64} = zeros(mesh.dofs),
                               K_u = elastic_initialization!(fe, mesh),
                               kwargs...)

    for incr in 1:nincr
        incr_time = @elapsed begin
            du = displacement_load(a, K = K_u)
            converged = false
            for iter in 1:miter
                displacement_apply!(fe, du = du)
                constitutive_law_apply!(cl, fe, kwargs...)
                K_u, f_int = elastic_system_assemble!(fe, mesh)
                u += du
                displacement_apply!(fe, u = u)
                a.f_int .= f_int
                a.R .= a.f_ext - a.f_int[a.free_dofs]
                du = displacement_converge(a, K = K_u)
                if norm(du) / norm(u) <= a.tolerance && norm(a.R) / norm(a.f_int) <= a.tolerance
                    converged = true
                    push!(rc.iteration, iter)
                    print("Incriment $incr converged at $iter", " th iteration. ")
                    break
                end
            end
            if !converged
                push!(rc.iteration, miter)
                @warn "The solution did not converge within $miter iterations."
            end
            push!(rc.convergence, converged)
            record!(rc, a, fe, u, i = incr)
        end
        push!(rc.time, incr_time)
        println("Time comsume: $incr_time", "s.")
    end
end

function gradient_displacement_control!(rc::AbstractRecorder,
                                        a::simple_analysis,
                                        fe::AbstractFem,
                                        mesh::AbstractMesh,
                                        g::AbstractConstitutiveLaw;
                                        nincr::Int,
                                        miter::Int = 30,
                                        u::Vector{Float64} = zeros(mesh.dofs),
                                        K_u = elastic_initialization!(fe, mesh),
                                        K_g = gradient_stiffness_assemble(g, fe, mesh)
                                        )

    for incr in 1:nincr
        incr_time = @elapsed begin
            du = displacement_load(a, K = K_u)
            converged = false
            for iter in 1:miter
                displacement_apply!(fe, du = du)
                constitutive_law_apply!(g, fe, K_g, mesh)
                K_u, f_int = elastic_system_assemble!(fe, mesh)
                u += du
                displacement_apply!(fe, u = u)
                a.f_int .= f_int
                a.R .= a.f_ext - a.f_int[a.free_dofs]
                du = displacement_converge(a, K = K_u)
                if norm(du) / norm(u) <= a.tolerance && norm(a.R) / norm(a.f_int) <= a.tolerance
                    converged = true
                    push!(rc.iteration, iter)
                    print("Incriment $incr converged at $iter", " th iteration. ")
                    break
                end
            end
            if !converged
                push!(rc.iteration, miter)
                @warn "The solution did not converge within $miter iterations."
            end
            push!(rc.convergence, converged)
            record!(rc, a, fe, u, i = incr)
        end
        push!(rc.time, incr_time)
        println("Time comsume: $incr_time", "s.")
    end
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

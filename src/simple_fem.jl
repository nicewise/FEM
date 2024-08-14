struct simple_fem{T<:FE} <: AbstractFem{T}
    element_table::Matrix{Int}
    dof_per_element::Int
    number_of_quadrature::Int
    ref_coords::Matrix{Float64}
    weights::Vector{Float64}
    N::Vector{Vector{Float64}}
    E::Float64
    ν::Float64
    planestrain::Bool# 3d情况下，这一项没有意义
end

function simple_fem_gen(mode, e, x...; kwargs...)
    params = Dict(
        :planestrain => false,
        :E => 0.0,
        :ν => 0.0
    )
    for (key, value) in kwargs
        if haskey(params, key)
            params[key] = value
        else
            throw(ArgumentError("Invalid keyword argument: $key"))
        end
    end

    planestrain = params[:planestrain]
    E = params[:E]
    ν = params[:ν]

    if E == 0.0 || ν == 0.0
        throw(ArgumentError("Material properties E and ν must be non-zero."))
    end

    types_dict = Dict(
        "t1" => fet1, "t2" => fet2, "t3" => fet3,
        "q1" => feq1, "q2" => feq2, "q3" => feq3,
        "T1" => feT1, "T2" => feT2, "T3" => feT3,
        "H1" => feH1, "H2" => feH2, "H3" => feH3,
        "P1" => feP1, "P2" => feP2, "P3" => feP3
    )

    # 获取类型
    T = get(types_dict, mode, nothing)
    if T === nothing
        throw(ArgumentError("Invalid mode: $mode"))
    end

    is_2d = islowercase(mode[1])

    # 验证单元类型对应的节点数
    valid_node_counts = Dict(
        "t1" => 3, "t2" => 6, "t3" => 10,
        "q1" => 4, "q2" => 9, "q3" => 16,
        "T1" => 4, "T2" => 8, "T3" => 27,
        "H1" => 8, "H2" => 8, "H3" => 27,
        "P1" => 6, "P2" => 8, "P3" => 27
    )

    node_per_element = size(e, 1)
    expected_nodes = valid_node_counts[mode]

    if node_per_element != expected_nodes
        throw(ArgumentError("For mode '$mode', expected $expected_nodes nodes per element, but found $node_per_element."))
    end

    # 2D 和 3D 的自由度
    dof_per_element = (is_2d ? 2 : 3) * node_per_element

    ref_coords, weights = quad_form(mode[1], x...)
    number_of_quadrature = length(weights)
    N = [N_gen(T, ref_coords[:, i]) for i in 1:number_of_quadrature]

    simple_fem{T}(e,
                  dof_per_element,
                  length(weights),
                  ref_coords,
                  weights,
                  N,
                  E,
                  ν,
                  planestrain)
end

function element_stiffness_gen(fe::simple_fem{T}, mesh::mesh,  element_index) where T <: FE
    node_index = fe.element_table[:, element_index]
    coords = mesh.nodes[:, node_index]
    cl = if T <: fe3d
        linear_elastic{fe3d}(fe.E, fe.ν)
    elseif fe.planestrain
        linear_elastic{fe2d}(fe.E, fe.ν; planestrain = true)
    else
        linear_elastic{fe2d}(fe.E, fe.ν; planestrain = false)
    end
    V = zeros(fe.dof_per_element, fe.dof_per_element)
    for i in 1:fe.number_of_quadrature
        B, JxW = B_gen(T, fe.ref_coords[:, i], fe.weights[i], coords)
        V += B' * cl.D * B * JxW * mesh.thick
    end
    dof = hcat(2*node_index .- 1, 2*node_index)'[:]
    I = repeat(dof, fe.dof_per_element)
    J = repeat(dof, inner=fe.dof_per_element)
    return(I, J, V[:])
end

function stiffness_assemble(fe::simple_fem{<:FE}, mesh::mesh)
    n = fe.dof_per_element^2
    element_number = size(fe.element_table, 2)
    Is = zeros(n*element_number)
    Js = zeros(n*element_number)
    Vs = zeros(n*element_number)
    for i in 1:element_number
        I, J, V = element_stiffness_gen(fe, mesh, i)
        Is[n*(i-1)+1:n*i] = I
        Js[n*(i-1)+1:n*i] = J
        Vs[n*(i-1)+1:n*i] = V
    end
    sparse(Is, Js, Vs, mesh.dofs, mesh.dofs)
end

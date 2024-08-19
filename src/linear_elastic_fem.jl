struct linear_elastic{C<:AbstractCellType, D<:AbstractDimension} <: AbstractFem{C, D}
    element_table::Matrix{Int}
    dof_per_element::Int
    number_of_quadrature::Int
    ref_coords::Matrix{Float64}
    weights::Vector{Float64}
    N::Vector{Vector{Float64}}
    E::Float64
    ν::Float64
end

function linear_elastic(T::Type, e, x...; kwargs...)
    node_per_element = assert_vertex_number(T, e)

    params = Dict(
        :isplanestrain => nothing,
        :E => 0.0,
        :ν => 0.0
    )

    # 合并传入的关键字参数，并确保没有无效的关键字
    params = merge(params, Dict(kwargs))
    unknown_keys = setdiff(keys(kwargs), keys(params))
    @assert isempty(unknown_keys) "Invalid keyword arguments: $(collect(unknown_keys))"

    isplanestrain = params[:isplanestrain]
    E = params[:E]
    ν = params[:ν]

    # 检查 `isplanestrain` 是否为布尔值或为空
    @assert isnothing(isplanestrain) || isa(isplanestrain, Bool) "Keyword argument `isplanestrain` must be a Bool. Provided: $isplanestrain"

    @assert isa(E, Float64) && E != 0.0 "Material property `E` must be a non-zero Float64. Provided: $E"
    @assert isa(ν, Float64) && ν != 0.0 "Material property `ν` must be a non-zero Float64. Provided: $ν"

    # 确定维度类型
    D = isnothing(isplanestrain) ? Dim3 :
        isplanestrain ? PlaneStrain : PlaneStress

    # 检查2D问题是否正确设置了 D
    is2d = T <: Cell2
    @assert !(is2d && D <: Dim3) "For $T cell, specify whether it's planestrain or not."

    # 计算自由度
    dof_per_element = (is2d ? 2 : 3) * node_per_element

    # 生成积分点和形函数
    ref_coords, weights = quad_form(T, x...)
    number_of_quadrature = length(weights)
    N = [N_gen(T, ref_coords[:, i]) for i in 1:number_of_quadrature]

    # 返回有限元模型对象
    return linear_elastic{T, D}(e,
                                dof_per_element,
                                number_of_quadrature,
                                ref_coords,
                                weights,
                                N,
                                E,
                                ν)
end
function element_stiffness_gen(fe::linear_elastic{C, D}, mesh, element_index) where {C<:AbstractCellType, D<:AbstractDimension}
    node_index = fe.element_table[:, element_index]
    coords = mesh.nodes[:, node_index]
    cl = constitutive_linear_elastic{D}(fe.E, fe.ν)
    V = zeros(fe.dof_per_element, fe.dof_per_element)
    for i in 1:fe.number_of_quadrature
        B, JxW = B_gen(C, D, fe.ref_coords[:, i], fe.weights[i], coords)
        V += B' * cl.D * B * JxW * mesh.thick
    end
    dof = vec(hcat(2*node_index .- 1, 2*node_index)')
    I = repeat(dof, fe.dof_per_element)
    J = repeat(dof, inner=fe.dof_per_element)
    return(I, J, V[:])
end

function stiffness_assemble(fe::linear_elastic{<:AbstractCellType, <:AbstractDimension}, mesh::AbstractMesh)
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

struct quadrature_point{T<:FE} <: AbstractQuadraturePoint{T}
    position::Vector{Float64}
    ε::Vector{Float64}
    dε::Vector{Float64}
    σ::Vector{Float64}
    D::Matrix{Float64}
    B::Matrix{Float64}
    JxW::Float64
    statev::Vector{Float64}
end
function quadrature_point_gen(T::Type, position, B, JxW, nstatev)
    dim = T <: fe2d ? 2 : 3
    nstress = T <: fe2d ? 3 : 6

    ε = zeros(nstress)
    dε = zeros(nstress)
    σ = zeros(nstress)
    D = zeros(nstress, nstress)
    statev = zeros(nstatev)

    quadrature_point{T}(position, ε, dε, σ, D, B, JxW, statev)
end

function quadrature_points_gen(T::Type,
                               coords,
                               number_of_quadrature,
                               ref_coords,
                               weight,
                               N,
                               nstatev)
    quadrature_points = Vector{quadrature_point{T}}(undef, number_of_quadrature)

    for i in 1:number_of_quadrature
        B, JxW = B_gen(T, ref_coords[:, i], weight[i], coords)
        quadrature_points[i] = quadrature_point_gen(T, coords * N[i], B, JxW, nstatev)
    end

    return quadrature_points
end

function quadrature_points_gen(T::Type,
                               element_table,
                               nodes,
                               number_of_quadrature,
                               ref_coords,
                               weight,
                               N,
                               nstatev)
    element_number = size(element_table, 2)
    quadrature_points = Matrix{quadrature_point{T}}(undef,
                                                    number_of_quadrature,
                                                    element_number)

    @threads for i in 1:element_number
        quadrature_points[:, i] = quadrature_points_gen(T,
                                                        nodes[:, element_table[:, i]],
                                                        number_of_quadrature,
                                                        ref_coords,
                                                        weight,
                                                        N,
                                                        nstatev)
    end
    return vec(quadrature_points)
end
struct element{T<:FE} <: AbstractElement{T}
    quadrature_points::Vector{quadrature_point{<:FE}}
    dof::Vector{Int}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{Float64}
    F::Vector{Float64}
end
function element_gen(T::Type, node_index, i, nq, quadrature_points)
    dof = vec(hcat(2node_index.-1, 2node_index)')
    ndof = length(dof)
    I = repeat(dof, outer = ndof)
    J = repeat(dof, inner = ndof)
    V = similar(I, Float64)
    F = zeros(ndof)
    element{T}(quadrature_points[(i-1)*nq+1:i*nq],
               collect(dof),
               I,
               J,
               V,
               F)
end
function elements_gen(T::Type,
                      element_table,
                      number_of_quadrature,
                      quadrature_points)
    elements = Vector{element}(undef, size(element_table, 2))
    @threads for i in axes(element_table, 2)
        elements[i] = element_gen(T,
                                  element_table[:, i],
                                  i,
                                  number_of_quadrature,
                                  quadrature_points)
    end
    return elements
end

struct fem{T<:FE} <: AbstractFem{T}
    element_table::Matrix{Int}
    dof_per_element::Int
    number_of_quadrature::Int
    ref_coords::Matrix{Float64}
    weights::Vector{Float64}
    N::Vector{Vector{Float64}}
    elements::Vector{element}
    quadrature_points::Vector{quadrature_point}
    E::Float64
    ν::Float64
    planestrain::Bool# 3d情况下，这一项没有意义
end

function fem_gen(mode, e, n, x...; kwargs...)
    params = Dict(
        :nstatev => 0,
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

    nstatev = params[:nstatev]
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

    quadrature_points = quadrature_points_gen(T, e, n, number_of_quadrature,
                                              ref_coords, weights, N, nstatev)
    elements = elements_gen(T, e, number_of_quadrature, quadrature_points)
    fem{T}(e,
           dof_per_element,
           number_of_quadrature,
           ref_coords,
           weights,
           N,
           elements,
           quadrature_points,
           E,
           ν,
           planestrain)
end

function element_stiffness_gen!(element::AbstractElement{<:FE}, dpe::Int, thick::Float64)
    ndof = length(element.dof)
    V = zeros(dpe, ndof)
    for qpoint in element.quadrature_points
        V += qpoint.B' * qpoint.D * qpoint.B * qpoint.JxW * thick
    end
    element.V .= vec(V)
    return nothing
end

function elements_stiffness_gen!(fe::AbstractFem{T}, thick::Float64) where T<:FE
    thick = T <: fe2d ? thick : 1.0
    @threads for element in fe.elements
        element_stiffness_gen!(element, fe.dof_per_element, thick)
    end
end
function K_gen(fe::fem{<:FE}, dofs::Int)
    n = fe.dof_per_element^2
    element_number = size(fe.element_table, 2)
    Is = Vector{Int}(undef, n*element_number)
    Js = Vector{Int}(undef, n*element_number)
    Vs = Vector{Float64}(undef, n*element_number)
    @threads for i in 1:element_number
        @inbounds begin
            Is[n*(i-1)+1:n*i] = fe.elements[i].I
            Js[n*(i-1)+1:n*i] = fe.elements[i].J
            Vs[n*(i-1)+1:n*i] = fe.elements[i].V
        end
    end
    sparse(Is, Js, Vs, dofs, dofs)    
end
stiffness_assemble(fe::fem{<:FE}, dofs::Int) = K_gen(fe, dofs)
function elastic_initialization!(fe::AbstractFem{T}, mesh::mesh, E::Float64, ν::Float64) where T<:FE
    cl = if T <: fe3d
        linear_elastic{fe3d}(E, ν)
    elseif fe.planestrain
        linear_elastic{fe2d}(E, ν; planestrain = true)
    else
        linear_elastic{fe2d}(E, ν; planestrain = false)
    end

    constitutive_law_apply!(fe, cl)
    elements_stiffness_gen!(fe, mesh.thick)
    K_gen(fe, mesh.dofs)
end
elastic_initialization!(fe::AbstractFem{<:FE}, mesh::mesh) = elastic_initialization!(fe, mesh, fe.E, fe.ν)

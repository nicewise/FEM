struct quadrature_point{C<:AbstractCellType, Dim<:AbstractDimension} <: AbstractPoint{C, Dim}
    position::Vector{Float64}
    ε::Vector{Float64}
    dε::Vector{Float64}
    σ::Vector{Float64}
    D::Matrix{Float64}
    D4::Matrix{Float64} # only valid for 2d model
    B::Matrix{Float64}
    JxW::Float64
    statev::Vector{Float64}
end
function simple_point_gen(; C::Type = t1, Dim::Type = Dim3, nstatev = 0)
    return quadrature_point_gen(C,
                                Dim,
                                Float64[],
                                zeros(0, 0),
                                0.0,
                                nstatev)
end
function quadrature_point_gen(C::Type, Dim::Type, position, B, JxW, nstatev)
    dim = Dim <: Dim2 ? 2 : 3
    nstress = Dim <: Dim2 ? 3 : 6

    ε = zeros(nstress)
    dε = zeros(nstress)
    σ = zeros(nstress)
    D = zeros(nstress, nstress)
    D4 = zeros(4, 4)
    statev = zeros(nstatev)

    quadrature_point{C, Dim}(position, ε, dε, σ, D, D4, B, JxW, statev)
end

function quadrature_points_gen(C::Type,
                               D::Type,
                               coords,
                               number_of_quadrature,
                               ref_coords,
                               weight,
                               N,
                               nstatev)
    quadrature_points = Vector{quadrature_point{C, D}}(undef, number_of_quadrature)

    for i in 1:number_of_quadrature
        B, JxW = B_gen(C, D, ref_coords[:, i], weight[i], coords)
        quadrature_points[i] = quadrature_point_gen(C, D, coords * N[i], B, JxW, nstatev)
    end

    return quadrature_points
end

function quadrature_points_gen(C::Type,
                               D::Type,
                               element_table,
                               nodes,
                               number_of_quadrature,
                               ref_coords,
                               weight,
                               N,
                               nstatev)
    element_number = size(element_table, 2)
    quadrature_points = Matrix{quadrature_point{C, D}}(undef,
                                                       number_of_quadrature,
                                                       element_number)

    @threads for i in 1:element_number
        quadrature_points[:, i] = quadrature_points_gen(C,
                                                        D,
                                                        nodes[:, element_table[:, i]],
                                                        number_of_quadrature,
                                                        ref_coords,
                                                        weight,
                                                        N,
                                                        nstatev)
    end
    return vec(quadrature_points)
end


struct element{C<:AbstractCellType, D<:AbstractDimension} <: AbstractCell{C, D}
    quadrature_points::Vector{quadrature_point{C, D}}
    dof::Vector{Int}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{Float64}
    F::Vector{Float64}
end
function element_gen(C::Type, D::Type, node_index, i, nq, quadrature_points)
    dof = D <: Dim2 ? vec(hcat(2node_index.-1, 2node_index)') : vec(hcat(3node_index.-2, 3node_index.-1, 3node_index)')
    ndof = length(dof)
    I = repeat(dof, outer = ndof)
    J = repeat(dof, inner = ndof)
    V = similar(I, Float64)
    F = zeros(ndof)
    element{C, D}(quadrature_points[(i-1)*nq+1:i*nq],
                  collect(dof),
                  I,
                  J,
                  V,
                  F)
end
function elements_gen(C::Type,
                      D::Type,
                      element_table,
                      number_of_quadrature,
                      quadrature_points)
    elements = Vector{element}(undef, size(element_table, 2)) # FIXME
    @threads for i in axes(element_table, 2)
        elements[i] = element_gen(C,
                                  D,
                                  element_table[:, i],
                                  i,
                                  number_of_quadrature,
                                  quadrature_points)
    end
    return elements
end

struct elastoplastic{C<:AbstractCellType, D<:AbstractDimension} <: AbstractFem{C, D}
    element_table::Matrix{Int}
    dof_per_element::Int
    number_of_quadrature::Int
    ref_coords::Matrix{Float64}
    weights::Vector{Float64}
    N::Vector{Vector{Float64}}
    elements::Vector{element{C, D}}
    quadrature_points::Vector{quadrature_point{C, D}}
    E::Float64
    ν::Float64
end

function elastoplastic(T::Type, e, n, x...; kwargs...)
   node_per_element = assert_vertex_number(T, e)

    params = Dict(
        :nstatev => 0,
        :planestrain => false,
        :E => 0.0,
        :ν => 0.0
    )

    # 合并传入的关键字参数，并确保没有无效的关键字
    params = merge(params, Dict(kwargs))
    unknown_keys = setdiff(keys(kwargs), keys(params))
    @assert isempty(unknown_keys) "Invalid keyword arguments: $(collect(unknown_keys))"

    nstatev = params[:nstatev]
    isplanestrain = params[:isplanestrain]
    E = params[:E]
    ν = params[:ν]

    # 检查 `isplanestrain` 是否为布尔值或为空
    @assert isnothing(isplanestrain) || isa(isplanestrain, Bool) "Keyword argument `isplanestrain` must be a Bool. Provided: $isplanestrain"

    @assert isa(E, Float64) && E != 0.0 "Material property `E` must be a non-zero Float64. Provided: $E"
    @assert isa(ν, Float64) && ν != 0.0 "Material property `ν` must be a non-zero Float64. Provided: $ν"
    @assert isa(nstatev, Int64) "Number of state variable `nstatev` must be a Int64. Provided: $nstatev"

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

    quadrature_points = quadrature_points_gen(T, D, e, n, number_of_quadrature,
                                              ref_coords, weights, N, nstatev)
    elements = elements_gen(T, D, e, number_of_quadrature, quadrature_points)
    elastoplastic{T, D}(e,
                        dof_per_element,
                        number_of_quadrature,
                        ref_coords,
                        weights,
                        N,
                        elements,
                        quadrature_points,
                        E,
                        ν)
end

function element_stiffness_gen!(element::AbstractCell, dpe::Int, thick::Float64)
    ndof = length(element.dof)
    V = zeros(dpe, ndof)
    for qpoint in element.quadrature_points
        V += qpoint.B' * qpoint.D * qpoint.B * qpoint.JxW * thick
    end
    element.V .= vec(V)
    return nothing
end

function elements_stiffness_gen!(fe::AbstractFem{<:AbstractCellType, D}, thick::Float64) where D<:AbstractDimension
    thick = D <: Dim2 ? thick : 1.0
    @threads for element in fe.elements
        element_stiffness_gen!(element, fe.dof_per_element, thick)
    end
end
function K_gen(fe::AbstractFem, dofs::Int)
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
stiffness_assemble(fe::AbstractFem, dofs::Int) = K_gen(fe, dofs)

function constitutive_law_apply!(F::AbstractConstitutiveLaw, fe::AbstractFem)
    @threads for point in fe.quadrature_points
        constitutive_law_apply!(F, point)
    end
    return nothing
end
function elastic_initialization!(fe::AbstractFem{<:AbstractCellType, D}, mesh, E::Float64, ν::Float64) where D<:AbstractDimension
    cl = constitutive_linear_elastic{D}(E, ν)
    constitutive_law_apply!(cl, fe)
    elements_stiffness_gen!(fe, mesh.thick)
    K_gen(fe, mesh.dofs)
end
elastic_initialization!(fe::AbstractFem, mesh) = elastic_initialization!(fe, mesh, fe.E, fe.ν)

function displacement_apply!(e::element, u::Vector{Float64})
    for qpoint in e.quadrature_points
        qpoint.ε .= qpoint.B * u
    end
end
function delta_displacement_apply!(e::element, u::Vector{Float64})
    for qpoint in e.quadrature_points
        qpoint.dε .= qpoint.B * u
    end
end
function displacement_apply!(fe::AbstractFem; u::Vector{Float64}=Float64[], du::Vector{Float64}=Float64[])
    if isempty(u)
        @threads for element in fe.elements
            delta_displacement_apply!(element, du[element.dof])
        end
    elseif isempty(du)
        @threads for element in fe.elements
            displacement_apply!(element, u[element.dof])
        end
    else
        @threads for element in fe.elements
            displacement_apply!(element, u[element.dof])
            delta_displacement_apply!(element, du[element.dof])
        end
    end
end

function gradient_stiffness_gen(g::AbstractConstitutiveLaw{D}, fe::AbstractFem{C, D}, mesh, element_index)  where {C<:AbstractCellType, D<:AbstractDimension}
    dof = fe.element_table[:, element_index]
    dof_per_element = length(dof)
    I = repeat(dof, outer=dof_per_element)
    J = repeat(dof, inner=dof_per_element)
    coords = mesh.nodes[:, dof]
    V = zeros(dof_per_element, dof_per_element)
    for i in 1:fe.number_of_quadrature
        B, JxW = B_gen(C, fe.ref_coords[:, i], fe.weights[i], coords)
        V += JxW * (fe.N[i] * fe.N[i]' + g.a * B * B')
    end
    return(I, J, V[:])
end
function gradient_stiffness_assemble(g::AbstractConstitutiveLaw, fe::AbstractFem, mesh)
    n, element_number = size(fe.element_table)
    n = n^2
    dofs = n*element_number
    Is = zeros(dofs)
    Js = zeros(dofs)
    Vs = zeros(dofs)
    for i in 1:element_number
        I, J, V = gradient_stiffness_gen(g, fe, mesh, i)
        Is[n*(i-1)+1:n*i] = I
        Js[n*(i-1)+1:n*i] = J
        Vs[n*(i-1)+1:n*i] = V
    end
    sparse(Is, Js, Vs, size(mesh.nodes, 2), size(mesh.nodes, 2))
end

function gradient_rhs_gen(g::AbstractConstitutiveLaw, p::AbstractPoint, N)
    return N * p.statev[g.x]' * p.JxW
end
function gradient_rhs_gen(g::AbstractConstitutiveLaw, fe::AbstractFem, element_index::Int)
    n = size(fe.element_table, 1)
    o = zeros(n, length(g.x))
    for i in 1:fe.number_of_quadrature
        o += gradient_rhs_gen(g, fe.elements[element_index].quadrature_points[i], fe.N[i])
    end
    return o
end
function gradient_rhs_gen(g::AbstractConstitutiveLaw, fe::AbstractFem, mesh)
    o = zeros(size(mesh.nodes, 2), length(g.x))
    for i in axes(fe.element_table, 2)
        dof = fe.element_table[:, i]
        o[dof, :] += gradient_rhs_gen(g, fe, i)
    end
    return o
end

function gradient_apply!(g, fe, u)
    @threads for i in eachindex(fe.elements)
        gradient_apply!(g, fe.elements[i], fe.N, u[fe.element_table[:, i], :])
    end
end
function gradient_apply!(g, e, N, u)
    for i in eachindex(e.quadrature_points)
        for (j, x) in enumerate(u' * N[i])
            @inbounds e.quadrature_points[i].statev[g.y[j]] = x
        end
    end
end

function constitutive_law_apply!(g::AbstractConstitutiveLaw,
                                 fe::AbstractFem,
                                 K_g,
                                 mesh::AbstractMesh)
    constitutive_law_apply!(g, fe, before_gradient=true)
    rhs = K_g \ gradient_rhs_gen(g, fe, mesh)
    gradient_apply!(g, fe, rhs)
    constitutive_law_apply!(g, fe, before_gradient=false)
    return nothing
end

struct fem_mesh{D<:AbstractDimension} <: AbstractMesh{D}
    nodes::Matrix{Float64}
    dofs::Int
    thick::Float64
end
function fem_mesh_gen(nodes; thick = nothing)
    dimension, node_number = size(nodes) # 2xnnodes或3xnnodes
    # 维度检查和厚度验证
    @assert !(dimension == 2 && isnothing(thick)) "For 2D problem, `thick` must be provided."
    thick = dimension == 2 ? thick : 1.0
    dofs = dimension * node_number
    D = dimension == 1 ? Dim1 :
        dimension == 2 ? Dim2 :
        dimension == 3 ? Dim3 :
        nothing
    fem_mesh{D}(nodes, dofs, thick)
end
#=
struct pd_mesh{D<:AbstractDimension} <: AbstractMesh{D}
    nodes::Matrix{Float64}
    node_number::Int
    dofs::Int
    inter_crack::Function
    thick::Float64
    δ::Float64
    dx::Float64
end
function pd_mesh_gen(nodes, elements, mat; thick = 0.0, delta = 0.0, dx = 0.0, crack = zeros(0, 4))
    dimension, node_number = size(nodes) # 2xnnodes或3xnnodes

    # 维度检查和厚度验证
    if dimension == 2 && thick === 0.0
        throw(ArgumentError("For 2D problem, `thick` must be provided."))
    end

    dofs = dimension * node_number
    inter_crack = create_inter_crack(crack)
    mesh(nodes, elements, node_number, dimension, dofs, material_para(mat), inter_crack, thick, delta, dx)
end

function create_inter_crack(crack::Matrix{Float64})
    if isempty(crack)
        return (n1::Vector{Float64}, n2::Vector{Float64}) -> true
    end

    return (n1::Vector{Float64}, n2::Vector{Float64}) -> begin
        ax, ay = n1[1], n1[2]
        bx, by = n2[1], n2[2]

        for k in 1:size(crack, 1)
            cx, cy = crack[k, 1], crack[k, 2]
            dx, dy = crack[k, 3], crack[k, 4]

            u = (cx - ax) * (by - ay) - (bx - ax) * (cy - ay)
            v = (dx - ax) * (by - ay) - (bx - ax) * (dy - ay)
            w = (ax - cx) * (dy - cy) - (dx - cx) * (ay - cy)
            z = (bx - cx) * (dy - cy) - (dx - cx) * (by - cy)

            if u * v < 0 && w * z < 0
                return false
            end
        end

        return true
    end
end

function save_model_h5(mesh::mesh, fn::String)
    fn = fn * ".h5"
    h5open(fn, "w") do fid
        write(fid, "nodes", mesh.nodes)
        write(fid, "elements", mesh.elements)
    end
end
=#

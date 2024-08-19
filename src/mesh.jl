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
        Dimension == 2 ? Dim3 :
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
function vtk_gen_from_h5_binary(fn_base::String, i::Int)
    folder, fn = splitdir(fn_base)
    fn_h5 = fn_base * ".h5"
    @assert ispath(fn_h5) "HDF5 file $fn_h5 does not exist."

    folder_vtk = joinpath(folder, "vtk")

    if !isdir(folder_vtk)
        mkpath(folder_vtk)
    end

    fn_vtk = joinpath(folder_vtk, string(fn, "_", lpad(i, 4, '0'), ".vtk"))
    open(fn_vtk, "w") do fid
        # Write VTK file header
        write(fid, "# vtk DataFile Version 3.0\n")
        write(fid, "VTK from Julia\n")
        write(fid, "BINARY\n")

        write(fid, "DATASET UNSTRUCTURED_GRID\n")

        # Read and write point data
        nodes = h5read(fn_h5, "/nodes")
        num_nodes = size(nodes, 2)
        write(fid, "POINTS $num_nodes float\n")
        # Convert point data to Little Endian byte order and write to file
        nodes_with_zeros = vcat(nodes, zeros(1, num_nodes))
        write(fid, nodes_with_zeros)

        # Read and write cell data
        elements = h5read(fn_h5, "/elements")
        num_elements = size(elements, 2)
        num_indices = length(elements)
        write(fid, "CELLS $num_elements $(num_indices + num_elements)\n")
        # Prepare cell data and write to file
        cell_data = vcat(fill(4, 1, num_elements), elements .- 1)
        write(fid, cell_data)

        # Write cell types (9 represents quadrilateral)
        write(fid, "CELL_TYPES $num_elements\n")
        cell_types = fill(9, num_elements)
        write(fid, cell_types)

        # Read and write specific timestep data
        ux = h5read(fn_h5, string(i)*"/displacement/ux")
        uy = h5read(fn_h5, string(i)*"/displacement/uy")
        fx = h5read(fn_h5, string(i)*"/force/fx")
        fy = h5read(fn_h5, string(i)*"/force/fy")
        d = h5read(fn_h5, string(i)*"/damage")

        # Write point data - displacements
        write(fid, "POINT_DATA $num_nodes\n")
        write(fid, "VECTORS displacements float\n")
        # Prepare displacement data and write to file
        u = vcat(ux', uy', zeros(1, num_nodes))
        write(fid, u)

        # Write point data - node forces
        write(fid, "VECTORS force float\n")
        # Prepare force data and write to file
        force = vcat(fx', fy', zeros(1, num_nodes))
        write(fid, force)

        # Write point data - damage
        write(fid, "SCALARS damage float 1\n")
        write(fid, "LOOKUP_TABLE default\n")
        # Convert damage data to Little Endian byte order and write to file
        write(fid, d)
    end
    # println("Binary VTK file successfully generated")
end

function vtk_gen_from_h5_ascii(fn_base::String, i::Int)
    folder, fn = splitdir(fn_base)
    fn_h5 = fn_base * ".h5"
    @assert ispath(fn_h5) "HDF5 file $fn_h5 does not exist."

    folder_vtk = joinpath(folder, "vtk")

    if !isdir(folder_vtk)
        mkpath(folder_vtk)
    end

    fn_vtk = joinpath(folder_vtk, string(fn, "_", lpad(i, 4, '0'), ".vtk"))
    open(fn_vtk, "w") do fid
        # Write VTK file header
        write(fid, "# vtk DataFile Version 3.0\n")
        write(fid, "VTK from Julia\n")
        write(fid, "ASCII\n")

        write(fid, "DATASET UNSTRUCTURED_GRID\n")

        # Read and write point data
        nodes = h5read(fn_h5, "/nodes")
        num_nodes = size(nodes, 2)
        write(fid, "POINTS $num_nodes float\n")
        # Prepare point data and write to file
        nodes_with_zeros = vcat(nodes, zeros(1, num_nodes))
        write(fid, join(string.(nodes_with_zeros), " ") * "\n")

        # Read and write cell data
        elements = h5read(fn_h5, "/elements")
        num_elements = size(elements, 2)
        num_indices = length(elements)
        write(fid, "CELLS $num_elements $(num_indices + num_elements)\n")
        # Prepare cell data and write to file
        cell_data = vcat(fill(4, 1, num_elements), elements .- 1)
        write(fid, join(string.(cell_data), " ") * "\n")

        # Write cell types (9 represents quadrilateral)
        write(fid, "CELL_TYPES $num_elements\n")
        cell_types = fill(9, num_elements)
        write(fid, join(string.(cell_types), " ") * "\n")

        # Read and write specific timestep data
        ux = h5read(fn_h5, string(i)*"/displacement/ux")
        uy = h5read(fn_h5, string(i)*"/displacement/uy")
        fx = h5read(fn_h5, string(i)*"/force/fx")
        fy = h5read(fn_h5, string(i)*"/force/fy")
        d = h5read(fn_h5, string(i)*"/damage")

        # Write point data - displacements
        write(fid, "POINT_DATA $num_nodes\n")
        write(fid, "VECTORS displacements float\n")
        # Prepare displacement data and write to file
        u = vcat(ux', uy', zeros(1, num_nodes))
        write(fid, join(string.(u), " ") * "\n")

        # Write point data - node forces
        write(fid, "VECTORS force float\n")
        # Prepare force data and write to file
        force = vcat(fx', fy', zeros(1, num_nodes))
        write(fid, join(string.(force), " ") * "\n")

        # Write point data - damage
        write(fid, "SCALARS damage float 1\n")
        write(fid, "LOOKUP_TABLE default\n")
        write(fid, join(string.(d), " ") * "\n")
    end
    # println("Binary VTK file successfully generated")
end

function vtk_gen_from_h5_ascii2(fn_base::String, i::Int)
    folder, fn = splitdir(fn_base)
    fn_h5 = fn_base * ".h5"
    @assert ispath(fn_h5) "HDF5 file $fn_h5 does not exist."

    folder_vtk = joinpath(folder, "vtk")

    if !isdir(folder_vtk)
        mkpath(folder_vtk)
    end

    fn_vtk = joinpath(folder_vtk, string(fn, "_", lpad(i, 4, '0'), ".vtk"))
    open(fn_vtk, "w") do fid
        # Write VTK file header
        write(fid, "# vtk DataFile Version 3.0\n")
        write(fid, "VTK from Julia\n")
        write(fid, "ASCII\n")

        write(fid, "DATASET UNSTRUCTURED_GRID\n")

        # Read and write point data
        nodes = h5read(fn_h5, "/nodes")
        num_nodes = size(nodes, 2)
        write(fid, "POINTS $num_nodes float\n")
        # Prepare point data and write to file
        nodes_with_zeros = vcat(nodes, zeros(1, num_nodes))
        write(fid, join(string.(nodes_with_zeros), " ") * "\n")
        #write(fid, join(string.(nodes), " ") * "\n")

        # Read and write cell data
        elements = h5read(fn_h5, "/elements")
        num_elements = size(elements, 2)
        num_indices = length(elements)
        write(fid, "CELLS $num_elements $(num_indices + num_elements)\n")
        # Prepare cell data and write to file
        cell_data = vcat(fill(4, 1, num_elements), elements .- 1)
        write(fid, join(string.(cell_data), " ") * "\n")

        # Write cell types (9 represents quadrilateral)
        write(fid, "CELL_TYPES $num_elements\n")
        cell_types = fill(9, num_elements)
        write(fid, join(string.(cell_types), " ") * "\n")

        # Read and write specific timestep data
        ux = h5read(fn_h5, string(i)*"/displacement/ux")
        uy = h5read(fn_h5, string(i)*"/displacement/uy")
        fx = h5read(fn_h5, string(i)*"/force/fx")
        fy = h5read(fn_h5, string(i)*"/force/fy")
        d = h5read(fn_h5, string(i)*"/damage")

        # Write point data - displacements
        write(fid, "POINT_DATA $num_nodes\n")
        write(fid, "VECTORS displacements float\n")
        write(fid, join(string.(ux), " ") * "\n")
        write(fid, join(string.(uy), " ") * "\n")
        #write(fid, join(string.(zeros(num_nodes)), " ") * "\n")

        # Write point data - node forces
        write(fid, "VECTORS force float\n")
        write(fid, join(string.(fx), " ") * "\n")
        write(fid, join(string.(fy), " ") * "\n")

        # Write point data - damage
        write(fid, "SCALARS damage float 1\n")
        write(fid, "LOOKUP_TABLE default\n")
        write(fid, join(string.(d), " ") * "\n")
    end
    # println("Binary VTK file successfully generated")
end
=#

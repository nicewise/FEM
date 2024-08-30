function load_fem_data(jld2fn::String, i::Int)
    stringi = string(i)
    JLD2.jldopen(jld2fn, "r") do file
        u = file[stringi]["displacement"]
        f = file[stringi]["force"]
        σ = file[stringi]["stress"]
        ε = file[stringi]["strain"]
        statev = file[stringi]["statev"]
        return u, f, σ, ε, statev
    end
end
function extrapolate_data(rc::fem_recorder, mesh::AbstractMesh, data)
    itp = ScatteredInterpolation.interpolate(rc.interpolation_method, rc.points, vcat(data'...))
    return ScatteredInterpolation.evaluate(itp, mesh.nodes)'
end
function plot2vtu(rc::fem_recorder, mesh::AbstractMesh, i::Int, e...)
    jld2fn = joinpath(rc.output_dir, rc.fn*".jld2")
    u, f, σ, ε, statev = load_fem_data(jld2fn, i)
    nodes = mesh.nodes
    dim, node_number = size(nodes)
    u = reshape(u, :, node_number)
    f = reshape(f, :, node_number)
    nodes = vcat(nodes, zeros(3 - dim, node_number))
    u = vcat(u, zeros(3 - dim, node_number))
    f = vcat(f, zeros(3 - dim, node_number))
    σ = extrapolate_data(rc, mesh, σ)
    ε = extrapolate_data(rc, mesh, ε)
    pointdata = isnothing(rc.statev_index) ?
        ((u, "displacement"), (f, "force"), (σ, "stress"), (ε, "strain")) :
        ((u, "displacement"), (f, "force"), (σ, "stress"), (ε, "strain"), (extrapolate_data(rc, mesh, statev), "statev"))
    vtkfn = joinpath(rc.vtk_dir, rc.fn * "_" * lpad(i, 4, '0') * ".vtu")
    plot2vtu(nodes, e..., pointdata = pointdata, fn = vtkfn)
end
function plot2vtu(rc::fem_recorder, mesh::AbstractMesh, i::Union{StepRange{Int64, Int64}, UnitRange{Int64}, Vector{Int}}, e...)
    for j in i # FIXME 这里暂时不能用@threads，有时间优化一下
        plot2vtu(rc, mesh, j, e...)
    end
end
plot2vtu(rc::fem_recorder, mesh::AbstractMesh, e...) = begin
    vtkfn = joinpath(rc.vtk_dir, rc.fn*"_mesh.vtu")
    plot2vtu(mesh.nodes, e...; pointdata = (), fn = vtkfn)
end

const celltype = Dict(
    # (dimension, vertex_number) => type
    (2, 3) => 5,  # t1
    (2, 6) => 22, # t2
    (2, 4) => 9,  # q1
    (2, 8) => 23, # q2
    (3, 4) => 10, # T1
    (3,10) => 24, # T2
    (3, 8) => 12, # H1
    (3,20) => 25, # H2
    (3, 6) => 13, # P1
    (3,15) => 26  # P2
)
function plot2vtu(n::Matrix{Float64}, e...; pointdata, fn)
    connectivity = vcat([vec(ele) for ele in e]...)
    offset = Int[]
    type = Int[]
    for ei in e
        (n_per_e, n_e) = size(ei)
        append!(offset, fill(n_per_e, n_e))
        e_type = celltype[(size(n, 1), size(ei, 1))]
        append!(type, fill(e_type, n_e))
    end
    plot2vtu(connectivity .- 1, cumsum(offset), type, n, pointdata = pointdata, fn = fn)
end

function plot2vtu(connectivity::Vector{Int}, offset::Vector{Int}, type::Vector{Int}, nodes::Matrix{Float64}; pointdata, fn)
    node_number = size(nodes, 2)
    cell_number = length(offset)
    open(fn, "w") do io
        write(io, """
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1">
  <UnstructuredGrid>
    <Piece NumberOfPoints="$node_number" NumberOfCells="$cell_number">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
""")
        write( io, "          " * join(nodes, " ") * "\n" )
        write(io, """
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
""")
        write( io, "          " * join(connectivity, " ") * "\n" )
        write(io, """
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
""")
        write( io, "          " * join(offset, " ") * "\n" )
        write(io, """
        </DataArray>
        <DataArray type="Int32" Name="types" format="ascii">
""")
        write( io, "          " * join(type, " ") * "\n" )
        write(io, """
        </DataArray>
      </Cells>
""")
        if !isempty(pointdata)
            write(io, """
      <PointData>
""")
            for (data, name) in pointdata
                NumberOfComponents = size(data, 1)
                write(io, """
        <DataArray type="Float32" Name="$name" NumberOfComponents="$NumberOfComponents" format="ascii">
""")
                write( io, "          " * join(data, " ") * "\n" )
                write(io, """
        </DataArray>
""")
            end
            write(io, """
      </PointData>
""")
        end
        write(io, """
    </Piece>
  </UnstructuredGrid>
</VTKFile>
""")
    end
    return nothing
end

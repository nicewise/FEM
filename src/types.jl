"""
定义了：
- 三角形： triangle
- 四边形： quadrilateral
- 四面体： tetrahedron
- 六面体： hexahedron
"""
abstract type AbstractCellType end
abstract type AbstractBeam <: AbstractCellType end
abstract type AbstractBond <: AbstractCellType end
abstract type Cell2 <: AbstractCellType end
abstract type Cell3 <: AbstractCellType end
abstract type triangular <: Cell2 end
abstract type quadrilateral <: Cell2 end
abstract type Tetrahedral <: Cell3 end
abstract type Hexahedral <: Cell3 end

abstract type t1 <: triangular end
abstract type t2 <: triangular end
abstract type q1 <: quadrilateral end
abstract type T1 <: Tetrahedral end
abstract type H1 <: Hexahedral end


const VertexNumber = Dict(
    t1 => 3,
    t2 => 6,
    q1 => 4,
    T1 => 4,
    H1 => 8
)
function assert_vertex_number(T, e)
    # 验证单元类型对应的节点数
    node_per_element = size(e, 1)
    expected_nodes = VertexNumber[T]
    @assert node_per_element == expected_nodes "For cell type '$T', expected $expected_nodes nodes per element, but found $node_per_element."
    return node_per_element
end

abstract type AbstractDimension end
abstract type Dim1 <: AbstractDimension end
abstract type Dim2 <: AbstractDimension end
abstract type Dim3 <: AbstractDimension end
abstract type PlaneStress <: Dim2 end
abstract type PlaneStrain <: Dim2 end

abstract type AbstractConstitutiveLaw{D<:AbstractDimension} end
abstract type AbstractPoint{C<:AbstractCellType, D<:AbstractDimension} end
abstract type AbstractCell{C<:AbstractCellType, D<:AbstractDimension} end
abstract type AbstractFem{C<:AbstractCellType, D<:AbstractDimension} end

abstract type AbstractSolver end
# default solver wrapper
struct LeftDivisionSolver <: AbstractSolver end
solveit(solver::LeftDivisionSolver, A, b) =  A \ b

abstract type AbstractAnalysis end
abstract type AbstractMesh{D<:AbstractDimension} end

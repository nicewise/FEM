module FEM
using LinearAlgebra, SparseArrays, StaticArrays
using Base.Threads

include("types.jl")
export
    AbstractCellType,
    Cell2, Cell3,
    t1, t2,
    q1,
    T1,
    H1,
    VertexNumber,
    assert_vertex_number,
    AbstractDimension,
    Dim1, Dim2, Dim3,
    PlaneStress, PlaneStrain,
    AbstractConstitutiveLaw,
    AbstractQuadraturePoint,
    AbstractElement,
    AbstractFem,
    AbstractSolver,
    solveit,
    LeftDivisionSolver,
    AbstractAnalysis,
    AbstractMesh

include("MandelNotation.jl")
export  Mandel # some constant tensor

include("mesh.jl")
export
    fem_mesh,
    fem_mesh_gen

include("shape.jl")
export
    B_gen,
    N_gen

include("quadrature.jl")
export
    quad_form

include("linear_elastic_fem.jl")
export
    linear_elastic, # <: AbstractFem
    stiffness_assemble

include("general_elastoplastic_fem.jl")
export
    quadrature_point, # <: AbstractPoint
    simple_point_gen, # for constitutive law testing
    elastoplastic, # <: AbstractFem
    constitutive_law_apply!,
    elastic_initialization!,
    elements_stiffness_gen!,
    displacement_apply!

include("popular_constitutive_laws.jl")
export
    constitutive_linear_elastic

include("solvers.jl")

include("analysis_procedures.jl")
export
    simple_analysis,
    displacement_load,
    displacement_converge

end # module FEM

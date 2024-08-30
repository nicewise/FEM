module FEM
using LinearAlgebra, SparseArrays, StaticArrays
using Base.Threads
import ScatteredInterpolation, JLD2, Roots

include("types.jl")
export
    AbstractCellType,
    AbstractBeam,
    AbstractBond,
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
    AbstractPoint,
    AbstractCell,
    AbstractFem,
    AbstractSolver,
    solveit,
    LeftDivisionSolver,
    AbstractMesh,
    AbstractAnalysis,
    AbstractRecorder

include("MandelNotation.jl")
export  Mandel # some constant tensor

include("mesh.jl")
export
    fem_mesh,
    fem_mesh_gen

include("shape.jl")
export
    N_gen,
    B_gen,
    elastic_B_gen

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
    elastic_initialization!,
    displacement_apply!,
    constitutive_law_apply!,
    stiffness_assemble!,
    internal_force_assemble!,
    elastic_system_assemble!

include("gradient_enhanced_fem.jl")
export
    gradient_stiffness_assemble,
    gradient_rhs_gen,
    gradient_apply!

include("popular_constitutive_laws.jl")
export
    constitutive_linear_elastic,
    AbstractRd,
    Rd_gen,
    Rd,
    elastic_damage,
    gradient_elastic_damage,
    before_gradient!,
    after_gradient!

include("solvers.jl")

include("analysis_procedures.jl")
export
    simple_analysis,
    displacement_load,
    displacement_converge,
    record!,
    testcon, # <: AbstractRecorder
    testcon!,
    fem_recorder, # <: AbstractRecorder
    displacement_control!,
    gradient_displacement_control!

include("painters.jl")
export plot2vtu

end # module FEM

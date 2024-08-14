module FEM
using LinearAlgebra, SparseArrays, HDF5
using Base.Threads

include("mesh.jl")
export
    mesh,
    mesh_gen,
    save_model_h5,
    vtk_gen_from_h5_binary,
    vtk_gen_from_h5_ascii,
    vtk_gen_from_h5_ascii2

include("types.jl")
export
    FE,
    fe2d, fe3d,
    fet1, fet2, fet3,
    feq1, feq2, feq3,
    feT1, feT2, feT3,
    feH1, feH2, feH3,
    feP1, feP2, feP3,
    AbstractConstitutiveLaw,
    AbstractQuadraturePoint,
    AbstractElement,
    AbstractFem

include("shape.jl")
export
    B_gen,
    N_gen

include("quadrature.jl")
export
    quad_form

include("simple_fem.jl")
export
    simple_fem_gen,
    stiffness_assemble

include("general_elastoplastic_fem.jl")
export
    fem_gen,
    elastic_initialization!,
    elements_stiffness_gen!

include("popular_constitutive_laws.jl")
export
    constitutive_law_apply!,
    linear_elastic

include("solve.jl")
export
    implicit_static_analysis,
    load_attempt,
    R_gen!,
    save_step_data_h5

end # module FEM

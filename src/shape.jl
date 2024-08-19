"""
N_r --- the derivative of N to the reference coordinates
N_g --- the derivative of N to the global coordinates
JxW --- the product of the determinant of the Jacobian matrix and the weight of quadrature point
"""
function B_gen(C::Type, D::Type, ref_coords::Vector{Float64}, weight::Float64, coords::Matrix{Float64})
    N_r = ∂N_gen(C, ref_coords)
    N_g, JxW = B_gen(N_r, weight, coords)
    dim = D <: Dim2 ? 2 : 3
    dof = dim * size(coords, 2)
    B = dim == 2 ? zeros(3, dof) : zeros(6, dof)
    if dim == 2
        for i in axes(N_g, 1)
            B[1, 2i-1] = B[3,   2i] = N_g[i, 1]
            B[2,   2i] = B[3, 2i-1] = N_g[i, 2]
        end
    else
        for i in axes(N_g, 1)
            B[1, 3i-2] = B[5,   3i] = B[6, 3i-1] = N_g[i, 1]
            B[2, 3i-1] = B[4,   3i] = B[6, 3i-2] = N_g[i, 2]
            B[3,   3i] = B[4, 3i-1] = B[5, 3i-2] = N_g[i, 3]
        end
    end
    B[dim+1:end, :] ./= sqrt(2)
    return B, JxW
end
@inline function B_gen(N_r::Matrix{Float64}, weight::Float64, coords::Matrix{Float64})
    J = coords * N_r
    N_g = N_r / J
    return N_g, weight * det(J)
end


N_gen(::Type{t1}, ref_coords::Vector{Float64}) = ref_coords
∂N_gen(::Type{t1}, ref_coords::Vector{Float64}) = [-1 -1
                                                   1  0
                                                   0  1]

function N_gen(::Type{t2}, ref_coords::Vector{Float64})
    L1, L2, L3 = ref_coords
    [(2L1 - 1) * L1
     (2L2 - 1) * L2
     (2L3 - 1) * L3
     4L1*L2
     4L2*L3
     4L3*L1]
end
@inline function ∂N_gen(::Type{t2}, ref_coords::Vector{Float64})
    L1, L2, L3 = ref_coords
    [ 1-4L1     1-4L1
      4L2-1     0
      0         4L3-1
      4L1-4L2  -4L2
      4L3       4L2
     -4L3       4L1-4L3]
end

N_gen(::Type{T1}, ref_coords::Vector{Float64}) = ref_coords
∂N_gen(::Type{T1}, ref_coords::Vector{Float64}) = [-1 -1 -1
                                                   1  0  0
                                                   0  1  0
                                                   0  0  1]

function N_gen(::Type{q1}, ref_coords::Vector{Float64})
    ξ, η = ref_coords
    [(1-ξ)*(1-η)
     (1+ξ)*(1-η)
     (1+ξ)*(1+η)
     (1-ξ)*(1+η)] / 4
end
@inline function ∂N_gen(::Type{q1}, ref_coords::Vector{Float64})
    ξ, η = ref_coords
    [ η-1  ξ-1
      1-η -1-ξ
      1+η  1+ξ
     -1-η  1-ξ] / 4
end

function N_gen(::Type{H1}, ref_coords::Vector{Float64})
    ξ, η, ζ = ref_coords
    [(1-ξ)*(1-η)*(1-ζ)
     (1+ξ)*(1-η)*(1-ζ)
     (1+ξ)*(1+η)*(1-ζ)
     (1-ξ)*(1+η)*(1-ζ)
     (1-ξ)*(1-η)*(1+ζ)
     (1+ξ)*(1-η)*(1+ζ)
     (1+ξ)*(1+η)*(1+ζ)
     (1-ξ)*(1+η)*(1+ζ)] / 8
end
@inline function ∂N_gen(::Type{H1}, ref_coords::Vector{Float64})
    ξ, η, ζ = ref_coords
    [-(1-η)*(1-ζ) -(1-ξ)*(1-ζ) -(1-ξ)*(1-η)
      (1-η)*(1-ζ) -(1+ξ)*(1-ζ) -(1+ξ)*(1-η)
      (1+η)*(1-ζ)  (1+ξ)*(1-ζ) -(1+ξ)*(1+η)
     -(1+η)*(1-ζ)  (1-ξ)*(1-ζ) -(1-ξ)*(1+η)
     -(1-η)*(1+ζ) -(1-ξ)*(1+ζ)  (1-ξ)*(1-η)
      (1-η)*(1+ζ) -(1+ξ)*(1+ζ)  (1+ξ)*(1-η)
      (1+η)*(1+ζ)  (1+ξ)*(1+ζ)  (1+ξ)*(1+η)
     -(1+η)*(1+ζ)  (1-ξ)*(1+ζ)  (1-ξ)*(1+η)] / 8
end

function B_gen(::Type{feq1}, ref_coords::Vector{Float64}, weight::Float64, coords)
    ξ, η = ref_coords
    N_xieta = [η-1  1-η 1+η -1-η
               ξ-1 -1-ξ 1+ξ  1-ξ] / 4
    J = N_xieta * coords' #FIXME
    N_xy = J \ N_xieta
    B = [N_xy[1,1] 0         N_xy[1,2] 0         N_xy[1,3] 0         N_xy[1,4] 0
         0         N_xy[2,1] 0         N_xy[2,2] 0         N_xy[2,3] 0         N_xy[2,4]
         N_xy[2,1] N_xy[1,1] N_xy[2,2] N_xy[1,2] N_xy[2,3] N_xy[1,3] N_xy[2,4] N_xy[1,4]]
    JxW = weight * det(J)
    return B, JxW
end
function N_gen(::Type{feq1}, ref_coords::Vector{Float64})
    ξ, η = ref_coords
    N1 = (1-ξ)*(1-η)/4
    N2 = (1+ξ)*(1-η)/4
    N3 = (1+ξ)*(1+η)/4
    N4 = (1-ξ)*(1+η)/4
    return [N1, N2, N3, N4]
end

function B_gen(::Type{feH1}, ref_coords::Vector{Float64}, weight::Float64, coords)
    (fe::fedata{feH1}, coors::Matrix{Float64}, gi::Int)
    ξ, η, ζ = ref_coords
    N_local = [-(1-η)*(1-ζ)  (1-η)*(1-ζ)  (1+η)*(1-ζ) -(1+η)*(1-ζ) -(1-η)*(1+ζ)  (1-η)*(1+ζ) (1+η)*(1+ζ) -(1+η)*(1+ζ)
               -(1-ξ)*(1-ζ) -(1+ξ)*(1-ζ)  (1+ξ)*(1-ζ)  (1-ξ)*(1-ζ) -(1-ξ)*(1+ζ) -(1+ξ)*(1+ζ) (1+ξ)*(1+ζ)  (1-ξ)*(1+ζ)
               -(1-ξ)*(1-η) -(1+ξ)*(1-η) -(1+ξ)*(1+η) -(1-ξ)*(1+η)  (1-ξ)*(1-η)  (1+ξ)*(1-η) (1+ξ)*(1+η)  (1-ξ)*(1+η)] / 8
    J = N_local * coors # 3x8 x 8x3 = 3x3
    N_global = J \ N_local # 3x8
    B = 
    return(N_global, det(J)*weight)
end
function N_gen(::Type{feH1}, ref_coords::Vector{Float64})
    ξ, η, ζ = ref_coords
    N1 = (1-ξ)*(1-η)*(1-ζ)/8
    N2 = (1+ξ)*(1-η)*(1-ζ)/8
    N3 = (1+ξ)*(1+η)*(1-ζ)/8
    N4 = (1-ξ)*(1+η)*(1-ζ)/8
    N5 = (1-ξ)*(1-η)*(1+ζ)/8
    N6 = (1+ξ)*(1-η)*(1+ζ)/8
    N7 = (1+ξ)*(1+η)*(1+ζ)/8
    N8 = (1-ξ)*(1+η)*(1+ζ)/8
    return [N1, N2, N3, N4, N5, N6, N7, N8]
end

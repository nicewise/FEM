constitutive_law_apply!(qpoint::AbstractQuadraturePoint{S}, F::AbstractConstitutiveLaw{T}) where {T<:FE, S<:T} = F(qpoint)
function constitutive_law_apply!(fe::AbstractFem{S}, F::AbstractConstitutiveLaw{T}) where {T<:FE, S<:T}
    @threads for qpoint in fe.quadrature_points
        F(qpoint)
    end
    return nothing
end

struct linear_elastic{T<:FE} <: AbstractConstitutiveLaw{T}
    E::Float64
    ν::Float64
    D::Matrix{Float64}
end
function linear_elastic{T}(E, ν; planestrain::Bool = false) where T<:FE
    D = 0;
    if T <: fe3d
        D = [1-ν ν   ν   0    0    0
             ν   1-ν ν   0    0    0
             ν   ν   1-ν 0    0    0
             0   0   0   1-2ν 0    0
             0   0   0   0    1-2ν 0
             0   0   0   0    0    1-2ν]*E/(1+ν)/(1-2ν)
    elseif planestrain
        D = [1-ν ν   0
             ν   1-ν 0
             0   0   1-2ν]*E/(1+ν)/(1-2ν)
    else
         D = [1 ν 0
              ν 1 0
              0 0 (1-ν)/2]*E/(1-ν^2)
    end
    linear_elastic{T}(E, ν, D)
end
function (F::linear_elastic{T})(qpoint::quadrature_point{S}) where {T<:FE, S<:T}
    ε = qpoint.ε + qpoint.dε
    qpoint.σ .= F.D * ε
    qpoint.D .= F.D
end

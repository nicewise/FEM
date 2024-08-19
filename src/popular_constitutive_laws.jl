struct constitutive_linear_elastic{T<:AbstractDimension} <: AbstractConstitutiveLaw{T}
    E::Float64
    ν::Float64
    D::Matrix{Float64}

    function constitutive_linear_elastic{T}(E, ν) where T
        D = 0
        if T <: Dim3
            D = [1-ν ν   ν   0    0    0
                 ν   1-ν ν   0    0    0
                 ν   ν   1-ν 0    0    0
                 0   0   0   1-2ν 0    0
                 0   0   0   0    1-2ν 0
                 0   0   0   0    0    1-2ν]*E/(1+ν)/(1-2ν)
        elseif T <: PlaneStrain
            D = [1-ν ν   0
                 ν   1-ν 0
                 0   0   1-2ν]*E/(1+ν)/(1-2ν)
        else
            D = [1 ν 0
                 ν 1 0
                 0 0 (1-ν)/2]*E/(1-ν^2) # FIXME
        end
        return new{T}(E, ν, D)
    end
end

function constitutive_law_apply!(F::constitutive_linear_elastic, p::AbstractPoint)
    ε = p.ε + p.dε
    p.σ .= F.D * ε
    p.D .= F.D
end

struct vm_iso{T<:Dim3} <: AbstractConstitutiveLaw{T}
    E::Float64
    ν::Float64
    σy0::Float64
    H::Float64
    tol::Float64
    G::Float64
    G3::Float64
    D::Matrix{Float64}

    function vm_iso{T}(E, ν, σy0, H) where T
        k = E/3/(1-2ν) # 体积模量
        G = E/2/(1+ν)  # 剪切模量
        D = 3k * Mandel.J + 2G * Mandel.K
        return new{T}(E, ν, σy0, H, σy0*1e-6, G, 3G, D)
    end
end

function constitutive_law_apply!(F::vm_iso{T}, p::AbstractPoint{<:AbstractCellType, T}) where T<:Dim3
    εp = p.statev[1:6] # plastic strain
    γp = p.statev[7] # effective plastic strain

    σ = F.D * (p.ε + p.dε - εp) # trial stress
    s = Mandel.K * σ
    J2 = sum(σ.^2)
    q = sqrt(3/2) * J2
    f = q - (F.σy0 + F.H*γp)
    function yield!()
        Δγ = f / (F.G3 + F.H)
        N = s / J2 # 6x1
        Nb = sqrt(2/3) * N # \bar{N}
        σ -= F.G3 * Δγ * Nb
        εp += Δγ * N
        γp += Δγ
        c1 = 6*F.G^2/(F.G3+H)
        c2 = 6*G^2*Δγ/q
        D = F.D - (c1-c2)*Nb*Nb' - c2*Mandel.K
        #D = F.D - c1 * Nb*Nb'
        return D
    end
    D = f >= F.tol ? yield!() : F.D
    p.statev[1:6] .= εp
    p.statev[7] = γp
    p.D .= D
    p.σ .= σ
    return nothing
end

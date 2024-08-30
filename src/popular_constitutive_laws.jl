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
    return nothing
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

abstract type AbstractRd end
struct R1 <: AbstractRd
    rc::Float64
    dc::Float64
end
Rd(R::R1, d) = 4d*R.dc*R.rc/(d + R.dc)^2

function Rd_gen(filter::Int, p::Float64...)
    R_map = Dict(
        1 => R1(p[1], p[2])
        # TODO
    )
    return R_map[filter]
end

#FIXME delute跑不出正确的本构曲线
# 对于pcw，体积模量和剪切模量不劣化到小于0已经解决了问题
# 但是delute不知道问题在哪
struct elastic_damage{T<:Dim3} <: AbstractConstitutiveLaw{T}
    a::Float64
    x::Int
    y::Int
    k3::Float64
    μ2::Float64
    weaken_ratio::Dict{Bool, NTuple{4, Function}}
    isopen::Function
    Rd::AbstractRd
    tol::Float64

    function elastic_damage{T}(S::Symbol, E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, l::Float64 = 0) where T
        β1 = 16(1 - ν^2) / 9(1 - 2ν)
        β2 = 32(1 - ν) / 15(2 - ν)
        ϑ = (5 - ν) / 3
        weaken_ratio = Dict(
            :DELUTE => Dict(
                true => (
                    d -> max(1e-10, 1 - β1 * d),
                    d -> - β1,

                    d -> max(1e-10, 1 - β2 * ϑ * d),
                    d -> - β2 * ϑ
                ),
                false => (
                    d -> 1,
                    d -> 0,

                    d -> max(1e-10, 1 - β2 * d),
                    d -> - β2
                )
            ),
            :MT => Dict(
                true => (
                    d -> 1 / (1 + β1 * d),
                    d -> - β1 / (1 + β1 * d)^2,

                    d -> 1 / (1 + β2 * ϑ * d),
                    d -> - β2 * ϑ / (1 + β2 * ϑ * d)^2
                ),
                false => (
                    d -> 1,
                    d -> 0,

                    d -> 1 / (1 + β2 * d),
                    d -> - β2 / (1 + β2 * d)^2
                )
            ),
            :PCW => Dict(
                true => (
                    d -> max(1e-10, 1 - 16(1 - ν^2)d / ( 9(1 - 2ν) + 16/3*(1 + ν)^2 * d )),
                    d -> 1296(1 - 2ν)*(ν^2 - 1)/(16d*(ν + 1)^2 - 54ν + 27)^2,

                    d -> max(1e-10, 1 - 480(1 - ν)ϑ * d / ( 225(2 - ν) + 64(4 - 5ν)ϑ * d )),
                    d -> 324000(ν - 5)*(ν - 2)*(ν - 1)/(64d*(ν - 5)*(5ν - 4) - 675ν + 1350)^2
                ),
                false => (
                    d -> 1,
                    d -> 0,

                    d -> max(1e-10, 1 - 480(1 - ν)d / ( 225(2 - ν) + 64(4 - 5ν)d )),
                    d -> -108000(ν - 2)*(ν - 1)/(320d*ν - 256d + 225ν - 450)^2
                )
            )
        )

        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        D = k3 * Mandel.J + μ2 * Mandel.K
        isopen = x -> sum(x .* Mandel.δ) / 3 >= 0
        return new{T}(l, 1, 2, k3, μ2, weaken_ratio[S], isopen, Rd, tol)
    end
end

function constitutive_law_apply!(F::elastic_damage{T}, p::AbstractPoint{<:AbstractCellType, T}; before_gradient::Union{Nothing, Bool}=nothing) where T<:Dim3
    if isnothing(before_gradient)
        d = p.statev[1]
        ε = p.ε + p.dε
        kd, dkd, μd, dμd = F.weaken_ratio[F.isopen(p.D * ε)]
        k3 = F.k3
        μ2 = F.μ2
        g(x) = -ε' * (dkd(d+x) * k3 * Mandel.J + dμd(d+x) * μ2 * Mandel.K) * ε / 2 - Rd(F.Rd, d+x)
        if g(0.) > F.tol
            d += Roots.find_zero(g, 0.0)
        end
        Chom = kd(d) * k3 * Mandel.J + μd(d) * μ2 * Mandel.K
        p.D .= Chom
        p.σ .= Chom * ε
        p.statev[1] = d
        return nothing
    elseif before_gradient
        d = p.statev[1]
        ε = p.ε + p.dε
        _, dkd, _, dμd = F.weaken_ratio[F.isopen(p.D * ε)]
        k3 = F.k3
        μ2 = F.μ2
        g(x) = -ε' * (dkd(d+x) * k3 * Mandel.J + dμd(d+x) * μ2 * Mandel.K) * ε / 2 - Rd(F.Rd, d+x)
        if g(0.) > F.tol
            d += Roots.find_zero(g, 0.0)
        end
        p.statev[1] = d
        return nothing
    else
        d = p.statev[2]
        ε = p.ε + p.dε
        kd, _, μd, _ = F.weaken_ratio[F.isopen(p.D * ε)]
        Chom = kd(d) * k3 * Mandel.J + μd(d) * μ2 * Mandel.K
        p.D .= Chom
        p.σ .= Chom * ε
        return nothing
    end
end



struct PDMT end
struct PDPCW end
struct AEDPCW end




struct APDPCW end

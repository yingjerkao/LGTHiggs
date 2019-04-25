# code for Higgs fields
const AdjointHiggsVar = Vector{Float64}
export AdjointHiggsVar


struct AdjointHiggsField{d} <: AbstractField{d}
    field::Array{AdjointHiggsVar, d}
    lattice::AbstractLattice{d}

    AdjointHiggsField{d}(st::AbstractArray{AdjointHiggsVar,d}, lat::AbstractLattice) where d = new(st, lat)
end
export AdjointHiggsField

function AdjointHiggsField(st::AbstractArray{AdjointHiggsVar,d}, lat::AbstractLattice{d}) where d
    AdjointHiggsField{d}(st, lat)
end

# for now this just creates ϕ=(0,0,1) everywhere (Unitary Gauge)
function AdjointHiggsField(lat::AbstractLattice, n::Integer=3)
    st=fieldsites(AdjointHiggsVar, lat)
    for idx ∈ eachindex(st)
        st[idx] = [ (0.0 for i in 1:(n-1))... , 1.0]
    end
    AdjointHiggsField(st, lat)
end
function (ϕ::AdjointHiggsField)(x::Tuple)
    idx = getindex(ϕ.lattice, x)
    getindex(ϕ.field,idx...)
end
(ϕ::AdjointHiggsField)(x::CartesianIndex) = ϕ(Tuple(x))
(ϕ::AdjointHiggsField)(x...) = ϕ(Tuple(x))



function AdjointLink(U::GaugeFieldSU2, ν:: Int, x::CartesianIndex)

    Alink=zeros(ComplexF64,3,3)
    for α in 1:3
        for β in 1:α
            Alink[α,β]=0.5*tr(U(x)[ν]*σ[α]*U(x)[ν]'*σ[β])
            Alink[β,α]=Alink[α,β]
        end
    end
    return Alink
end
export AdjointLink

function higgslagrangian(U::AbstractGaugeField{d}, ϕ::AbstractField{d},βh::
    Float64, x::CartesianIndex) where d


    ℒ = 0.0
    for ν in 1:d
        Alink=AdjointLink(U,ν,x)
        ℒ +=ϕ(x)'* Alink *ϕ(shift(x,ν))
    end
    βh*ℒ
end
export higgslagrangian
function higgsaction(U::AbstractGaugeField{d}, ϕ::AbstractField{d},βh::Float64) where d
    S = 0.0
    for x ∈ eachspacetimeindex(U)
        S += higgslagrangian(U, ϕ, βh, x)
    end
    S/U.lattice.a^(4 - d) # this factor to fix the powers of a
end
export higgsaction

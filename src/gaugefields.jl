# code for gauge fields

# TODO for now the wilson lines will be `Vector`s of SU(2) matrices
# TODO looks like for arbitrary dimensions we'll have lots of splatting inefficiency...
#   this should be ok because most of this stuff should just be for initialization...
#   in the future can fix splatting inefficiency with macros

const WilsonLineVar = Vector{UnitaryMatrix}
export WilsonLineVar

wilsonlineunity(N::Integer, d::Integer) = UnitaryMatrix[UnitaryMatrix(I, N,N) for i ∈ 1:d]
wilsonlineunity(N::Integer, lat::AbstractLattice{d}) where d = wilsonlineunity(N, d)
export wilsonlineunity

randwilsonlineSU2(d::Integer) = UnitaryMatrix[randSU2() for i ∈ 1:d]
randwilsonlineSU2(lat::AbstractLattice{d}) where d = randwilsonlineSU2(d)
export randwilsonlineSU2

"""
    GaugeFieldSU2{d} <: AbstractGaugeField

An SU(2) gauge field.  Note that we take the fundamental field variables to be Wilson lines emanating
from lattice site x in direction μ. Here `d` is the total number of dimensions `D+1`.
"""
struct GaugeFieldSU2{d} <: AbstractGaugeField{d}
    field::Array{WilsonLineVar, d}
    lattice::AbstractLattice{d}

    GaugeFieldSU2{d}(st::AbstractArray{WilsonLineVar,d}, lat::AbstractLattice) where d = new(st, lat)
end
export GaugeFieldSU2

function GaugeFieldSU2(st::AbstractArray{WilsonLineVar,d}, lat::AbstractLattice{d}) where d
    GaugeFieldSU2{d}(st, lat)
end

# for now this just creates a "vacuum" with A=0 everywhere
function GaugeFieldSU2(lat::AbstractLattice)
    # note inefficiency due to splatting
    st = fieldsites(WilsonLineVar, lat)
    for idx ∈ eachindex(st)
        st[idx] = wilsonlineunity(2, lat)
    end
    GaugeFieldSU2(st, lat)
end


# TODO ensure that splatting is really ok here!
# gets the gauge field at spacetime location x1, t
function (U::GaugeFieldSU2)(x::Tuple)
    idx = getindex(U.lattice, x)
    getindex(U.field, idx...)
end
(U::GaugeFieldSU2)(x::CartesianIndex) = U(Tuple(x))
(U::GaugeFieldSU2)(x...) = U(Tuple(x))

# returns spacetime function of U[μ]
getindex(U::GaugeFieldSU2, μ::Integer) = x -> U(x)[μ]

eachspacetimeindex(U::AbstractGaugeField{d}) where d = eachspacetimeindex(U.lattice)

function rand(::Type{GaugeFieldSU2}, lat::AbstractLattice)
    st= fieldsites(WilsonLineVar, lat)
    for idx ∈ eachindex(st)
        st[idx] = randwilsonlineSU2(lat)
    end
    GaugeFieldSU2(st, lat)
end


"""
    untracedplaquette(U::AbstractGaugeField, μ::Int, ν::Int, x::Tuple)

Compute the value of a plaquette before taking the trace.  In other words
\$\$ U_{\\mu}(x)U_{\\nu}(x+e_{\\mu})U^{\\dagger}_{\\mu}(x+e_{\\nu})U^{\\dagger}_{\\nu}(x) \$\$
"""
function untracedplaquette(U::AbstractGaugeField, μ::Integer, ν::Integer, x::CartesianIndex)
    U[μ](x)*U[ν](shift(x,μ))*U[μ](shift(x,ν))'*U[ν](x)'
end
export untracedplaquette

"""
    plaquette(U::AbstractGaugeField, μ::Integer, ν::Integer, x::CartesianIndex)

Compute the value of the plaquette at location `x` over the gauge field.

Note that we follow the definition in the "Lattice QCD For Novices" paper and define the paquette
to be the real part of the trace.

*TODO*: Note that ideally the normalization factor here should depend on the gauge group.
"""
plaquette(U::AbstractGaugeField, μ::Integer, ν::Integer, x::CartesianIndex) = real(tr(untracedplaquette(U,μ,ν,x)))/2.0
export plaquette
"""
    staple(U::AbstractGaugeField, μ::Integer, x::CartesianIndex)

Compute the staple connecting the link U[μ](x)

"""
function staple(U::AbstractGaugeField{d}, μ::Integer, x::CartesianIndex) where d
    st=UnitaryMatrix(zeros(2,2))
    for ν in  1:d
        if ν != μ
            st+=(U[ν](shift(x,μ))*U[μ](shift(x,ν))'*U[ν](x)'
            +U[ν](shift(shift(x,μ),-ν))'*U[μ](shift(x,-ν))'*U[ν](shift(x,μ)))
        end
    end
    return st
end
staple(U::AbstractGaugeField{d}, μ::Integer, x::NTuple{d, Int}) where d =
    staple(U, μ, CartesianIndex(x))

export staple
"""
    wilsonlagrangian(U::AbstractGaugeField, β::AbstractFloat, x::CartesianIndex)

Compute the lowest order Wilson plaquette Lagrangian at location `x` and inverse temperature `β`.

Note that this has a factor of a^4 in front of it relative to the usual G^2 Yang-Mills kinetic term.

Note that we neglect the constant term altogether.
"""
function wilsonlagrangian(U::AbstractGaugeField{d}, β::AbstractFloat, x::CartesianIndex) where d
    ℒ = 0.0
    for μ ∈ 1:d
        for ν ∈ 1:(μ-1)
            ℒ -= plaquette(U, μ, ν, x)
        end
    end
    β*ℒ
end
export wilsonlagrangian


"""
    wilsonaction(U::AbstractGaugeField, β::AbstractFloat)

Compute the Wilson plaquette action at inverse temperature `β`.

Again, the constant term is neglected.
"""
function wilsonaction(U::AbstractGaugeField{d}, β::AbstractFloat) where d
    S = 0.0
    for x ∈ eachspacetimeindex(U)
        S += wilsonlagrangian(U, β, x)
    end
    S/U.lattice.a^(4 - d) # this factor to fix the powers of a
end
export wilsonaction

function renormalize_GaugeFieldSU2!(U::GaugeFieldSU2{d}) where d
    for x in eachspacetimeindex(U)
        for μ in 1:d
            project_SU2!(U(x)[μ])
        end
    end
end

export renormalize_GaugeFieldSU2!

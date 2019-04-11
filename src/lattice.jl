# Code to generate lattices

#Base.getindex(lat::AbstractLattice, idx::Tuple) = idx
#Base.getindex(lat::AbstractLattice, idx...) = getindex(lat, idx)



"""
    fieldsites(::Type{T}, lat::AbstractLattice)

Creates an uninitialized array with elements of type T for every point on the lattice.
"""
function fieldsites(::Type{T}, lat::AbstractLattice{d}) where {T, d}
    Array{T, d}(undef,(lat.L for i ∈ 1:(d-1))..., lat.T)
end
export fieldsites


nsites(l::AbstractLattice{d}) where d = l.T*(l.L^(d-1))
export nsites


# this is an implementation that doesn't involve splatting


size(lat::AbstractLattice{d}) where d = tuple((lat.L for i ∈ 1:(d-1))..., lat.T)

sizenocheck(lat::AbstractLattice{d}, idx::Integer) where d = (idx == d) ? lat.T : lat.L

"""
    eachspacetimeindex(lat::AbstractLattice)

Return an iterator over every spacetime index in the lattice `lat`.  Spacetime indices will be
in the form of `CartesianIndex`s.
"""
function eachspacetimeindex(lat::AbstractLattice{d}) where d
    CartesianIndices(size(lat))
end
export eachspacetimeindex

"""
    checkerboardspacetimeindex(lat::AbstractLattice)

Return a tuple of iterators over spacetime index with odd and even parity in
the lattice `lat`.  Spacetime indices will be in the form of `CartesianIndex`s.
"""
function checkerboardspacetimeindex(lat::AbstractLattice{d}) where d

    T=CartesianIndices(size(lat))
    Parity_Odd=[]
    Parity_Even=[]
    for i in T
        isodd(sum(i.I)) ?  push!(Parity_Odd, i) : push!(Parity_Even, i)
    end
    return Parity_Odd, Parity_Even
end
export checkerboardspacetimeindex


function size(lat::AbstractLattice{d}, idx::Integer) where d
    if idx == d
        return lat.T
    elseif 1 ≤ idx < d
        return lat.L
    else
        ArgumentError("Lattice only has $d spacetime dimensions.")
    end
end


"""
    LatticeToroidal <: AbstractLattice

A type for storing parameters of a simple Euclidean rectilinear lattice with D spatial dimensions
and toroidal topology in all directions.
"""
struct LatticeToroidal{d} <: AbstractLattice{d}
    a::Float64  # lattice spacing in GeV-1
    L::Int  # lattice length in units of a
    T::Int  # lattice duration in units of a
end
export LatticeToroidal

LatticeToroidal(d::Int, a::Float64, L::Integer, T::Integer) = LatticeToroidal{d}(a, L, T)
LatticeToroidal(d::Int, a::Float64, L::Integer) = LatticeToroidal(d, a, L, L)

# compute index
@inline function _loopedidx(i::Integer, L::Integer)
    if i < 1
        L + (i % L)
    else
        ((i-1) % L) + 1
    end
end

# these duplications are to avoid splatting
# WARNING right now these don't check lattice dimensions
# for rank-2
function Base.getindex(lat::LatticeToroidal, i::Integer, j::Integer)
    _loopedidx(i, lat.L), _loopedidx(j, lat.T)
end
# for rank-3
function Base.getindex(lat::LatticeToroidal, i::Integer, j::Integer, k::Integer)
    _loopedidx(i, lat.L), _loopedidx(j, lat.L), _loopedidx(k, lat.T)
end
# for rank-4
function Base.getindex(lat::LatticeToroidal, i::Integer, j::Integer, k::Integer, l::Integer)
    _loopedidx(i, lat.L), _loopedidx(j, lat.L), _loopedidx(k, lat.L), _loopedidx(l, lat.L)
end
Base.getindex(lat::LatticeToroidal, idx...) = getindex(lat, idx)
Base.getindex(lat::LatticeToroidal, idx::Tuple) = getindex(lat, idx...)
Base.getindex(lat::LatticeToroidal, idx::CartesianIndex) = getindex(lat, Tuple(idx))

# TODO compiler will not know size of these!
spacetimebasisvec(μ::Integer, N::Integer) = CartesianIndex{N}(Tuple(i ≠ μ ? 0 : 1 for i ∈ 1:N))
spacetimebasisvec(μ::Integer, lat::AbstractLattice{d}) where d = spacetimebasisvec(μ, d)
export spacetimebasisvec

function shift(x::CartesianIndex{N}, μ::Integer) where N
    _μ=abs(μ)
    @assert _μ<=N && _μ>0 "Index μ=$μ out of range"
    x+=sign(μ)*CartesianIndex{N}(Tuple(i == _μ ? 1 : 0 for i in 1:N))
end

export shift

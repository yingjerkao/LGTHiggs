__precompile__(true)
module LGTHiggs

# abstract types, d is number of spacetime dimensions
abstract type AbstractLattice{d} end
# the type parameters give the spacetime dimensionality
abstract type AbstractField{d} end
abstract type AbstractGaugeField{d} <: AbstractField{d} end
export AbstractLattice, AbstractField, AbstractGaugeField

import Base.rand
import Base: size, getindex, dec

include("unitary.jl")
include("lattice.jl")
include("gaugefields.jl")
include("metropolis.jl")
include("measurement.jl")
end # module

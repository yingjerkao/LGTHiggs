# code for generating uniform random unitary matrices
# TODO for now we'll just do SU(2)

using LinearAlgebra
const UnitaryMatrix = Matrix{Complex{Float64}}

#Pauli Matrices
const σ₁ = Complex{Float64}[0 1; 1 0]
const σ₂ = Complex{Float64}[0 -im; im 0]
const σ₃ = Complex{Float64}[1 0; 0 -1]
#Idensity Matrix
const one2by2 = UnitaryMatrix(I, 2, 2)

"""
    SU2(a::AbstractVector)
    Create a SU2 matrix from four points on S4
"""

function SU2(a::AbstractVector)
    @assert abs(norm(a)-1)<1e-6
    one2by2*a[1]+im*(a[2]*σ₁+a[3]*σ₂+a[4]*σ₃)
end

export SU2



"""
Marsaglia method to generate a random SU(2) matrix

G. Marsaglia, “Choosing a Point from the Surface of a Sphere,” Ann. Math. Statist. 43, 645 (1972).

"""
function randSU2()
    u0=1.
    u1=1.
    u2=1.
    u3=1.
    while(u0^2+u1^2>1)
        u0=rand()*2.0-1
        u1=rand()*2.0-1
    end
    while(u2^2+u3^2>1)
        u2=rand()*2.0-1
        u3=rand()*2.0-1
    end

    R=sqrt((1. -u0^2-u1^2)/(u2^2+u3^2))
    u2=u2*R
    u3=u3*R
    U=SU2([u0,u1,u2,u3])
end
export randSU2

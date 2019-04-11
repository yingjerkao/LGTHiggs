# code for implementing the metropolis algorithm (right now only for gauge fields)

# See Appendix A in PhD thesis by Mark B. Wurtz
# The variables follow exactly the notations there
function heatbath_gauge!(U::GaugeFieldSU2{d}, β::AbstractFloat, MonteCarloHits::Integer) where d
    X=UnitaryMatrix(zeros(2,2)) ## placeholder for Higgs fields
    for x in eachspacetimeindex(U)
        for μ in 1:d
            for i in 1:MonteCarloHits # Update each link MonteCarloHits time
                V=staple(U,μ,x)
                W=β*V+X
                α=sqrt(real(det(W)))
                a=a_generator(β,α) # Generate vectors on a sphere to genearte a new SU(2) matrix
                U(x)[μ]=W/α*SU2(a) # Update U
            end
        end
    end
end
export heatbath_gauge!

"""
    a_generator(β::Float64, α::Float64)

    Generate a 4-vector according to the Kennedy-Pendleton method
"""

function a_generator(β::Float64, α::Float64)
    a0=2.0

    while abs(a0)>1
        y=[1.0 1.0]

        while norm(y)>1
            y=rand(2)*2 .- 1
        end
        normy_sqr=norm(y)^2
#    println(normy_sqr)
        x1=rand()
        a0=1.0+(log(x1)+y[1]^2/normy_sqr*log(normy_sqr))/α # Eq. (A.33)
#    println(a0)
    end
    r=sqrt(1.0-a0^2)
    x=[1.0 1.0]
    while norm(x)>1
        x=rand(2)* .-1
    end
    normx_sqr=norm(x)^2
    a1=2*x[1]*sqrt(1-normx_sqr)*r
    a2=2*x[2]*sqrt(1-normx_sqr)*r
    a3=(1-2*normx_sqr)*r

    return [a0, a1, a2, a3 ]
end


function a_generator_kp(β::Float64, α::Float64)
    generated=false
    a0=0.
    while ! generated
        x=rand(4)
        a0=1+(log(x[1])+cos(2*pi*x[3])^2*log(x[2]))/α
        if a0 >= (2*x[4]^2-1)
            generated=true
        end
    end
    r=sqrt(1.0-a0^2)
    x=[1.0 1.0]
    while norm(x)>1
        x=rand(2)* .-1
    end
    normx_sqr=norm(x)^2
    a1=2*x[1]*sqrt(1-normx_sqr)*r
    a2=2*x[2]*sqrt(1-normx_sqr)*r
    a3=(1-2*normx_sqr)*r
    return [a0, a1, a2, a3]
end

function overrelaxation_gauge!(U::GaugeFieldSU2{d},β::AbstractFloat) where d
    X=UnitaryMatrix(zeros(2,2)) ## placeholder for Higgs fields
    for x in eachspacetimeindex(U)
        for μ in 1:d
            V=staple(U,μ,x)
            W=β*V+X
            α=det(W)
            U(x)[μ]=W'*U(x)[μ]'*W'/det(W)
        end
    end
end
export overrelaxation_gauge!

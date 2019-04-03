# code for implementing the metropolis algorithm (right now only for gauge fields)

# See Appendix A in PhD thesis by Mark B. Wurtz
function heatbath_U(U::GaugeFieldSU2{d}, μ::Integer, x::SVector, β::AbstractFloat) where d
        X=UnitaryMatrix(zeros(2,2)) ## placeholder for Higgs fields
        V=staple(U,μ,x)
        W=β*V+X
        α=sqrt(det(W))

        a=a_generator(β,α)


end


function a_generator(β::Float64, α::Float64)
    y=[1.0 1.0]

    while norm(y)>1
        y=rand(2)* .- 1
    end
    normy_sqr=norm(y)^2

    a0=1.0+(log(y[1])+y[1]^2/normy_sqr*log(normy_sqr))/α # Eq. (A.33)
    r=sqrt(1.0-a0^2)
    x=[1.0 1.0]
    while norm(x)>1
        x=rand(2)* .-1
    end
    normx_sqr=norm(x)^2
    a1=2*x[1]*sqrt(1-normx_sqr)/r
    a2=2*x[2]*sqrt(1-normx_sqr)/r
    a3=(1-2*normx_sqr)/r
    return [a0, a1, a2, a3 ]
end

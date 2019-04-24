
include("src/LGTHiggs.jl")
using .LGTHiggs

latt=LatticeToroidal(3,1.0, 6, 4)
A=rand(GaugeFieldSU2,latt)

heatbath_gauge!(A,10.0,5)
overrelaxation_gauge!(A,10.0)

function test_wislonloop()
    result=0.0
    count=0
    for μ in 1:3, ν in 1:3
    if μ != ν
        count+=1
        result+=plaquette(A,μ,ν,CartesianIndex(1,1,1))
    end
    end
    println(result/count)
    println(wilsonloop(A,CartesianIndex(1,1,1),1,1))
end
test_wislonloop()

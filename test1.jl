
include("src/LGTHiggs.jl")
using .LGTHiggs

latt=LatticeToroidal(3,1.0, 6, 4)
A=rand(GaugeFieldSU2,latt)

heatbath_gauge!(A,10.0,5)
overrelaxation_gauge!(A,10.0)

# Code for measurements

#





function wilsonloop(U::GaugeFieldSU2{d}, origin::CartesianIndex{d}, width::Int, length::Int) where d

    result=0.0
    count=0
    for μ in 1:d, ν in 1:d
        if μ != ν
            count+=1
            top_product=UnitaryMatrix(I, 2,2)
            bottom_product=UnitaryMatrix(I, 2,2)
            left_product=UnitaryMatrix(I, 2,2)
            right_product=UnitaryMatrix(I, 2,2)
            corner=shift(shift(origin,μ,width),ν,length)

            for x in 1:width
            #    println("$x, $μ, $ν")
                top_product*=U[μ](shift(origin,μ,x-1))
                bottom_product *= U[μ](shift(corner,μ,-x))'
            end
            for y in 1:length
            #    println("$y, $μ, $ν")
                right_product*=U[ν](shift(corner,ν,y-length-1))
                left_product*=U[ν](shift(origin,ν,length-y))'
            end
            result+=tr(top_product*right_product*bottom_product*left_product)
        end
    end # μ,ν iteration
    result=real(result/2.0/count)
end

export wilsonloop


function polyakovloop(U::GaugeFieldSU2{d}) where d

# slicing the spatial dimensions

    spaceindices=CartesianIndices(size(U.lattice)[1:end-1])
    poly=0.0
    for x in spaceindices
        product=UnitaryMatrix(I, 2,2)
        origin=CartesianIndex((x.I...,1))
        for t in 1:size(U.lattice)[end]
            product*=U[d](shift(origin,d,1))
        end
        poly+=0.5*tr(product)
    end
    return poly/=length(spaceindices) # spatial volume

end
export polyakovloop

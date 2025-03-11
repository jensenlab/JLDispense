

struct CobraHead <: TransferHead
    nozzles::AbstractArray{Nozzle}
end 

const cobra_nozzle = Nozzle(0.3u"µL",40u"µL",800u"µL",20u"µL",1.1,false,true)


const cobra_head = fill(cobra_nozzle,4,1)



mask(head::CobraHead,labware::JLIMS.WP96) = (4,96)



function mask(head::CobraHead,labware::JLIMS.WP96) 
    


    

    











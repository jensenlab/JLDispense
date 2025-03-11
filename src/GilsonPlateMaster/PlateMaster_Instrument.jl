

const platemaster_nozzle = Nozzle(2u"µL",200u"µL",200u"µL",0u"µL",1,false,false)



struct PlateMasterHead<: TransferHead
    nozzles::AbstractArray{Nozzle}
end 

const platemaster_head =fill(platemaster_nozzle,8,12)



mask_shape(head::PlateMasterHead,labware::JLIMS.Labware) = nothing 


mask_shape(head::PlateMasterHead,labware::JLIMS.WP96) = (1,1) 

mask_shape(head::PlateMasterHead,labware::JLIMS.WP384) = (2,2)



function mask(head::PlateMasterHead,labware::JLIMS.WP96,r::Integer,c::Integer) 

    shape= mask_shape(head,labware) 

    return [(1,1)]
end 


function mask(head::PlateMasterHead,labware::JLIMS.WP384,r::Integer,c::Integer)
    shape=mask_shape(head,labware)

    row= shape[1] -(r % shape[1]) 
    col = shape[2] -(c % shape[2]) 
    return [(row,col)]
end 




    


    

    



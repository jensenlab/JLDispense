struct Instrument 

    deck::Deck
    head::Head 
end 





struct DeckPosition 
    name::AbstractString
    slots::Tuple{Integer,Integer}
    compatible_containers::Set{Type{<:JLIMS.Labware}}
end 


Deck = AbstractArray{DeckPosition}

abstract type Head end 

abstract type TransferHead <:Head

mask_shape(::TransferHead,::JLIMS.Labware)=() 






struct Nozzle
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    is_discrete::Bool
    multidispense::Bool
end 



 


function mask(::TransferHead,::JLIMS.Labware,::Integer,::Integer) = [()]






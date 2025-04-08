





abstract type Head end 

abstract type TransferHead <: Head end 

abstract type FixedTransferHead <: TransferHead end 

abstract type AdjustableTransferHead <: TransferHead end 

nozzles(x::TransferHead)=x.nozzles


abstract type ReadHead <: Head end 

abstract type MovementHead <: Head end 




abstract type Tool end 

abstract type TransferTool <: Tool end 

abstract type MovementTool <: Tool end 

abstract type ReadTool <: Tool end 



struct Nozzle <: TransferTool
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    is_discrete::Bool
    multidispense::Bool
end 




abstract type DeckPosition end 
 
struct TransferPosition <: DeckPosition 
    can_aspirate::Bool
    can_dispense::Bool
    slots::Tuple{Integer,Integer}
    compatible_labware::Set{Type{<:Labware}}
end 

Deck = AbstractArray{<:DeckPosition}

struct Instrument{H<:Head} 
    head::H
    deck::Deck
end 
    


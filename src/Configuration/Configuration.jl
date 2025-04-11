
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

abstract type Nozzle <: TransferTool 

struct ContinuousNozzle <: Nozzle 
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    multidispense::Bool
end 


struct DiscreteNozzle <: Nozzle
    minVol::Unitful.Volume
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
end 


 
abstract type DeckPosition end 


struct UnconstrainedDeckPosition <:DeckPosition 
end 

const single_channel_position = DeckPosition(Set([Labware]))



labware(x::DeckPosition) = x.labware


function can_place(l::Labawre,d::DeckPosition) 
    return any(map(x->typeof(l) <: x),labware(d))
  end

function can_place(l::Labware,d::UnconstrainedDeckPosition)
    return true 
end 



function can_aspirate(h::Head, d::DeckPosition,l::Labware) 
    return false
  end
  function can_dispense(h::CobraHead,d::CobraDeckPosition,l::Labware) 
    return false
  end
  function can_move(h::Head,d::DeckPosition,l::Labware)
    return false 
  end
  function can_read(h::Head,d::DeckPosition,l::Labware)
    return false
  end 


Deck{T} = Vector{T<:DeckPosition}

abstract type InstrumentSettings end 

struct Configuration{H<:Head,D<:Deck,S<:InstrumentSettings}  
    head::H
    deck::D
    settings::S
end 
    
settings(x::Configuration)=x.settings
head(x::Configuration)=x.head
deck(x::Configuration)=x.deck

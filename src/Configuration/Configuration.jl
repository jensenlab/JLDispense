
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

abstract type Nozzle <: TransferTool end 

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
  name::String
end 

struct SBSPosition <: DeckPosition 
  name::String 
  labware::Set{<:Type{Labware}}
  slots::Tuple{Int,Int}
end 






SlottingDict = Dict{Labware,Tuple{DeckPosition,Integer}}

length(x::Labware) = prod(shape(x))

labware(x::DeckPosition) = x.labware
slots(x::DeckPosition) = x.slots 



function can_place(l::Labware,d::DeckPosition) 
    for x in collect(labware(d))
      if typeof(l) <: x 
        return true 
      end 
    end 
    return false 
  end

function can_place(l::Labware,d::UnconstrainedDeckPosition)
    return true 
end 



function can_aspirate(h::Head, d::DeckPosition,l::Labware) 
    return false
  end
  function can_dispense(h::Head,d::DeckPosition,l::Labware) 
    return false
  end
  function can_move(h::Head,d::DeckPosition,l::Labware)
    return false 
  end
  function can_read(h::Head,d::DeckPosition,l::Labware)
    return false
  end 


Deck{T} = AbstractArray{T} where T<:DeckPosition

abstract type InstrumentSettings end 

struct Configuration{H<:Head,D<:Deck,S<:InstrumentSettings}  
    head::H
    deck::D
    settings::S
end 
    
settings(x::Configuration)=x.settings
head(x::Configuration)=x.head
deck(x::Configuration)=x.deck

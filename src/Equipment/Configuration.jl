
abstract type Head end 

abstract type TransferHead <: Head end 

abstract type FixedTransferHead <: TransferHead end 

abstract type AdjustableTransferHead <: TransferHead end 

channels(x::TransferHead)=x.channels


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


struct EmptyPosition <: DeckPosition 
  name::String
end 


struct UnconstrainedPosition <:DeckPosition 
  name::String
  aspirate::Bool
  dispense::Bool
  move::Bool
  read::Bool
  plotting_shape::String
end 




struct SBSPosition <: DeckPosition 
  name::String 
  labware::Set{Type{<:Labware}}
  slots::Tuple{Int,Int}
  aspirate::Bool
  dispense::Bool
  move::Bool
  read::Bool
  plotting_shape::String
end 


struct StackPosition <: DeckPosition 
  name::String 
  labware::Set{Type{<:Labware}}
  slots::Tuple{Int,Int}
  aspirate::Bool
  dispense::Bool
  move::Bool
  read::Bool
  plotting_shape::String
end 













name(x::DeckPosition) = x.name 
name(::EmptyPosition)= "" 
length(x::Labware) = prod(JLIMS.shape(x))

labware(x::DeckPosition) = x.labware
slots(x::DeckPosition) = x.slots 
labware(::UnconstrainedPosition) = Set([JLIMS.Labware])
labware(::EmptyPosition) = Set{Type{<:Labware}}()
slots(::EmptyPosition)= ()
slots(::UnconstrainedPosition) = (4,6)




function can_place(l::Labware,d::DeckPosition) 
    for x in collect(labware(d))
      if typeof(l) <: x 
        return true 
      end 
    end 
    return false 
  end

function can_place(l::Labware,d::UnconstrainedPosition)
    return true 
end 

function can_place(l::Labware,d::EmptyPosition)
  return false 
end 


function can_aspirate(h::TransferHead,d::DeckPosition,l::Labware)
  return can_place(l,d) && d.aspirate == true 
end 

function can_dispense(h::TransferHead,d::DeckPosition,l::Labware)
  return can_place(l,d)  && d.dispense==true 
end 

function can_move(h::MovementHead,d::DeckPosition,l::Labware)
  return can_place(l,d)  && d.move ==true 
end 

function can_read(h::ReadHead,d::DeckPosition,l::Labware)
  return can_place(l,d) && d.read==true 
end 

can_aspirate(h::Head,d::DeckPosition,l::Labware) = false 
can_dispense(h::Head,d::DeckPosition,l::Labware) = false 
can_move(h::Head,d::DeckPosition,l::Labware) = false 
can_read(h::Head,d::DeckPosition,l::Labware) = false 


SlottingDict = Dict{Labware,Tuple{DeckPosition,Integer}}

Deck{T} = AbstractArray{T} where T<:DeckPosition

abstract type InstrumentSettings end 

struct Configuration{H<:Head,D<:Deck,S<:InstrumentSettings}  
    name::String
    head::H
    deck::D
    settings::S
end 
    
settings(x::Configuration)=x.settings
head(x::Configuration)=x.head
deck(x::Configuration)=x.deck
name(x::Configuration)=x.name





function masks(h::Head,l::Labware)

  Ma(w::Integer,p::Integer,c::Integer)= false 
  return Ma,Ma,(0,0),(0,0)
end


function can_place(labware::Labware,deck::Deck)

  for pos in deck 
    if can_place(labware,pos)
      return true 
    end 
  end
  return false 
end 
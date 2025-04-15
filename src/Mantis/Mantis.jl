





#=
struct Nozzle
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    is_discrete::Bool
    multidispense::Bool
end=#
const mantis_lv_nozzle = DiscreteNozzle(0.1u"µL",10u"µL",1.1)
const mantis_hv_nozzle = DiscreteNozzle(1u"µL",25u"µL",1.1)

abstract type MantisHead <: TransferHead end 

struct MantisLVHead <: MantisHead 
    nozzles::Nozzle
    MantisLVHead() = new(mantis_lv_nozzle)
end 

struct MantisHVHead <: MantisHead 
    nozzles::Nozzle
    MantisHVHead() =new(mantis_hv_nozzle)
end 

struct MantisSettings <:InstrumentSettings
    MantisSettings() = new()
end 
abstract type MantisDeckPosition <: DeckPosition end 
struct MantisLC3Position <: MantisDeckPosition
    labware::Set{Type{<:Labware}}
    MantisLC3Position()=new(Set([JLConstants.TipReservior,JLConstants.Conical15]))
end 

struct MantisHVPosition <: MantisDeckPosition
    labware::Set{Type{<:Labware}}
    MantisHVPosition()=new(Set([JLConstants.Conical]))
end 

struct MantisMainPosition <: MantisDeckPosition
    labware::Set{Type{<:Labware}}
    MantisMainPosition()=new(Set([JLConstants.MicroPlate,JLConstants.BreakawayPCRPlate]))
end  




const mantis_deck = vcat(fill(MantisLC3Position(),24),fill(MantisHVPosition(),4),MantisMainPosition())


MantisConfiguration = Configuration{MantisHead,Deck,MantisSettings}

const mantis_lv = MantisConfiguration(MantisLVHead(),mantis_deck,MantisSettings())
const mantis_hv = MantisConfiguration(MantisHVHead(),mantis_deck,MantisSettings())

### define deck access functions 


function can_aspirate(h::MantisHead, d::MantisDeckPosition,l::Labware) 
    return false 
end 
function can_dispense(h::MantisHead,d::MantisDeckPosition,l::Labware)
    return false
end 
function can_dispense(h::MantisHead, d::MantisMainPosition,l::Labware) 
    return can_place(l,d)
end 

function can_aspirate(h::MantisLVHead, d::MantisLC3Position,l::Labware) 
    return can_place(l,d)
end 
function can_aspirate(h::MantisHVHead,d::MantisHVPosition,l::Labware)
    return can_place(l,d)
end
function can_move(h::MantisHead,d::DeckPosition,l::Labware)
  return false 
end
function can_read(h::CobraHead,d::DeckPosition,l::Labware)
  return false
end 



function masks(h::MantisHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj = Wi,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # Mantis cannot aspirate from well plates 
        return false
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    return Ma,Md
end 

function masks(h::MantisLVHead,l::Union{JLConstants.Conical15,JLConstants.TipReservior}) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # Mantis cannot aspirate from well plates 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        return false 
    end 
    return Ma,Md
end 

function masks(h::MantisHVHead,l::JLConstants.Conical50) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # Mantis cannot aspirate from well plates 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        return false 
    end 
    return Ma,Md
end 




"""
    dispense(config::MantisConfiguration, design::DataFrame,directory::AbstractString,protocol_name::AbstractString,destination::Labware)

Create Mantis dipsense instructions for SBS plate dispensing

## Arguments 
* `config`: A Configuration object defining the Mantis
* `design`: a (# of sources) x (# of destinations) dataframe containing the volume of each transfer in µL.
* `directory`: output file directory
* `protocol_name`: protocol name  
* `destinations`: A JLIMS.Labware destination object  
"""
function dispense(config::MantisConfiguration, design::DataFrame,directory::AbstractString,protocol_name::AbstractString,destination::Labware)
        R,C=destination.shape
        source_names=names(design)
        full_dir=joinpath(directory,protocol_name)
        if ~isdir(full_dir)
            mkdir(full_dir)
        end 
        filename=joinpath(full_dir,protocol_name*".dl.txt")
        n_stocks=ncol(dispenses)
        delay_header=vcat(n_stocks,repeat([0,""],n_stocks))
        outfile=open(filename,"w")
        platefilename=mantis_names[typeof(destination)]
        print(outfile,join(["[ Version: 5 ]"],'\t'),"\r\n")
        print(outfile,platefilename,"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")
        print(outfile,join([1],'\t'),"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")

        for i in 1:n_stocks 
            vols=ustrip.(uconvert((u"µL",),design[:,i]))
            vols=round.(reshape(vols,R,C),digits=1)
            print(outfile,join([names(design)[i],"","Normal"],'\t'),"\r\n")
            print(outfile,join(["Well",1],'\t'),"\r\n")
            for r in 1:R
                print(outfile,join(vols[r,:],'\t'),"\r\n")
            end 
        end 
        close(outfile)
        write(joinpath(full_dir,"config.json"),JSON.json(robot))
        return nothing 
end 


"""
    dispense(config::MantisConfiguration, design::DataFrame,directory::AbstractString,protocol_name::AbstractString,destination::Type{<:Labware})

Create Mantis dipsense instructions for SBS plate dispensing

## Arguments 
* `config`: A Configuration object defining the Mantis
* `design`: a (# of sources) x (# of destinations) dataframe containing the volume of each transfer in µL.
* `directory`: output file directory
* `protocol_name`: protocol name  
* `destinations`: A JLIMS.Labware destination subtype  
"""
function dispense(config::MantisConfiguration,design::DataFrame, directory::AbstractString,protocol_name::AbstractString,destination::Type{<:Labware})
    dest=destination(1,"example_instance")
    return dispense(config,desing,directory,protocol_name,dest)
end 



#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]

protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,mantis_default)

=# 
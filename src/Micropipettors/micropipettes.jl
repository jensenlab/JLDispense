#=
struct ContinuousNozzle
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    multidispense::Bool
end=#
const p1000_nozzle = ContinuousNozzle(100u"µL",1000u"µL",1000u"µL",0u"µL",1,false)
const p200_nozzle = ContinuousNozzle(20u"µL",200u"µL",200u"µL",0u"µL",1,false)
const p20_nozzle = ContinuousNozzle(2u"µL",20u"µL",20u"µL",0u"µL",1,false)
const p2_nozzle = ContinuousNozzle(0.2u"µL",2u"µL",2u"µL",0u"µL",1,false)
const p100_nozzle = ContinuousNozzle(10u"µL",100u"µL",100u"µL",0u"µL",1,false)
const p10_nozzle = ContinuousNozzle(1u"µL",10u"µL",10u"µL",0u"µL",1,false)



struct SingleChannelHead <: FixedTransferHead
    nozzles::Nozzle
end 

const p1000_head =SingleChannelHead(p1000_nozzle)
const p200_head =SingleChannelHead(p200_nozzle)
const p20_head =SingleChannelHead(p20_nozzle)
const p2_head =SingleChannelHead(p2_nozzle)

function masks(h::SingleChannelHead,l::JLIMS.Labware)
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj =shape(l)
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 












struct EightChannelHead <: FixedTransferHead 
    nozzles::Vector{Nozzle}
end 

const p1000_multichannel_head =EightChannelHead(fill(p1000_nozzle,8))
const p100_multichannel_head = EightChannelHead(fill(p100_nozzle,8))
const p10_multichannel_head = EightChannelHead(fill(p10_nozzle,8))



function masks(h::EightChannelHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj = 1,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return  wj == pn 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 

function masks(h::EightChannelHead,l::JLConstants.WP384) 
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj = 2,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi % 2 == 2-pm && wj == pn 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 






struct HumanProperties <: RobotProperties
    isdiscrete::Bool
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    minMass::Unitful.Mass 
    maxMass::Unitful.Mass
    compatible_stocks::Vector{DataType}
end 

    
mutable struct HumanConfiguration <: RobotConfiguration 
    comment::AbstractString
end

struct Human <:Robot
    name::AbstractString 
    properties::HumanProperties 
    configuration::HumanConfiguration
end 








human_default=Human("Default Human",
HumanProperties(false,0.1u"µL",100u"L",0.1u"mg",100u"kg",[JLIMS.SolidStock,JLIMS.LiquidStock,JLIMS.EmptyStock]),
HumanConfiguration("n/a")
)

omnipotent_robot= Human("Omnipotent_Robot",
HumanProperties(false,0u"µL",Inf*u"µL",0u"g",Inf * u"g",[JLIMS.SolidStock,JLIMS.LiquidStock,JLIMS.EmptyStock]),
HumanConfiguration("n/a")
)



function human(directory::AbstractString, protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Human) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    nrow(design) ==length(sources) ? nothing : error("number of rows in the design must equal the number of source stocks")
    ncol(design)== length(destinations) ? nothing : error("number of columns in the design must equal the number of destination stocks")

    transfer_table=make_transfer_table(sources,destinations,design)

    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
        mkdir(full_dir)
    end 
    CSV.write(joinpath(full_dir,"$(protocol_name).csv"),transfer_table)
    write(joinpath(full_dir,"config.json"),JSON.json(robot))
    return transfer_table
end 


function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Human;kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    design=dispense_solver(sources,destinations,robot,minimize_overdrafts!,minimize_labware!,minimize_sources!,minimize_transfers!;pad=1.1,kwargs...)
    protocol_name=random_protocol_name()
    transfer_table=human(directory,protocol_name,design,sources,destinations,robot)
    return protocol_name, transfer_table
end 




    
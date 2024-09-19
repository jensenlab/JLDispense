
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


#= Test Code


using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]

protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,human_default)




=#

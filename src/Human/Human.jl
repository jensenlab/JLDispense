
struct HumanProperties <: RobotProperties
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








human_default=Human("",
HumanProperties(0.1u"µL",100u"L",0.1u"mg",100u"kg",[JLIMS.SolidStock,JLIMS.LiquidStock,JLIMS.EmptyStock]),
HumanConfiguration("n/a")
)




function human_instructor(directory::AbstractString, protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Human) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
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


function human(directory::AbstractString,protocol_name::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Human;kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    design=dispense_solver(sources,destinations,robot,minimize_overdrafts!,minimize_labware!,minimize_sources!,minimize_transfers!;pad=1.1,kwargs...)
    transfer_table=human_instructor(directory,protocol_name,design,sources,destinations,robot)
    return transfer_table
end 


#= Test Code


stocks=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]


tt,pn=human("/Users/BDavid/Desktop/",stocks,stocks,human_default;quiet=true)




=#

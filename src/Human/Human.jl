
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


function human_mixer(sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Human;quiet=false,timelimit=10)
 
    transfers=DataFrame()

    for dest in destinations 
        allowed_ingredients=ingredients(dest.composition)
        srcs=filter(y->all(map(x->in(x,allowed_ingredients),ingredients(y.composition))),sources)

        concentrations=stock_concentration_array(srcs)
        target_quantities=stock_quantity_array(dest)
        concentrations=concentrations[:,DataFrames.names(target_quantities)]
        
        concentrations=Float64.(ustrip.(Matrix(concentrations)))
  
        target_quantities=ustrip.(Vector(target_quantities[1,:]))
        
        S,I=size(concentrations)

        model=Model(Gurobi.Optimizer)
        if quiet 
            set_silent(model)
        end 
        set_attribute(model,"TimeLimit",timelimit)

        
        @variable(model,q[1:S] >=0)
        @variable(model,Is[1:S],Bin)

        @constraint(model, q'*concentrations .== target_quantities') # find the quantity of each stock needed to make the target 

        for s in 1:S
            @constraint(model, !Is[s]=> {q[s]==0})
        end 

        @objective(model, Min,sum(Is))
        optimize!(model)

        quants=JuMP.value.(q)


        
        quantities=Any[]
        for stock in sources
            un=preferred_stock_quantity(stock)
            if stock in srcs 
                idx=findfirst(x->x==stock,srcs)
                val=quants[idx]*un
                if isa(val,Unitful.Mass)
                    if 0u"g" <val < robot.properties.minMass
                        error("The $val mass required for the stock in well #$(stock.well.id) is too small for a $(typeof(robot)) to accurately transfer")
                    elseif  val > robot.properties.maxMass
                        error("The $val mass required for the stock in well #$(stock.well.id) is too large for a $(typeof(robot)) to accurately transfer")
                    end 
                elseif isa(val,Unitful.Volume)
                    if 0u"L" <val < robot.properties.minVol
                        error("The $val volume required for the stock in well #$(stock.well.id) is too small for a $(typeof(robot)) to accurately transfer")
                    elseif  val > robot.properties.maxVol
                        error("The $val volume required for the stock in well #$(stock.well.id) is too large for a $(typeof(robot)) to accurately transfer")
                    end 
                end 
                push!(quantities,quants[idx]*un) # commit the transfer to the transfers vector for this destination stock 
                
            else
                push!(quantities,0*un) # commit a 0 to the transfers vector
            end 
        end 
        transfers[:,Symbol("Well$(dest.well.id)")]=quantities
    end 

    transfer_table=DataFrame(Source=Integer[],Destination=Integer[],Quantity=Real[],Unit=AbstractString[])
    r=nrow(transfers)
    c=ncol(transfers)
    for col in 1:c
        for row in 1:r 
            val=transfers[row,col]
            quantity=ustrip(val)
            if quantity==0 
                continue 
            else 
                source=sources[row].well.id 
                destination=destinations[col].well.id
                un=string(unit(val))
                push!(transfer_table,(source,destination,quantity,un))
            end 
        end 
    end 
    return transfer_table , (transfers,sources,destinations)

end 



function human_instructor(directory::AbstractString, protocol_name::AbstractString,design::DataFrame,sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Human)
    nrow(design) ==length(sources) ? nothing : error("number of rows in the design must equal the number of source stocks")
    ncol(design)== length(destinations) ? nothing : error("number of columns in the design must equal the number of destination stocks")

    transfer_table=DataFrame(Source=Integer[],Destination=Integer[],Quantity=Real[],Unit=AbstractString[])
    r=nrow(design)
    c=ncol(design)
    for col in 1:c
        for row in 1:r 
            val=design[row,col]
            quantity=ustrip(val)
            if quantity==0 
                continue 
            else 
                source=sources[row].well.id 
                destination=destinations[col].well.id
                un=string(unit(val))
                push!(transfer_table,(source,destination,quantity,un))
            end 
        end 
    end 
    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
        mkdir(full_dir)
    end 
    CSV.write(joinpath(full_dir,"$(protocol_name).csv"),transfer_table)
    write(joinpath(full_dir,"config.json"),JSON.json(robot))

end 


function human(directory::AbstractString,sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Human;kwargs...)
    tt,protocol=human_mixer(sources,destinations,robot;kwargs...)
    protocol_name=random_protocol_name()
    human_instructor(directory,protocol_name,protocol[1],protocol[2],protocol[3],robot)
    return tt ,protocol_name
end 


#= Test Code


stocks=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]


tt=human("/Users/BDavid/Desktop/",stocks,stocks,human_default;quiet=true)




=#

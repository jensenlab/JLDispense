
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


function human_mixer(sources::Vector{T},destinations::Vector{U},robot::Human;quiet=false,timelimit=10) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
 
    transfers=DataFrame()
    allowed_ingredients=unique(vcat(map(x->ingredients(x.composition),destinations)...))
    srcs=filter(y->all(map(x->in(x,allowed_ingredients),ingredients(y.composition))),sources)

    concentrations=stock_concentration_array(srcs)
    target_quantities=stock_quantity_array(destinations)
    concentrations=concentrations[:,DataFrames.names(target_quantities)]
    
    concentrations=Float64.(ustrip.(Matrix(concentrations)))

    target_quantities=Float64.(ustrip.(Matrix(target_quantities)))
    src_quants=map(x->x.quantity,srcs)
    pref_units=preferred_stock_quantity.(srcs)
    src_quants=uconvert.(pref_units,src_quants)
    sq=ustrip.(src_quants)
    S,I=size(concentrations)
    D,I=size(target_quantities)
    model=Model(Gurobi.Optimizer)
    if quiet 
        set_silent(model)
    end 
    set_attribute(model,"TimeLimit",timelimit)

        
    @variable(model,q[1:S,1:D] >=0)
    @variable(model,Iq[1:S,1:D],Bin)

    @constraint(model, q'*concentrations .== target_quantities) # find the quantity of each stock needed to make the target 
    for s in 1:S
        for d in 1:D
            @constraint(model, !Iq[s,d]=> {q[s,d]==0})
        end 
    end 
    for s in 1:S
        @constraint(model, sum(q[s,:]) <= sq[s])
    end 

    @objective(model, Min,sum(Iq))
    optimize!(model)

    quants=JuMP.value.(q)


    transfers=DataFrame()    
    
    for dest in destinations
        quantities=Any[]
        for src in sources
            un=preferred_stock_quantity(src)
        if src in srcs 
            s=findfirst(x->x==src,srcs)
            d=findfirst(x->x==dest,destinations)
            val=quants[s,d]*un
            if isa(val,Unitful.Mass)
                if 0u"g" <val < robot.properties.minMass
                    error("The $val mass required for the stock in well #$(dest.well.id) is too small for a $(typeof(robot)) to accurately transfer")
                elseif  val > robot.properties.maxMass
                    error("The $val mass required for the stock in well #$(dest.well.id) is too large for a $(typeof(robot)) to accurately transfer")
                end 
            elseif isa(val,Unitful.Volume)
                if 0u"L" <val < robot.properties.minVol
                        error("The $val volume required for the stock in well #$(dest.well.id) is too small for a $(typeof(robot)) to accurately transfer")
                elseif  val > robot.properties.maxVol
                        error("The $val volume required for the stock in well #$(dest.well.id) is too large for a $(typeof(robot)) to accurately transfer")
                end 
            end 
            push!(quantities,val) # commit the transfer to the transfers vector for this destination stock 
                
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



function human_instructor(directory::AbstractString, protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Human) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
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


function human(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Human;kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    tt,protocol=human_mixer(sources,destinations,robot;kwargs...)
    protocol_name=random_protocol_name()
    human_instructor(directory,protocol_name,protocol[1],protocol[2],protocol[3],robot)
    return tt ,protocol_name
end 


#= Test Code


stocks=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]


tt,pn=human("/Users/BDavid/Desktop/",stocks,stocks,human_default;quiet=true)




=#

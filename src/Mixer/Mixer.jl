

function make_transfer_table(sources::Vector{T},destinations::Vector{U},design::DataFrame) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
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
    return transfer_table
end 




function mixer(protocol_name::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Robot;directory=pwd(),kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}

    robot_type=typeof(robot) 
    tt=DataFrame() 
    if robot_type == Human 
        tt=human(directory,protocol_name,sources,destinations,robot;kwargs...) 
    else
        error("$robot_type is unsupported by the mixer")
    end 
     
    return tt
end 


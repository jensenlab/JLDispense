function is_compatible_source(stock::JLIMS.Stock,robot::Robot)
    a = typeof(stock) in robot.properties.compatible_stocks  # make sure the robot can work with the stocktype 
    compatible_positions=filter(x->x.can_aspirate==true,robot.properties.positions)
    if any(ismissing.(map(x->x.compatible_containers,compatible_positions)))
        return a 
    else 
        compatible_containers=unique(vcat(map(t->t.compatible_containers,compatible_positions)...))
        b= JLIMS.well(stock).container in compatible_containers # make sure the robot can use the container that the stock is in 
        return a && b
    end 
end 

function is_compatible_source(culture::JLIMS.Culture,robot::Robot)
    return is_compatible_source(culture.media,robot)
end 

function is_compatible_destination(stock::JLIMS.Stock,robot::Robot)
    a = typeof(stock) in robot.properties.compatible_stocks  # make sure the robot can work with the stocktype 
    compatible_positions=filter(x->x.can_dispense==true,robot.properties.positions)
    if any(ismissing.(map(x->x.compatible_containers,compatible_positions)))
        return a 
    else 
        compatible_containers=unique(vcat(map(t->t.compatible_containers,compatible_positions)...))
        b= JLIMS.well(stock).container in compatible_containers # make sure the robot can use the container that the stock is in 
        return a && b
    end 
end 

function is_compatible_destination(culture::JLIMS.Culture,robot::Robot)
    return is_compatible_destination(culture.media,robot)
end 


function cartesian(x::AbstractArray,i)
    return Tuple(CartesianIndices(x)[i])
end 

function linear(x::AbstractArray,idxs...)
    return LinearIndeices(x)[idxs...]
end 


function transfer_table(source::Labware,destination::Labware,design::DataFrame)
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
                source=JLIMS.location_id(source[cartesian(source,row)])
                destination=JLIMS.location_id(destination[cartesian(destination,col)])
                un=string(unit(val))
                push!(transfer_table,(source,destination,quantity,un))
            end 
        end 
    end 
    return transfer_table
end 
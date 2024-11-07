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
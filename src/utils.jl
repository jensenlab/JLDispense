function is_compatible(stock::JLIMS.Stock,robot::Robot;source=true)
    a = typeof(stock) in robot.properties.compatible_stocks  # make sure the robot can work with the stocktype 
    
    

    compatible_containers=unique(vcat(map(t->t.compatible_containers,filter(x->x.is_source==source && !ismissing(x.compatible_containers),robot.properties.positions))...))
    b= JLIMS.well(stock).container in compatible_containers # make sure the robot can use the container that the stock is in 

    return a && b
end 

function is_compatible(culture::JLIMS.Culture,robot::Robot;source=true)
    a = typeof(culture) in robot.properties.compatible_stocks

    compatible_containers=unique(vcat(map(t->t.compatible_containers,filter(x->x.is_source==source && !ismissing(x.compatible_containers),robot.properties.positions))...))
    b= JLIMS.well(culture).container in compatible_containers

    return a && b
end 
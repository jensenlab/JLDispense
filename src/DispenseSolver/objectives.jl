function minimize_transfers!(model) # only use in the context of dispense solver, this function relies on model variables defined in dispense solver
    qi=model[:qi] # the transfer indicator in the main model
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(qi))
    optimize!(model)
    ti=JuMP.value.(qi) 
    @constraint(model,sum(qi)<=sum(ti))
end 

function minimize_sources!(model) # only use in the context of dispense_solver, this function relies on model variables defined in dispense solver
    S,D,R=size(model[:shots])
    shots=model[:shots]
    @variable(model,source_indicator[1:S],Bin) # Define an indicator for whether a source s is active
    for s in 1:S
        @constraint(model, !source_indicator[s]=>{sum(shots[s,:,:]) == 0}) # turn the indicator on if there is a nonzero transfer from source s to any destination
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(source_indicator))
    optimize!(model)
    si=JuMP.value.(source_indicator)
    @constraint(model,sum(source_indicator)<=sum(si))
end 

function minimize_labware!(model)
    SL,DL=size(model[:lw])
    lw=model[:lw]
    @variable(model,labware_indicator[1:SL],Bin) # Define an indicator for whether a source s is active
    for l in 1:SL
        @constraint(model, !labware_indicator[l]=>{sum(lw[l,:]) == 0}) # turn the indicator on if there is a nonzero transfer from source s to any destination
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(labware_indicator))
    optimize!(model)
    li=JuMP.value.(labware_indicator)
    @constraint(model,sum(labware_indicator)==sum(li))
end 


function minimize_overdrafts!(model) 
    S=length(model[:caps])
    q=model[:q]
    S,D,M=size(q)
    caps=model[:caps]
    @variable(model, overdraft_indicator[1:S],Bin)
    for s in 1:S 
        @constraint(model, !overdraft_indicator[s]=>{sum(q[s,:,1:(M-1)])<=caps[s]})
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(overdraft_indicator))
    optimize!(model)
    oi=JuMP.value.(overdraft_indicator)
    @constraint(model,sum(overdraft_indicator)==sum(oi))
end 


function minimize_robots!(model) 
    shots=model[:shots]
    S,D,M=size(shots)
    R=M-1
    @variable(model, robot_indicator[1:S,1:R],Bin)
    for r in 1:R 
        for s in 1:S
            @constraint(model,!robot_indicator[s,r]=> {sum(shots[s,:,r]) == 0})
        end 
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(robot_indicator))
    optimize!(model)
    ri=JuMP.value.(robot_indicator)
    @constraint(model,sum(robot_indicator)==sum(ri))
end 
    

function enforce_maxShots!(model)
    S,D=size(model[:q])
    q=model[:q]
    for s in 1:S
        @constraint(model,q[s,:,:] .<= model[:maxShots][s,:])
    end 
    optimize!(model) 
end

function minimize_overshots!(model)
    shots=model[:shots]
    S,D,M=size(shots)
    maxShots=model[:maxShots]
    R=M-1
    @variable(model,overshot_indicator[1:S,1:D,1:R],Bin)
    for s in 1:S
        for d in 1:D
            for r in 1:R
                @constraint(model,!overshot_indicator[s,d,r] => {q[s,d,r]<=maxShots[s,r]})
            end
        end 
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(overshot_indicator))
    optimize!(model) 
    oi=JuMP.value.(overshot_indicator)
    @constraint(model,sum(overshot_indicator)==sum(oi))
end

function minimize_labware_crossover!(model)
    SL,DL=size(model[:lw])
    lw=model[:lw]
    @variable(model,crossover_indicator[1:SL,1:DL],Bin) # Define an indicator for whether a source s is active
    for sl in 1:SL
        for dl in 1:DL
            @constraint(model, !crossover_indicator[sl,dl]=>{lw[sl,dl] == 0}) # turn the indicator on if there is a nonzero transfer from source s to any destination
        end
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(crossover_indicator))
    optimize!(model)
    li=JuMP.value.(crossover_indicator)
    @constraint(model,sum(crossover_indicator)==sum(li))
end 


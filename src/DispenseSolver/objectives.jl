#=
# General format for new objective definitions. New objectives transform the model and should use the ! to indicate as such. 

function new_objective!(model)
#1) define any new cost-tracking variables
#2) write constraints to track to the new variables in terms of Q or V
#3) set the model objective
#4) optimize the model
#5) constrain future solutions to ensure they are optimal for this
objective
return model
end

=#



function min_operations!(model;numerical_tolerance::Real=1e-8,kwargs...) 
    V = model[:V] 
    M,N=size(V)
    VI = @variable(model, VI[1:M,1:N],Bin)
    for m in 1:M 
            for n in 1:N 
                    @constraint(model, VI[m,n] --> {V[m,n] >= numerical_tolerance})
            end 
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min, sum(VI))
    optimize!(model)
    VIval = sum(JuMP.value.(VI))
    @constraint(model,sum(VI) <= VIval)
    optimize!(model) # solve once again to avoid downstream optimize not called errors
    return model 
end 



function min_sources!(model;numerical_tolerance::Real = 1e-8,kwargs...)
    Q = model[:Q]
    W,W =size(Q)
    @variable(model,QsI[1:W],Bin)
    for i in 1:W 
            @constraint(model,QsI[i] --> {sum(Q[i,:]) >= numerical_tolerance})
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min, sum(QsI))
    optimize!(model)
    QsIval= sum(JuMP.value.(QsI))
    @constraint(model,sum(QsI) <= QsIval)
    optimize!(model)
    return model 
end 


function min_labware_crossover!(model;numerical_tolerance::Real=1e-8,kwargs...)
    Qlw = model[:Qlw]
    L,L = size(Qlw)
    @variable(model,QlwI[1:L,1:L],Bin)
    for i in 1:L 
        for j in 1:L 
            @constraint(model, QlwI[i,j] --> {Qlw[i,j] >= numerical_tolerance})
        end 
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min, sum(QlwI))
    optimize!(model)
    QlwIval= sum(JuMP.value.(QlwI))
    @constraint(model,sum(QlwI) <= QlwIval)
    optimize!(model)
    return model 
end

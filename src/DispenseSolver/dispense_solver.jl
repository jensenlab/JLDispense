

struct MixingError <:Exception
    msg::AbstractString
    constraints::Vector{JuMP.ConstraintRef}
 end
    
struct OverdraftError <: Exception 
    msg::AbstractString
    balances::Dict{JLIMS.Stock,Unitful.Quantity}
end 

struct MissingIngredientError <: Exception 
    msg::AbstractString
    ings::Vector{AbstractString}
end 

struct InsufficientIngredientError <: Exception 
    msg::AbstractString
    ings::Vector{AbstractString}
end 

struct ContainerError <: Exception 
    msg::AbstractString
end 

struct StockCompatibilityError <: Exception 
    msg::AbstractString 
end 


Base.showerror(io::IO , e::MixingError) = begin 
    println(io, e.msg) 
    for con in e.constraints 
        println(io,con)
    end 
end 
Base.showerror(io,::IO, e::OverdraftError)= print(io,e.msg)
Base.showerror(io,::IO, e::MissingIngredientError)= print(io,e.msg)
Base.showerror(io,::IO, e::InsufficientIngredientError)= print(io,e.msg)
Base.showerror(io,::IO, e::ContainerError)= print(io,e.msg)
Base.showerror(io,::IO, e::StockCompatibilityError)= print(io,e.msg)



function preferred_ingredient_quantity(ing::JLIMS.Ingredient)
    pref_ing_quant=Dict(:solid=>u"mg",:liquid=>u"µL",:organism=>u"OD*µL")
    return pref_ing_quant[ing.class]
end 

function preferred_stock_quantity(stock::JLIMS.Stock)
    pref_stock_quant=Dict(JLIMS.SolidStock=> u"mg",JLIMS.LiquidStock => u"µL")
    return pref_stock_quant[typeof(stock)]
end 



""" 
    concentration(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the concentration of an ingredient in a stock using the preferred units for that ingredient and stock 
    

"""
function concentration(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient) 
    if ingredient in ingredients(stock.composition) 
        return uconvert(preferred_ingredient_quantity(ingredient)/preferred_stock_quantity(stock),stock.composition.ingredients[ingredient])
    else
        return 0*preferred_ingredient_quantity(ingredient)/preferred_stock_quantity(stock)
    end 
end 

""" 
    quantity(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the quantity of an ingredient in a stock using the preferred units for that ingredient
    

"""
function quantity(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)
    if ingredient in ingredients(stock.composition)
        return uconvert(preferred_ingredient_quantity(ingredient),stock.composition.ingredients[ingredient]*stock.quantity)
    else
        return 0*preferred_ingredient_quantity(ingredient)
    end 
end 

function stock_ingredient_array(stocks::Vector{T},ingredients::Vector{JLIMS.Ingredient};measure=concentration) where T<: JLIMS.Stock # can return concentration or quantity 
    S=length(stocks)
    out=DataFrame()
        for i in ingredients 
            vals=Any[]
            for s in stocks
                push!(vals,measure(s,i))
            end 
            out[:,Symbol(i.name)]=vals
        end 

    return out
end 

            
function stock_ingredient_array(stock::JLIMS.Stock,ingredients::Vector{JLIMS.Ingredient})
    return stock_ingredient_array([stock],ingredients)
end 

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


function minimize_transfers!(model) # only use in the context of dispense solver, this function relies on model variables defined in dispense solver
    qi=model[:qi] # the transfer indicator in the main model
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(qi))
    optimize!(model)
    ti=JuMP.value.(qi) 
    @constraint(model,sum(qi)==sum(ti))
end 

function minimize_sources!(model) # only use in the context of dispense_solver, this function relies on model variables defined in dispense solver
    S,D=size(model[:q])
    q=model[:q]
    @variable(model,source_indicator[1:S],Bin) # Define an indicator for whether a source s is active
    for s in 1:S
        @constraint(model, !source_indicator[s]=>{sum(q[s,:]) == 0}) # turn the indicator on if there is a nonzero transfer from source s to any destination
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(source_indicator))
    optimize!(model)
    si=JuMP.value.(source_indicator)
    @constraint(model,sum(source_indicator)==sum(si))
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
    caps=model[:caps]
    @variable(model, overdraft_indicator[1:S],Bin)
    for s in 1:S 
        @constraint(model, !overdraft_indicator[s]=>{sum(q[s,:])<=caps[s]})
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(overdraft_indicator))
    optimize!(model)
    oi=JuMP.value.(overdraft_indicator)
    @constraint(model,sum(overdraft_indicator)==sum(oi))
end 




function enforce_maxShots!(model)
    S,D=size(model[:q])
    q=model[:q]
    for s in 1:S
        @constraint(model,q[s,:] .<= model[:maxShots][s])
    end 
    optimize!(model) 
end

function minimize_overshots!(model)
    S,D=size(model[:q])
    q=model[:q]
    maxShots=model[:maxShots]
    @variable(model,overshot_indicator[1:S,1:D],Bin)
    for s in 1:S
        for d in 1:D
            @constraint(model,!overshot_indicator[s,d] => {q[s,d]<=maxShots[s]})
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







function dispense_solver(sources::Vector{T},destinations::Vector{U},robot::Robot,secondary_objectives...;quiet::Bool=true,timelimit::Real=10,pad::Real=1.25,slack_tolerance::Real=0,overdraft_tolerance::Real=1e-8,require_nonzero::Bool=true,return_model::Bool=false,obj_tolerance=1e-6,obj_cutoff=1e-3,priority::Dict{JLIMS.Ingredient,UInt64}=Dict{JLIMS.Ingredient,UInt64}(),kwargs...) where {T<: JLIMS.Stock,U<:JLIMS.Stock}
   
    
    
    # check inputs for issues 
    pad >= 1 ? nothing : error("padding factor must be greater than or equal to 1")
    0 <= overdraft_tolerance ? nothing : error("overdraft tolerance must be nonnegative")
    0 <= slack_tolerance ? nothing : error("Slack tolerance must be nonnegative")
    # check sources and destination for compatibility issues. If sources have an issue, ignore them. If the destinations have an issue, throw an error. 
    available_sources= filter(x->in(typeof(x),robot.properties.compatible_stocks),sources)
    if hasproperty(robot.properties,:positions) 
        
        stock_compatibility=map(x->in(typeof(x),robot.properties.compatible_stocks),destinations)
       
        for d in eachindex(destinations) 
            if !stock_compatibility[d]
                throw(StockCompatibilityError("the robot $(robot.name) is not compatible with a $(typeof(destinations[d]))."))
            end 
        end 
        if !any(map(x->typeof(x.compatible_containers)==Missing,filter(x->x.is_source==false,robot.properties.positions)))
            compatible_destination_containers=unique(vcat(map(x->x.compatible_containers,filter(x->x.is_source==false && typeof(x.compatible_containers) != Missing,robot.properties.positions))...))
            container_compatibility=map(x->in(x.well.container,compatible_destination_containers),destinations)
            for d in eachindex(destinations)
                if !container_compatibility[d]
                    throw(ContainerError("the robot $(robot.name) is not compatible with the container $(destinations[d].well.container.name) as a destination"))
                end 
            end 
        end 

        if !any(map(x->typeof(x.compatible_containers)==Missing,filter(x->x.is_source==true,robot.properties.positions)))
            compatible_source_containers=unique(vcat(map(x->x.compatible_containers,filter(x->x.is_source==true,robot.properties.positions))...))
            available_sources= filter(x->in(x.well.container,compatible_source_containers),available_sources) 
        end 
    end
    if length(available_sources)==0
        throw(error("No compatible sources exist for $(robot.name)"))
    end 

    minVol= ustrip(uconvert(u"µL",hasproperty(robot.properties,:minVol) ? robot.properties.minVol : 0u"µL"))
    maxVol=ustrip(uconvert(u"µL",hasproperty(robot.properties,:maxVol) ? robot.properties.maxVol : Inf*u"µL"))
    minMass=ustrip(uconvert(u"mg",hasproperty(robot.properties,:minMass) ? robot.properties.minMass : 0u"g"))
    maxMass=ustrip(uconvert(u"mg",hasproperty(robot.properties,:maxMass) ? robot.properties.maxMass : Inf*u"g"))
    maxShotDict=Dict(JLIMS.SolidStock => maxMass, JLIMS.LiquidStock => maxVol)
    minShotDict=Dict(JLIMS.SolidStock => minMass, JLIMS.LiquidStock => minVol)
    minshotvals=map(x->minShotDict[typeof(x)],available_sources)
    maxshotvals=map(x->maxShotDict[typeof(x)],available_sources)
    
    # Gather all ingredients contained in the sources, destinations, and priority list 
    source_ingredients= unique(vcat(map(x->ingredients(x.composition),available_sources)...))
    destination_ingredients= unique(vcat(map(x->ingredients(x.composition),destinations)...))
    all_ingredients = unique(vcat(source_ingredients,destination_ingredients,collect(keys(priority))))
    # group all available_sources by labware 
    src_labware=unique(map(x->x.well.labwareid,available_sources))
    SL=length(src_labware)
    slw_idxs=[findall(x->x.well.labwareid==src_labware[l],available_sources) for l in 1:SL ]


    # group all destinations by labware 
    dest_labware=unique(map(x->x.well.labwareid,destinations))
    DL=length(dest_labware)
    dlw_idxs=[findall(x->x.well.labwareid==dest_labware[l],destinations) for l in 1:DL ]

    # get and convert the source stock quantities for overdraft constants 
    source_quantities = ustrip.(map(x->uconvert(preferred_stock_quantity(x),x.quantity),available_sources))

    # Create the concentration array for the sources and the target quantity array for the destinations
    source_concentrations=stock_ingredient_array(available_sources,all_ingredients) 
    sc = Float64.(Matrix(Unitful.ustrip.(source_concentrations)))





    #create the target array for the destinations 
    destination_quantities=stock_ingredient_array(destinations,all_ingredients;measure=quantity)
    dq = Float64.(Matrix(Unitful.ustrip.(destination_quantities)))
    # grab the source quantities in the same units as the ingredient arrays 

    S=length(available_sources)
    D=length(destinations)
    I=length(all_ingredients) 
    # check for missing source ingredients needed to complete the destination
    missing_ingredients=filter(x-> sum(Unitful.ustrip.(source_concentrations[:,Symbol(x.name)]))==0,destination_ingredients)

    if length(missing_ingredients) > 0 
            ings=map(x->x.name,missing_ingredients)
            throw(MissingIngredientError("No valid source of $(join(ings,",")) available to complete the dispenses",ings))
    end 

    # Update priorities: priority increases with decreasing level

    # Level 0: All level 0 ingredeients must hit their target concentration exactly, unless the robot is a discrete dispenser, then we allow for slack up to the objective cutoff. 
    # Level 1+: All other ingredients are scheduled sequentially by priority. we allow slack for all nonzero priorities, where we try to minimize the slack for each ingredient and then constrain the slack for future priority levels. 
    # Level 2^(64)-1: Maximum allowable priority level, typically a solvent like water that we would use to back fill will be assigned maximum priority level. 

    for ing in destination_ingredients
        if !in(ing,keys(priority))
            val=UInt64(0) # any ingredient needed for any of the destination stocks is assigned 0 level priority unless the user has already explicitly defined a priority for that ingredient
            if robot.properties.isdiscrete
                val=UInt64(1) # we can't hit targets exactly with discrete instruments, so we give destination target ingredients a priority level 1 
            end 
            priority[ing]=val 
        end 
    end 
    for ing in source_ingredients
        if !in(ing,keys(priority)) # if the user hasn't specified a source ingredient priority or included it in the design (see above), assume it is priority 0. These ingredients will have a target concentration of 0 that must be hit exactly (they will be blocked, regardless of a discrete or continuous robot
            priority[ing]=UInt64(0) 
        end 
    end 

    capvals=source_quantities/pad
    maxshot_param=[Inf for s in 1:S]
    minshot_param=[0 for s in 1:S]
    if robot.properties.isdiscrete
        maxshot_param=maxshotvals ./ minshotvals # convert all masses and volumes to shots for the discrete problem 
        minshot_param=minshotvals ./ minshotvals 
        capvals = capvals ./ minshotvals 
        sc = minshotvals .* sc
    else 
        maxshot_param = maxshotvals 
        minshot_param= minshotvals
    end 



    # initialize the JuMP model 
    model=Model(Gurobi.Optimizer) 
    if quiet 
        set_silent(model)
    end 
    set_attribute(model,"TimeLimit",timelimit)

    # Define constants 
    @variable(model, minShots[1:S] in Parameter.(minshot_param))
    @variable(model, maxShots[1:S] in Parameter.(maxshot_param))
    @variable(model, caps[1:S] in Parameter.(capvals)) # save a parameter for the volume of each available source

    # Define Model variables 
    @variable(model, q[1:S,1:D]>=0) # q[s,d] = quantity of stock s transfered to stock d. q is the volume/mass for continous problems and shots for discrete problems 
    @variable(model, qi[1:S,1:D],Bin) #  Indicator for whether a transfer in q is active
    @variable(model,slacks[1:D,1:I]) # slack variables to measure the difference between the dispenses and the target quantities for each ingredient of each destination
    @variable(model, lw[1:SL,1:DL]>=0) # continuous variable to measure the total volume dispensed from and to each labware 
     

    @constraint(model, q'*sc .- slacks .== dq ) # create the destinations with the sources, allowing for some slack


    if robot.properties.isdiscrete 
        println("discrete problem")
        set_integer.(q) # q is an integer for discrete problems 
        @constraint(model, q .>= qi) # enforce activity constraint
    else 
        println("semi-continuous problem")
        for s in 1:S 
            @constraint(model, q[s,:] .>= minshotvals[s]*qi[s,:]) # if a dispense is active it must be larger than the minimum shot volume for the robot. This formulates q as a semicontinuous variable from (minshotval, Inf)
        end 
    end
    for i in 1:I
        if priority[all_ingredients[i]] == 0
            for d in 1:D
                @constraint(model, slacks[d,i]==0) # priority 0 ingredients must hit the target exactly for each destination. The slack must be zero. We can only enforce this for semicontinuous dispensers. 
            end 
        end 
    end 
    
    for s in 1:S
        for d in 1:D 
            @constraint(model, !qi[s,d] => {q[s,d]== 0}) # tie the transfer indicator to the transfers
        end 
    end 

    if require_nonzero 
        for i in 1:I
            sources_with_ingredient=sc[:,i] .> 0 
            for d in 1:D
                if dq[d,i] > 0 
                    @constraint(model,sum(qi[:,d] .* sources_with_ingredient) >=1) # check that at least one transfer happens if an ingredient is needed in a destination, even if the optimial solution is to not dispense anything. 
                end 
            end 
        end 
    end 

     

    for sl in 1:SL 
        for dl in 1:DL
            @constraint(model, lw[sl,dl]==sum(q[slw_idxs[sl],dlw_idxs[dl]])) # constrain the labware volume to track the dispenses from each source labware to each destination labware
        end
    end 
    priority_levels= sort(unique(collect(values(priority))))
    for level in priority_levels  # pass through all priority levels from lowest to higest level


        weights= falses(I)
        for i in 1:I
            if priority[all_ingredients[i]] == level 
                weights[i]=true # activate the weight term for this ingredient on this pass
            end 
        end 
        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        if level==UInt(0)
            set_attribute(model,"Cutoff",obj_cutoff* sum(weights'.*dq.^2)) # reject any solutions that are larger than the objective tolerance, which is a percentage of the sum squared target quantities.
            set_attribute(model, "BestObjStop",obj_tolerance*sum(weights'.*dq.^2)) 
        end 
        @objective(model, Min,sum(weights'.*(slacks.^2))) # penalize large slacks, search for q that minimize slacks.
        optimize!(model)

        # Check for optimality and feasibility 
        term = termination_status(model)
        if term == MOI.OPTIMAL || term == MOI.OBJECTIVE_LIMIT
            if primal_status(model)==MOI.FEASIBLE_POINT
                println("Optimal Solution Found For Level $level")
            elseif primal_status(model)==MOI.NO_SOLUTION
                throw(error("No solution exists that is better than the set objective cutoff $(obj_cutoff*100)%."))
            end 
        elseif term == MOI.TIME_LIMIT
            if primal_status(model) == MOI.FEASIBLE_POINT
                @warn "a solution was found for level $level, but it may be sub-optimal because the solver stopped due to reaching its $(timelimit)s  time limit."
            else 
                throw(error("the solver was unable to find a feasible solution for level $level in the $(timelimit)s time limit."))
            end 
        elseif term == MOI.INFEASIBLE || term == MOI.INFEASIBLE_OR_UNBOUNDED
            println("infeasible solution")
            compute_conflict!(model)
            out_cons=ConstraintRef[]
            if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
                cons=all_constraints(model; include_variable_in_set_constraints=false)
                for con in cons 
                    try get_attribute(con,MOI.ConstraintConflictStatus())
                        if get_attribute(con,MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                            push!(out_cons,con)
                        end 
                    catch 
                    end 
                end 
                throw(MixingError("the requested destination stocks cannot be made using this combination of sources and robot \n the following constraints are in conflict:",out_cons))
            end 
        end 
        current_slacks=abs.(JuMP.value.(slacks))
         
        for i in 1:I
            if weights[i]
                delta=slack_tolerance * dq[:,i] # delta is the tolerance we give to updating the slack in higher priority levels, it is some fraction of the dispense target quantity for every destination. Wiggle room for the slack in future iterations 
                current_slack=current_slacks[:,i]
                @constraint(model, slacks[:,i] .>= -current_slack .- delta)
                @constraint(model, slacks[:,i] .<= current_slack .+ delta)
            end 
        end 
       
        
    end 
    #reset the optimizer to remove the objective cutoff since we are switching objectives 
    set_optimizer(model,Gurobi.Optimizer)
    set_attribute(model,"TimeLimit",timelimit)
    optimize!(model) # resolve one last time with the final slack constraints -> we need to optimize before querying results for the secondary objectives 

    


    #solve for secondary objectives
    for obj! in secondary_objectives # secondary objectives only improve the solution quality after we have found a valid solution to the dispense problem 
        obj!(model) # execute the secondary objectives sequentially, objectives update the model in place. The general outline for a secondary objective function is:  Define variable to measure objective -> Optimize -> Constrain variable 
        optimize!(model) # re-optimize before qeurying again
    end 


    cap_vals=JuMP.value.(caps) 
    quants=JuMP.value.(q)

    sources_needed=sum(quants,dims=2)

    overdrafts = sources_needed .- cap_vals
    if any(overdrafts .>overdraft_tolerance) 
        overdraft_dict=Dict{JLIMS.Stock,Unitful.Quantity}()
        for i in findall(x-> x > 0 ,overdrafts)
            od=overdrafts[i]
            if robot.properties.isdiscrete
                od = od * minshotvals[i]
            end 
            overdraft_dict[available_sources[i]] = od * preferred_stock_quantity(available_sources[i])
        end 
        throw(OverdraftError("Refills are needed for $(length(collect(keys(overdraft_dict)))) stocks:",overdraft_dict))
    end 


    transfers=DataFrame()    
    for dest in destinations
        quantities=Any[]
        d=findfirst(x->x==dest,destinations)
        for src in sources
            un=preferred_stock_quantity(src)
            if in(src,available_sources)
                s=findfirst(x->x==src,available_sources)
                val=quants[s,d]
                if robot.properties.isdiscrete
                    val=Int64(round(val))*minshotvals[s]
                end 
                push!(quantities,val*un) # commit the transfer to the transfers vector for this destination stock 
            else
                push!(quantities,0*un) # commit a 0 to the transfers vector
            end 
        end 
        transfers[:,Symbol("Well$(dest.well.id)")]=quantities
    end 
    println("Solution off target by $(sum(JuMP.value.(slacks).^2)/sum(dq.^2)*100)%")

    if return_model
        return transfers,model
    else 
        return transfers
    end 
end 
    







#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
t,m=dispense_solver(sources,destinations,cobra_default;return_model=true)
t,m=dispense_solver(sources,destinations,mantis_default;return_model=true)
=#
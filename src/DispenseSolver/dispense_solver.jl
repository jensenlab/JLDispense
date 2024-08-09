

struct MixingError <:Exception
    msg::AbstractString
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


Base.showerror(io::IO , e::MixingError) = print(io, e.msg)
Base.showerror(io,::IO, e::OverdraftError)= print(io,e.msg)
Base.showerror(io,::IO, e::MissingIngredientError)= print(io,e.msg)
Base.showerror(io,::IO, e::InsufficientIngredientError)= print(io,e.msg)
Base.showerror(io,::IO, e::ContainerError)= print(io,e.msg)
Base.showerror(io,::IO, e::StockCompatibilityError)= print(io,e.msg)



function preferred_ingredient_quantity(ing::JLIMS.Ingredient)
    pref_ing_quant=Dict(:solid=>u"g",:liquid=>u"µL",:organism=>u"OD*µL")
    return pref_ing_quant[ing.class]
end 

function preferred_stock_quantity(stock::JLIMS.Stock)
    pref_stock_quant=Dict(JLIMS.SolidStock=> u"g",JLIMS.LiquidStock => u"µL")
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


function minimize_transfers!(model) # only use in the context of dispense solver, this function relies on model variables defined in dispense solver
    S,D=size(model[:q])
    q=model[:q]
    @variable(model,transfer_indicator[1:S,1:D],Bin) # Define an indicator for whether a transfer (q) is active
    for s in 1:S
        for d in 1:D
            @constraint(model, !transfer_indicator[s,d]=>{q[s,d] == 0}) # turn the indicator on if there is a nonzero transfer
        end 
    end 
    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    @objective(model, Min,sum(transfer_indicator))
    optimize!(model)
    ti=JuMP.value.(transfer_indicator) 
    @constraint(model,sum(transfer_indicator)==sum(ti))
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
    L=length(model[:lw])
    lw=model[:lw]
    @variable(model,labware_indicator[1:L],Bin) # Define an indicator for whether a source s is active
    for l in 1:L
        @constraint(model, !labware_indicator[l]=>{lw[l] == 0}) # turn the indicator on if there is a nonzero transfer from source s to any destination
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




function dispense_solver(sources::Vector{T},destinations::Vector{U},robot::Robot,secondary_objectives...;quiet::Bool=true,timelimit::Real=10,pad::Real=1.25,slack_tolerance::Real=0,overdraft_tolerance::Real=1e-8,priority::Dict{JLIMS.Ingredient,UInt64}=Dict{JLIMS.Ingredient,UInt64}(),kwargs...) where {T<: JLIMS.Stock,U<:JLIMS.Stock}
    # check inputs for issues 
    pad >= 1 ? nothing : error("padding factor must be greater than or equal to 1")
    0 <= tolerance ? nothing : error("Slack tolerance must be nonnegative")
    # check sources and destination for compatibility issues. If sources have an issue, ignore them. If the destinations have an issue, throw an error. 
    available_sources= filter(x->in(typeof(x),robot.properties.compatible_stocks),sources)
    if hasproperty(robot.properties,:positions)
        compatible_containers=unique(map(x->x.compatible_containers,robot.properties.positions))
        stock_compatibility=map(x->!in(typeof(x),robot.properties.compatible_stocks),destinations)
        container_compatibility=map(x->!in(x.well.container,compatible_containers),destinations)
        for d in eachindex(destinations) 
            if !stock_compatibility[d]
                throw(StockCompatibilityError("the robot $(robot.name) is not compatible with a $(typeof(destinations[d]))."))
            end 
            if !container_compatibility[d]
                throw(ContainerError("the robot $(robot.name) is not compatible with the container $(destinations[d].well.container.name)"))
            end 
        end 
    
        available_sources= filter(x->in(x.well.container,compatible_containers),avialable_sources) 
    end
    
    # Gather all ingredients contained in the sources, destinations, and priority list 
    source_ingredients= unique(vcat(map(x->ingredients(x.composition),available_sources)...))
    destination_ingredients= unique(vcat(map(x->ingredients(x.composition),destinations)...))
    all_ingredients = unique(vcat(source_ingredients,destination_ingredients,collect(keys(priority))))

    # group all available_sources by labware 
    src_labware=unique(map(x->x.well.labwareid,available_sources))
    L=length(src_labware)
    lw_idxs=[findall(x->x.well.labwareid==src_labware[l],available_sources) for l in 1:L ]

    # get and convert the source stock quantities for overdraft constants 
    source_quantities = ustrip.(map(x->uconvert(preferred_stock_quantity(x),x.quantity),available_sources))

    # Create the concentration array for the sources and the target quantity array for the destinations
    source_concentrations=stock_ingredient_array(available_sources,all_ingredients) 
    sc = Float64.(Matrix(Unitful.ustrip.(source_concentrations)))
    
    destination_quantities=stock_ingredient_array(destinations,all_ingredients;measure=quantity)
    dq = Float64.(Matrix(Unitful.ustrip.(destination_quantities)))
    # grab the source quantities in the same units as teh ingredient arrays 

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

    # Level 0: All level 0 ingredeients must hit their target concentration exactly 
    # Level 1+: All other ingredients are scheduled sequentially by priority. we allow slack for all nonzero priorities, where we try to minimize the slack for each ingredient and then constrain the slack for future priority levels. 
    # Level 2^(64)-1: Maximum allowable priority level, typically a solvent like water that we would use to back fill will be assigned maximum priority level. 

    for ing in destination_ingredients
        if !in(ing,keys(priority))
            priority[ing]=UInt64(0) # any ingredient needed for any of the destination stocks is assigned 0 level priority unless the user has already explicitly defined a priority for that ingredient
        end 
    end 
    for ing in source_ingredients
        if !in(ing,keys(priority)) # if the user hasnt specified a source ingredient priority or included it in the design (see above), assume it is priority 0. These ingredients will have a target concentration of 0 that must be hit exactly (they will be blocked)
            priority[ing]=UInt64(0) 
        end 
    end 


    # initialize the JuMP model 
    model=Model(Gurobi.Optimizer) 
    if quiet 
        set_silent(model)
    end 
    set_attribute(model,"TimeLimit",timelimit)

    @variable(model, q[1:S,1:D]>=0) # q[s,d] = quantity of stock s transfered to stock d 
    @variable(model,slacks[1:D,1:I]) # slack variables to measure the difference between the dispenses and the target quantities for each ingredient of each destination
    @variable(model, lw[1:L]>=0) # continuous variable to measure the total volume dispensed from each piece of labware 
    @variable(model, caps[i=1:S] in Parameter(source_quantities[i]/pad)) # save a parameter for the volume of each available source 
    
    @constraint(model, q'*sc .- slacks .== dq ) # create the destinations with the sources, allowing for some slack

    for l in 1:L 
        @constraint(model, lw[l]==sum(q[lw_idxs[l],:])) # constrain the labware volume to track the dispenses 
    end 

    for i in 1:I
        if priority[all_ingredients[i]] == 0
            for d in 1:D
                @constraint(model, slacks[d,i]==0) # priority 0 ingredients must hit the target exactly for each destination. The slack must be zero.
            end 
        end 
    end 
    
    priority_levels= sort(unique(collect(values(priority))))
    for level in priority_levels  # pass through all priority levels from lowest to higest level


        weights= falses(D,I)
        for i in 1:I
            if priority[all_ingredients[i]] == level 
                weights[:,i].=true # activate the weight term for this ingredient on this pass
            end 
        end 

        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        @objective(model, Min,sum(weights.*(slacks.^2))) # penalize large slacks, search for q that minimize slacks.
        optimize!(model)
        if !JuMP.is_solved_and_feasible(model)
            throw(MixingError("the requested stocks cannot be made using this combination of sources"))
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

    optimize!(model) # resolve one last time with the final slack constraints -> we need to optimize before querying results for the secondary objectives 
    if !JuMP.is_solved_and_feasible(model)
        throw(MixingError("the requested stocks cannot be made using this combination of sources"))
    end

    for objective! in secondary_objectives # secondary objectives only improve the solution quality after we have found a valid solution to the dispense problem 
        objective!(model) # execute the secondary objectives in order, objectives update the model in place. 
        optimize!(model) # re-optimize before qeurying again
    end 


    cap_vals=JuMP.value.(caps) 
    quants=JuMP.value.(q)

    sources_needed=sum(quants,dims=2)

    overdrafts = sources_needed .- cap_vals
    if any(overdrafts .>overdraft_tolerance) 
        overdraft_dict=Dict{JLIMS.Stock,Unitful.Quantity}()
        for i in findall(x-> x > 0 ,overdrafts)
            overdraft_dict[available_sources[i]] = overdrafts[i] * preferred_stock_quantity(available_sources[i])
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
                val=quants[s,d]*un
                push!(quantities,val) # commit the transfer to the transfers vector for this destination stock 
            else
                push!(quantities,0*un) # commit a 0 to the transfers vector
            end 
        end 
        transfers[:,Symbol("Well$(dest.well.id)")]=quantities
    end 
    return transfers
end 
    


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



preferred_quantity_unit(ing::JLIMS.Solid) =u"mg"
preferred_quantit_unit(ing::JLIMS.Liquid) = u"µL"


preferred_quantity_unit(stock::JLIMS.Stock)= u"µL"
preferred_quantity_unit(stock::JLIMS.Mixture) =u"mg"


const chemical_access_dict=Dict(
    JLIMS.Solid => JLIMS.solids 
    JLIMS.Liquid => JLIMS.liquids
)


""" 
    concentration(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the concentration of an ingredient in a stock using the preferred units for that ingredient and stock 
    

"""
function concentration(stock::JLIMS.Stock,ingredient::JLIMS.Chemical) 
    if ingredient in stock 
        return uconvert(preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient]/quantity(stock))
    else
        return 0*preferred_quantity(ingredient)/preferred_quantity(stock)
    end 
end 


""" 
    quantity(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the quantity of an ingredient in a stock using the preferred units for that ingredient
    

"""
function quantity(stock::JLIMS.Stock,ingredient::JLIMS.Chemical)
    if ingredient in ingredients(stock.composition)
        return uconvert(preferred_quantity(ingredient),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient])
    else
        return 0*preferred_quantity(ingredient)
    end 
end 


function ingredient_array(stocks::Vector{<:JLIMS.Stock},ingredients::Vector{<:JLIMS.Chemical};measure::Function=concentration) # can return concentration or quantity 
    S=length(stocks)
    out=DataFrame()
        for i in ingredients 
            vals=Any[]
            for s in stocks
                push!(vals,measure(s,i))
            end 
            out[:,Symbol(JLIMS.name(i))]=vals
        end 

    return out
end 

        


function organism_array(stocks::Vector{<:JLIMS.Stock},orgs::Vector{JLIMS.Organism}) where T<:JLIMS.Culture
    pairs=Iterators.product(orgs,stocks) |> collect 
    out = Base.splat(in).(pairs)
    return out' # output expected bo be cultures by strains instead of strains by cultures
end  










function dispense_solver(sources::Vector{T},destinations::Vector{U},robots::Vector{V},source_compatibility::BitMatrix,secondary_objectives...;quiet::Bool=true,timelimit::Real=10,pad::Real=1.25,slack_tolerance::Real=0,overdraft_tolerance::Real=1e-8,require_nonzero::Bool=true,return_model::Bool=false,obj_tolerance=1e-6,obj_cutoff=1e-3,inoculation_quantity::Real=2, priority::Dict{JLIMS.Chemical,UInt64}=Dict{JLIMS.Chemical,UInt64}(),kwargs...) where {T<: JLIMS.Culture,U<:JLIMS.Culture,V<:Robot}
    # check the length of the inputs 
    S= length(sources) 
    D= length(destinations) 
    R= length(robots)

    
    # check that the source compatibility matrix is valid  
    a,b=size(source_compatibility)
    (a,b) == (S,R) ? nothing : error("compatibility matrix size ($a x $b) does not agree with the stocks and robots ($(S) x $(R))")
    src_pairs= Iterators.product(sources,robots) |> collect # iterate all pairs of sources and robots 
    all( source_compatibility .<= Base.splat(is_compatible_source).(src_pairs)) ? nothing : throw(StockCompatibilityError("at least one source -- robot pair is incompatible"))


    # calculate destination compatibility 
    dest_pairs= Iterators.product(destinations,robots) |> collect
    destination_compatibility = Base.splat(is_compatible_destination).(dest_pairs) 
    at_least_one_robot=sum(destination_compatibility,dims=2) .>= 1 
    all(at_least_one_robot) ? nothing : error("at least one destination has no compatible robots")

    
    # check keyword parameters  for issues 
    pad >= 1 ? nothing : error("padding factor must be greater than or equal to 1")
    0 <= overdraft_tolerance ? nothing : error("overdraft tolerance must be nonnegative")
    0 <= slack_tolerance ? nothing : error("Slack tolerance must be nonnegative")

    # create an S x R matrix of the minimum and maximum shot values for each stock -- robot pair 
    minshots = zeros(S,R)
    maxshots = Inf * ones(S,R)
    for s in 1:S
        for r in 1:R 
            robot=robots[r]
            stock= sources[s].media  
            minshots[s,r] = ustrip(uconvert(u"µL",hasproperty(robot.properties,:minVol) ? robot.properties.minVol : 0u"µL"))
            maxshots[s,r]=ustrip(uconvert(u"µL",hasproperty(robot.properties,:maxVol) ? robot.properties.maxVol : Inf*u"µL"))
            if stock isa JLIMS.SolidStock
                minshots[s,r]=ustrip(uconvert(u"mg",hasproperty(robot.properties,:minMass) ? robot.properties.minMass : 0u"g"))
                maxshots[s,r]=ustrip(uconvert(u"mg",hasproperty(robot.properties,:maxMass) ? robot.properties.maxMass : Inf*u"g"))
            end 
        end 
    end 
            


    # Gather all ingredients contained in the sources, destinations, and priority list 
    source_ingredients= unique(vcat(map(x->ingredients(JLIMS.composition(x)),sources)...))
    destination_ingredients= unique(vcat(map(x->ingredients(JLIMS.composition(x)),destinations)...))
    all_ingredients = unique(vcat(source_ingredients,destination_ingredients,collect(keys(priority))))
    I= length(all_ingredients)

    # gather all strains contained in the sources and destinations 
    all_strains=JLIMS.Strain[]
    source_strains = unique(vcat(filter(y->!ismissing(y),map(x->x.strains,sources))...))
    destination_strains = unique(vcat(filter(y->!ismissing(y),map(x->x.strains,destinations))...))
    all_strains=push!(all_strains,source_strains...)
    all_strains=push!(all_strains,destination_strains...)
    Y=length(all_strains)
    src_strain_array=strain_array(sources,all_strains)

    dest_strain_array=strain_array(destinations,all_strains)

    # group all available_sources by location 
    src_labware=unique(map(x->JLIMS.well(x).locationid,sources))
    SL=length(src_labware)
    slw_idxs=[findall(x->JLIMS.well(x).locationid==src_labware[l],sources) for l in 1:SL ]


    # group all destinations by location 
    dest_labware=unique(map(x->well(x).locationid,destinations))
    DL=length(dest_labware)
    dlw_idxs=[findall(x->well(x).locationid==dest_labware[l],destinations) for l in 1:DL ]

    # get and convert the source stock quantities for overdraft constants 
    source_quantities = ustrip.(map(x->uconvert(preferred_quantity(x),JLIMS.quantity(x)),sources))

    # Create the concentration array for the sources and the target quantity array for the destinations
    source_concentrations=ingredient_array(sources,all_ingredients) 
    sc = Float64.(Matrix(Unitful.ustrip.(source_concentrations)))

    #create the target array for the destinations 
    destination_quantities=ingredient_array(destinations,all_ingredients;measure=quantity)
    dq = Float64.(Matrix(Unitful.ustrip.(destination_quantities)))


    # check for missing source ingredients needed to complete the destination
    missing_ingredients=filter(x-> sum(Unitful.ustrip.(source_concentrations[:,Symbol(x.name)]))==0,destination_ingredients)

    if length(missing_ingredients) > 0 
            ings=map(x->x.name,missing_ingredients)
            throw(MissingIngredientError("No valid source of $(join(ings,", ")) available to complete the dispenses",ings))
    end 

    # Update priorities: priority increases with decreasing level

    # Level 0: All level 0 ingredeients are blocked from the design. They may not appear in the stock. 
    # Level 1+: All other ingredients are scheduled sequentially by priority. we allow slack for all nonzero priorities, where we try to minimize the slack for each ingredient and then constrain the slack for future priority levels. 
    # Level 2^(64)-1: Maximum allowable priority level, typically a solvent like water that we would use to back fill will be assigned maximum priority level. 

    for ing in destination_ingredients
        if !in(ing,keys(priority))
            priority[ing]=UInt64(1)
        end 
    end 
    for ing in source_ingredients
        if !in(ing,keys(priority)) # if the user hasn't specified a source ingredient priority or included it in the design (see above), assume it is priority 0. (These ingredients will be blocked, regardless of a discrete or continuous robot)
            priority[ing]=UInt64(0) 
        end 
    end 

    capvals=source_quantities/pad  # the maximum quantity of each source 

    maxshot_param=Inf * ones(S,R+1)
    minshot_param=zeros(S,R+1)
    shot_density = ones(S,R+1) # assume that one shot delivers one unit of source 
    for r in 1:R
        if robots[r].properties.isdiscrete
            maxshot_param[:,r].=maxshots[:,r] ./ minshots[:,r] # convert all masses and volumes to shots for the discrete problem 
            minshot_param[:,r].=minshots[:,r] ./ minshots[:,r]
            shot_density[:,r] .= deepcopy(minshots[:,r]) # for discrete robots, one shot delivers the minimum shot value
        else 
            maxshot_param[:,r] .= maxshots[:,r]
            minshot_param[:,r].= minshots[:,r]
        end 
    end 

    maxshot_param[:,(R+1)] .= source_quantities
    minshot_param[:,(R+1)].= source_quantities

    # initialize the JuMP model 
    model=Model(Gurobi.Optimizer) 
    if quiet 
        set_silent(model)
    end 
    set_attribute(model,"TimeLimit",timelimit)
    M=R+1
    # Define constants 
    @variable(model, minShots[1:S,1:M] in Parameter.(minshot_param))
    @variable(model, maxShots[1:S,1:M] in Parameter.(maxshot_param))
    @variable(model, caps[1:S] in Parameter.(capvals)) # save a parameter for the volume of each available source

    # Define Model variables 
    @variable(model, shots[1:S,1:D,1:M]>=0) # q[s,d,r] = shots of source s transfered to destination d using robot r . We add an extra 'null' robot that "transfers" contents of a single well to itself
    @variable(model,q[1:S,1:D,1:M]>=0)
    @variable(model,qt[1:S,1:D]>=0)
    @variable(model, qi[1:S,1:D,1:M],Bin) #  Indicator for whether a transfer in q is active
    @variable(model, stri[1:S,1:D],Bin) # strain transfer indicator 
    @variable(model,slacks[1:D,1:I]) # slack variables to measure the difference between the dispenses and the target quantities for each ingredient of each destination
    @variable(model, lw[1:SL,1:DL]>=0) # continuous variable to measure the total quantity dispensed from and to each labware 
    @variable(model,strain_slacks[1:D,1:Y])
    for s in 1:S
        for d in 1:D 
            for r in 1:M
                @constraint(model, q[s,d,r] == shots[s,d,r]*shot_density[s,r])
                @constraint(model, !qi[s,d,r] => {shots[s,d,r]== 0}) # tie the transfer indicator to the transfers
            end 
        end 
    end 
    for s in 1:S 
        for d in 1:D 
            @constraint(model,qt[s,d] == sum(q[s,d,:]))
        end 
    end 

    for s in 1:S 
        for d in 1:D 
            @constraint(model, !stri[s,d] => {sum(src_strain_array[s,:])*sum(shots[s,d,:])==0})
        end 
    end 

    
    @constraint(model, qt'*sc .- slacks .== dq ) # create the destinations with the sources, allowing for some slack. This looks messy but it is a mass/volume balance. 

    @constraint(model, stri'*src_strain_array .== dest_strain_array) # ensure that only strains meant to be dispensed are dispensed 
    @constraint(model, qt'*src_strain_array .- strain_slacks .== inoculation_quantity*dest_strain_array) # ensure that if a strain is dispensed, meaure the discrepancy with the strain slack

    for r in 1:R 
        if robots[r].properties.isdiscrete 
        set_integer.(shots[:,:,r])  
        @constraint(model, shots[:,:,r] .>= qi[:,:,r]) # enforce activity constraint
        else 
            for s in 1:S 
                @constraint(model, shots[s,:,r] .>= minShots[s,r]*qi[s,:,r]) # if a dispense is active it must be larger than the minimum shot volume for the robot. This formulates q as a semicontinuous variable from (minshotval, Inf)
            end 
        end
    end
    @constraint(model,shots[:,:,M] .>= minShots[:,M].*qi[:,:,M])

    for s in 1:S 
        for d in 1:D 
            if JLIMS.well(sources[s]) == JLIMS.well(destinations[d]) # if the source well is the same as the destination well, the entire quantity of the source must be transferred to the destination
                @constraint(model,q[s,d,R+1]==source_quantities[s]) # the total quantity must be the source quantity
                @constraint(model,shots[s,:,1:R].==0) # all dispenses on other robots must be zero for the source
                @constraint(model,shots[s,setdiff(1:D,d),R+1].==0) # the `null` robot cannot be used for other destinations
            else 
                @constraint(model,shots[s,d,R+1]==0)  # the `null` robot cannot be used for the source-dispense pair
            end 
        end 
    end 
 

    for s in 1:S
        for d in 1:D
            for r in 1:R
                if !(source_compatibility[s,r] && destination_compatibility[d,r])
                    @constraint(model, shots[s,d,r]==0) # don't allow transfers with incompatible robots
                end
            end 
        end 
    end 
    for y in 1:Y 
        for d in 1:D
            if dest_strain_array[d,y] == false 
                @constraint(model, strain_slacks[d,y]==0)
            end 
        end 
    end 
    for i in 1:I
        if priority[all_ingredients[i]] == UInt(0)
            for d in 1:D
                @constraint(model, slacks[d,i]==0) # priority 0 ingredients must hit the target exactly for each destination. The slack must be zero (because the delivered quantity must be zero)  
            end 
        end 
    end 


    if require_nonzero 
        for i in 1:I
            sources_with_ingredient=sc[:,i] .> 0 
            for d in 1:D
                if dq[d,i] > 0 
                    @constraint(model,sum(qi[:,d,:] .* sources_with_ingredient) >=1) # check that at least one transfer happens if an ingredient is needed in a destination, even if the optimial solution is to not dispense anything. 
                end 
            end 
        end 
    end 
     

    for sl in 1:SL 
        for dl in 1:DL
            @constraint(model, lw[sl,dl]==sum(qt[slw_idxs[sl],dlw_idxs[dl]])) # constrain the labware volume to track the dispenses from each source location to each destination location 
        end
    end 


    set_objective_sense(model, MOI.FEASIBILITY_SENSE)
    #target_bound=sum((inoculation_quantity*dest_strain_array).^2)
    #set_attribute(model,"Cutoff",obj_cutoff* target_bound) # reject any solutions that are larger than the objective tolerance, which is a percentage of the sum squared target quantities.
    #set_attribute(model, "BestObjStop",obj_tolerance*target_bound) 
    @objective(model, Min,sum(strain_slacks.^2)) # penalize large slacks, search for q that minimize slacks.
    optimize!(model)
    println(termination_status(model)) 
    current_strain_slacks=abs.(JuMP.value.(strain_slacks))
    for y in 1:Y
            delta=slack_tolerance * inoculation_quantity* dest_strain_array[:,y] # delta is the tolerance we give to updating the slack in higher priority levels, it is some fraction of the dispense target quantity for every destination. Wiggle room for the slack in future iterations 
            current_strain_slack=current_strain_slacks[:,y]
            @constraint(model, strain_slacks[:,y] .>= -current_strain_slack .- delta)
            @constraint(model, strain_slacks[:,y] .<= current_strain_slack .+ delta)
    end 
    #reset the optimizer to remove the objective cutoff since we are switching objectives 
    set_optimizer(model,Gurobi.Optimizer)
    set_attribute(model,"TimeLimit",timelimit)


    priority_levels= sort(unique(collect(values(priority))))
    for level in priority_levels  # pass through all priority levels from lowest to higest level


        weights= falses(I)
        for i in 1:I
            if priority[all_ingredients[i]] <= level 
                weights[i]=true # activate the weight term for this ingredient on this pass
            end 
        end 
        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        target_bound=sum(weights'.*dq.^2)
        set_attribute(model,"Cutoff",obj_cutoff* target_bound) # reject any solutions that are larger than the objective tolerance, which is a percentage of the sum squared target quantities.
        set_attribute(model, "BestObjStop",obj_tolerance*target_bound) 
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
            println(term)
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
                open("/Users/BDavid/Desktop/iis.txt","w") do file
                    for con in out_cons
                        println(file,con)
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
    quants=JuMP.value.(q)[:,:,1:R]

    
    

    sources_needed=sum(quants,dims=2) # the total quantity of each source across all robots and destinations.

    overdrafts = sources_needed .- cap_vals
    if any(overdrafts .>overdraft_tolerance) 
        overdraft_dict=Dict{JLIMS.Culture,Unitful.Quantity}()
        for i in findall(x-> x > 0 ,overdrafts)
            od=overdrafts[i]
            overdraft_dict[sources[i]] = od * preferred_quantity(sources[i])
        end 
        throw(OverdraftError("Refills are needed for $(length(collect(keys(overdraft_dict)))) sources:",overdraft_dict))
    end 
    transfers=DataFrame[]   
    
    for robot in robots 
        r = findfirst(x->x==robot,robots)
        tf = DataFrame()
        for dest in destinations
            quantities=Any[]
            d=findfirst(x->x==dest,destinations)
            for src in sources
                un=preferred_quantity(src)
                    s=findfirst(x->x==src,sources)
                    val=quants[s,d,r]
                    push!(quantities,val*un) # commit the transfer to the transfers vector for this destination stock 
            end 
            tf[:,Symbol("Well$(well(dest).id)")]=quantities
        end
        push!(transfers,tf)
    end  
    println("Solution off target by $(sum(JuMP.value.(slacks).^2)/sum(dq.^2)*100)%")

    if return_model
        return transfers,model
    else 
        return transfers
    end 
end 
    
function dispense_solver(sources::Any,destinations::Any,robots::Vector{V},secondary_objectives...;kwargs...) where  V<:Robot
        # calculate destination compatibility 
        src_pairs= Iterators.product(sources,robots) |> collect
        source_compatibility = Base.splat(is_compatible_source).(src_pairs) 
        return dispense_solver(sources,destinations,robots,source_compatibility,secondary_objectives...;kwargs...)
end 

function dispense_solver(sources::Vector{T},destinations::Any,robots::Vector{V},source_compatibility::BitMatrix,secondary_objectives...;kwargs...) where {T<:JLIMS.Stock,V<:Robot}
    srcs = convert.((JLIMS.Culture,),sources)
    return dispense_solver(srcs,destinations,robots,source_compatibility,secondary_objectives...;kwargs...)
end 

function dispense_solver(sources,destinations::Vector{U},robots::Vector{V},source_compatibility::BitMatrix,secondary_objectives...;kwargs...) where {U<:JLIMS.Stock,V<:Robot}
    dests = convert.((JLIMS.Culture,),destinations)
    return dispense_solver(sources,dests,robots,source_compatibility,secondary_objectives...;kwargs...)
end 

function dispense_solver(sources::Vector{T},destinations::Vector{T},robots::Vector{V},source_compatibility,secondary_objectives...;kwargs...) where {T<:JLIMS.Stock,V<:Robot}
    srcs = convert.((JLIMS.Culture,),sources)
    dests = convert.((JLIMS.Culture,),destinations)
    return dispense_solver(srcs,dests,robots,source_compatibility,secondary_objectives...;kwargs...)
end







#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
t,m=dispense_solver(sources,destinations,cobra_default;return_model=true)
t,m=dispense_solver(sources,destinations,mantis_default;return_model=true)
=#
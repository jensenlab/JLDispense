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

function feasibility(sources::Vector{T},destinations::Vector{U},robot::Robot;quiet=true,timelimit=30,pad=1.25) where {T <:JLIMS.Stock,U<:JLIMS.Stock} 
    # determine if the destinations can be made from the sources on this particular robot, if so,  find the set of sources that use the fewest pieces of LABWARE 
    available_sources=deepcopy(sources)
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
    
        available_sources= filter(x->in(x.well.container,compatible_containers),sources) 
    end 


    available_sources= filter(x->in(typeof(x),robot.properties.compatible_stocks),available_sources)

    destination_ingredients=unique(vcat(ingredients.(map(x->x.composition,destinations))...))
    available_sources=filter(y->all(map(x->in(x,destination_ingredients),ingredients(y.composition))),available_sources) # block any source that has extraneous ingredients 


    source_labware=unique(map(x->x.well.labwareid,available_sources))
    sl_idx=map(x->findall(y->y.well.labwareid==x,available_sources),source_labware)
    #source_containers=map(x->available_sources[x[1]].well.container,sl_idx)
  
    destination_labware=unique(map(x->x.well.labwareid,destinations))
    dl_idx=map(x->findall(y->y.well.labwareid==x,destinations),destination_labware)
    #destination_containers=map(x->destinations[x[1]].well.container,dl_idx)
    sl=length(source_labware)
    dl=length(destination_labware)

    
    target_quantities=stock_quantity_array(destinations)
    

    nt=DataFrames.names(target_quantities)
    
    concentrations=stock_concentration_array(available_sources)
    needed_ingredients=map(x->!in(x,DataFrames.names(concentrations)),nt)
    if any(needed_ingredients)
            ings=filter(x->!in(x,DataFrames.names(concentrations)),nt)
            throw(MissingIngredientError("No valid source of $(join(ings,",")) available to complete the dispenses",ings))
    end 

    concentrations=concentrations[:,nt]
    source_quantities= map(x->x.quantity,available_sources) 
    source_units= preferred_stock_quantity.(available_sources)
    source_quants=uconvert.(source_units,source_quantities)
    concentrations=Float64.(ustrip.(Matrix(concentrations)))
  
    target_quantities=Float64.(ustrip.(Matrix(target_quantities)))
    

    source_quants=ustrip.(source_quants)


    S,I=size(concentrations)
    D,I=size(target_quantities)
    
    model=Model(Gurobi.Optimizer) 
        if quiet 
            set_silent(model)
        end 
            set_attribute(model,"TimeLimit",timelimit)
  
        
            @variable(model,q[1:S,1:D] >=0)
            @variable(model,lw[1:sl],Bin)
            @constraint(model, q'*concentrations .== target_quantities) # find the quantity of each stock needed to make the target 



        for s in 1:sl 
            @constraint(model, !lw[s]=>{sum(q[sl_idx[s],:]) == 0})
        end 
         
        @objective(model, Min,sum(lw)) # minimize the number of labware needed to complete the dispenses  
        optimize!(model) # optimize to check that the source stocks could in theory make the destinations 
        if !JuMP.is_solved_and_feasible(model)
            throw(MixingError("the requested stocks cannot be made using this combination of sources"))
        end

        quants=pad*sum(JuMP.value.(q);dims=2) # find out how much we need 

        
        needed =Dict(available_sources .=> quants.*source_units)
        
        for s in 1:S
        @constraint(model, pad*sum(q[s,:]) <= source_quants[s])
        end 
        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        @objective(model, Min,sum(lw))
        optimize!(model)
        if !JuMP.is_solved_and_feasible(model)
            throw(OverdraftError("Refills are needed of the following stocks",needed))
        end
        println(JuMP.value.(q))
        return  true 

end 
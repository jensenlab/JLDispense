struct MixingError <:Exception
        msg::AbstractString
        constraints::Vector{JuMP.ConstraintRef}
     end
        
    struct OverdraftError <: Exception 
        msg::AbstractString
        balances::Dict{JLIMS.Stock,Unitful.Quantity}
    end 
    
    struct MissingChemicalError <: Exception 
        msg::AbstractString
    end 

    struct MissingOrganismError <: Exception 
        msg::AbstractString
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



function dispense_solver(sources::Vector{Well},targets::Vector{Well},configurations::Vector{<:Configuration};quiet::Bool=true,timelimit::Real=10,pad::Real=1.25,slack_tolerance::Real=0,overdraft_tolerance::Real=1e-8,require_nonzero::Bool=true,return_model::Bool=false,obj_tolerance=1e-2,obj_cutoff=1e-3,inoculation_quantity::Real=2, priority::Dict{JLIMS.Chemical,UInt64}=Dict{JLIMS.Chemical,UInt64}(),kwargs...) 

        all_labware = get_all_labware(sources,targets)

        all_chemicals =get_all_chemicals(sources,targets,priority) 

        all_organisms = get_all_organisms(sources,targets) 

        src_stocks, tgt_stocks, s_enforced,t_enforced= slot_stocks(sources,targets)

        src_quantities =ustrip.(map(x->uconvert(preferred_quantity(x),JLIMS.quantity(x)),src_stocks))

        configs = vcat(configurations,null_robot) # add in the null robot to handle the "hidden" transfers that go directly from a well to itself
        
        # Create the concentration array for the sources and the target quantity array for the destinations
        source_concentrations=chemical_array(src_stocks,all_ingredients) 
        sc = Float64.(Matrix(Unitful.ustrip.(source_concentrations)))
        so = strain_array(src_stocks,all_organisms)
        #create the target array for the destinations 
        target_quantities=chemical_array(destinations,all_ingredients;measure=quantity)
        tq = Float64.(Matrix(Unitful.ustrip.(target_quantities)))
        to = strain_array(tgt_stocks,all_organisms)

 
        # check for missing source ingredients needed to complete the destinations
        missing_chemicals=filter(x-> sum(Unitful.ustrip.(source_concentrations[:,Symbol(x.name)]))==0,get_all_chemicals(targets))       
        if length(missing_chemicals) > 0 
                throw(MissingChemicalError("No valid source of $(join(name.(missing_chemicals),", ")) available to complete the dispenses",missing_chemicals))
        end 

        # check for missing source organisms needed to create the targets 
        missing_organisms = filter(x-> sum(so[:,Symbol(JLIMS.name(x))])==0,get_all_organsims(targets))
        if length(missing_organisms)>0 
                throw(MissingOrgansimError("No valid source of $(join(name.(missing_organisms),", ")) available to complete the dispenses",missing_organisms))
        end

        # update the chemical priority dict to account for the sources and targets 
        update_priority!(sources,targets,priority)


        ######
        # start defining model 
        ######
        model=Model(Gurobi.Optimizer) 
        if quiet 
            set_silent(model)
        end 
        set_attribute(model,"TimeLimit",timelimit)
        W = length(src_stocks) # number of total wells in the model
        C = length(all_chemicals)
        O = length(all_organisms)
        L = length(all_labware)
        R = length(configs)
        lw_product = Iterators.product(all_labware,all_labware) # build a pairwise product of every labware to iterate over 
        lw_idxs = Iterators.product(1:L,1:L)
    
        # Define Model variables 
        @variable(model, Qw[1:W,1:W]>=0) # q(i,j) = quantity of material transferred from well i to well j 
        @variable(model,Qwa[1:W,1:W]>=0)
        Q = [@variable(model, [1:length(la),1:length(ld)],lower_bound=0) for (la,ld) in lw_product]
        Qa = [@variable(model, [1:length(la),1:length(ld)],lower_bound=0) for (la,ld) in lw_product]
        @variable(model,chem_slacks[1:W,1:C])
        @variable(model,org_slacks[1:W,1:O])
        for (aa,dd) in lw_idxs
                for (i,j) in CartesianIndices(Q[aa,dd])
                        @constraint(model,Qw[get_well_index(all_labware,aa,i),get_well_index(all_labware,dd,j)]== Q[aa,dd][i,j] ) # tie the wells segregated by labware to the WxW variable. Having both makes it easier to work with each type when indexing
                        @constraint(model,Qwa[get_well_index(all_labware,aa,i),get_well_index(all_labware,dd,j)]== Qa[aa,dd][i,j] ) 
                end
        end
        # constraints to ensure that we create the targets properly
        @constraint(model, Qw'*sc .- chem_slacks .== tq ) # create the targets with the sources, allowing for some slack. This is a mass/volume balance. 
        @constraint(model, Qw'*so .- org_slacks .== inoculation_quantity*to) # ensure that if a strain is dispensed, meaure the discrepancy with the strain slack
        @constraint(model,org_slacks .== 0 ) # organisms must be hit exactly. They trump all other priorities. 
    


        # measure the physical masked operations 
        Vd = [write_masked_operation_variables(model,c,la,ld,false) for  (la,ld) in lw_product, c in configs] # stores the quantity of material moved from position 
        Va = [write_masked_operation_variables(model,c,la,ld,true) for  (la,ld) in lw_product, c in configs]
        e = length(configs)
        for (aa,dd) in lw_idxs
                for c in eachindex(configs) 
                        cha= size(Va[aa,dd,c])[2] # aspirate channel id
                        for ch in 1:cha
                        @constraint(model, Va[aa,dd,c][:,ch,:,:] .>= pads(head(c))[ch]* Vd[aa,dd,c][:,ch,:,:] .+ deadvols(head(configs[c]))[ch] ) # add the appropriate padding factor and deadvolume of the aspirating channel to the aspirating volume 
                        end 
                end
                for (i,j) in CartesianIndices(Q[aa,dd])
                        if i != j 
                                @constraint(model, sum([sum([masks(head(configs[c]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[c]),all_labware[dd][2](i,pd,chd))*Vd[aa,dd,c][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,c])]) for c in eachindex(configs)[1:end]]) == Q[aa,dd][i,j]) # can use any of the real robots to do transfers
                                @constraint(model, sum([masks(head(configs[e]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[e]),all_labware[dd][2](i,pd,chd))*Vd[aa,dd,e][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,e])])  == 0) # null robot must not be used 
                                @constraint(model, sum([sum([masks(head(configs[c]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[c]),all_labware[dd][2](i,pd,chd))*Va[aa,dd,c][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,c])]) for c in eachindex(configs)[1:end]]) == Qa[aa,dd][i,j]) # can use any of the real robots to do transfers
                                @constraint(model, sum([masks(head(configs[e]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[e]),all_labware[dd][2](i,pd,chd))*Va[aa,dd,e][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,e])])  == 0) # null robot must not be used 
                                # This expression sums over all robot configurations c in the outer loop, and then sums over all mask position pairs (pa,ca,pd,cd) (p = position ,ch= channel) for the given labware and configuration pairings. 
                                # Ultimately, this is the constraint that converts masked operations into flows from well to well
                                # the main calculation is to multiply the aspirate and dispense masks Ma =masks()[1] and Md masks()[2]  by the dispensed volume flow Vd. 
                                # only flows that pass both masks count, the rest are irrelevant. 
                        else 
                                @constraint(model, sum([sum([masks(head(configs[c]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[c]),all_labware[dd][2](i,pd,chd))*Vd[aa,dd,c][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,c])]) for c in eachindex(configs)[1:end-1]]) == 0) # cannot use any of the real robots to do transfers, only the null robot
                                @constraint(model, sum([masks(head(configs[e]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[e]),all_labware[dd][2](i,pd,chd))*Vd[aa,dd,e][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,e])])  == Q[aa,dd][i,j]) # null robot must pick up all slack 
                                @constraint(model, sum([sum([masks(head(configs[c]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[c]),all_labware[dd][2](i,pd,chd))*Va[aa,dd,c][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,c])]) for c in eachindex(configs)[1:end-1]]) == 0) # cannot use any of the real robots to do transfers, only the null robot
                                @constraint(model, sum([masks(head(configs[e]),all_labware[aa])[1](i,pa,cha)*masks(head(configs[e]),all_labware[dd][2](i,pd,chd))*Va[aa,dd,e][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(Va[aa,dd,e])])  == Qa[aa,dd][i,j]) # null robot must pick up all slack 
                        end 
                end 
 
        end 

 
        for c in 1:C
                if priority[all_chemicals[c]] == UInt(0)
                        for d in 1:W
                                @constraint(model, chem_slacks[d,c]==0) # priority 0 ingredients must hit the target exactly for each destination. The slack must be zero (because the delivered quantity must be zero)  
                        end 
                end 
        end 


        if require_nonzero  # requiring nonzero for each needed destination well and chemical necesitates a Binary variable, which turns the problem into an MILP. 
                @variable(model,QwI[1:W,1:W],Bin) 
                for i in 1:W 
                        for j in 1:W
                                @constraint(model, !QwI[i,j] => {Q[i,j]== 0}) # tie the transfer indicator to the transfers
                        end
                end
                for c in 1:C
                sources_with_chemical=sc[:,c] .> 0 
                for d in 1:W
                        if tq[d,c] > 0 
                        @constraint(model,sum( QwI[:,d] .* sources_with_chemical) >= 1  ) # check that at least one transfer happens if an ingredient is needed in a destination, even if the optimial solution is to not dispense anything.  
                        end 
                end 
                end 
        end 

        ############################################################################################################################################################
        # SET UP OBJECTIVES

        ############################################################################################################################################################

        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        #target_bound=sum((inoculation_quantity*dest_strain_array).^2)
        #set_attribute(model,"Cutoff",obj_cutoff* target_bound) # reject any solutions that are larger than the objective tolerance, which is a percentage of the sum squared target quantities.
        #set_attribute(model, "BestObjStop",obj_tolerance*target_bound) 
        #=
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
                =# 
        target_weights = t_enforced 
        priority_levels= sort(unique(collect(values(priority))))
        for level in priority_levels  # pass through all priority levels from lowest to higest level


                chem_weights= falses(C)
                for c in 1:C
                        if priority[all_chemcials[c]] <= level 
                                chem_weights[i]=true # activate the weight term for this ingredient on this pass
                        end 
                end 
                weights = target_weights .* chem_weights' 
                set_objective_sense(model, MOI.FEASIBILITY_SENSE)
                #target_bound=sum(weights .*tq.^2) # activate only the slacks for ingredients with this priority level 
                #set_attribute(model,"Cutoff",obj_cutoff* target_bound) # reject any solutions that are larger than the objective tolerance, which is a percentage of the sum squared target quantities.
                #set_attribute(model, "BestObjStop",obj_tolerance*target_bound) 
                @objective(model, Min,sum(weights .*(chem_slacks.^2))) # penalize large slacks, search for Qw that minimize slacks.
                optimize!(model)

                # Check for optimality and feasibility 
                term = termination_status(model)
                if term == MOI.OPTIMAL || term == MOI.OBJECTIVE_LIMIT
                if primal_status(model)==MOI.FEASIBLE_POINT
                        println("Optimal Solution Found For Level $level")
                elseif primal_status(model)==MOI.NO_SOLUTION
                        throw(error("No solution exists For Level $level"))
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
                #=
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
                =#
                end 
                current_slacks=abs.(JuMP.value.(chem_slacks))
                
                for c in 1:C
                if chem_weights[c]
                        delta=slack_tolerance * tq[:,c] # delta is the tolerance we give to updating the slack in higher priority levels, it is some fraction of the dispense target quantity for every destination. Wiggle room for the slack in future iterations 
                        current_slack=current_slacks[:,c]
                        @constraint(model, chem_slacks[:,c] .>= -current_slack .- delta)
                        @constraint(model, chemslacks[:,c] .<= current_slack .+ delta)
                end 
                end 
        
                
        end 
        #reset the optimizer to remove the objective cutoff since we are switching objectives 
        set_optimizer(model,Gurobi.Optimizer)
        set_attribute(model,"TimeLimit",timelimit)
        optimize!(model) # resolve one last time with the final slack constraints -> we need to optimize before querying results for the secondary objectives 

        set_objective_sense(model, MOI.FEASIBILITY_SENSE)


        @objective(model, Min,sum(sum.(Va))+sum(sum.(Vd))) # minimize total masked operations 
        optimize!(model)

        # check for overdrafts 
        qa_quants = JuMP.value.(Qwa)
        quants_needed=sum(Qwa,dims=2) # the total quantity of each source across all robots and destinations.

        quants_needed=quants_needs .* s_enforced * u"µL"

        overdrafts = quants_needed .- map(x->uconvert(preferred_quantity_unit(x),quantity(x),stock.(sources)))
        if any(overdrafts .>overdraft_tolerance) 
                overdraft_dict=Dict{JLIMS.Well,Unitful.Quantity}()
                for i in findall(x-> x > 0 ,overdrafts)
                        od=overdrafts[i]
                        overdraft_dict[sources[i]] = od 
                end 
                throw(OverdraftError("Refills are needed for $(length(collect(keys(overdraft_dict)))) sources:",overdraft_dict))
        end 


        # check for solution optimality 

        c_slacks =JuMP.value.(chem_slacks)
        sol_quality=zeros(size(c_slacks))
        for (d,c) in CartesianIndices(c_slacks)
                slack= c_slacks[d,c]^2 
                target = tq[d,c]^2 
                if target == 0 && slack == 0 
                        sol_quality[d,c]= 0 
                elseif target == 0 && slack != 0 
                        error("stock target ($d), chemical $(all_chemcials[c]), unbounded error")
                else
                        sol_quality[d,c]= slack/target 
                end 
        end 
        if max(sol_quality) > obj_tolerance
                error("unacceptable solution: \n max error = $(max(sol_quality)*100)% \n tolerance limit = $(obj_tolerance*100)%")
        end 



        dispense_ops = map(x->JumP.value.(x),Vd)
        slotting_indicators = map(x-> sum(x) > 0, dispense_ops) 

        config_dispenses=Matrix{Real}[]

        for con in eachindex(configs)[1:end-1] 
                ops = dispense_ops[:,:,con]
                Q_config= zeros(W,W) 
                for la_idx in 1:L
                        la= all_labware[la_idx]
                        Ma,p1,p2,p3=masks(configs[con],la)
                        for ld in 1:L
                                ld=all_labware[ld_idx]
                                p1,Md,p2,p3=masks(configs[con],ld)
                                for ida in eachindex(children(la))
                                        i = get_well_index(la,la_idx,ida)
                                        for idd in eachindex(children(ld))
                                                j = get_well_index(ld,findfirst(x->x==ld,all_labware),idd)
                                                Q_config[i,j]=sum([Ma(i,pa,cha)*Md(j,pd,chd)*ops[la_idx,ld_idx][pa,cha,pd,chd] for (pa,cha,pd,chd) in CartesianIndices(ops[la_idx,ld_idx])])
                                        end
                                end
                        end
                end

                push!(config_dispenses,Q_config)
        end 




        if return_model 
                return config_dispenses,slotting_indicators,model 
        else
                return config_dispense,slotting_indicators
        end

end



function get_well_index(all_labware,lw_index,well_index)
        start = sum(x->length(x),all_labware[1:(lw_index-1)])
        return start+well_index
end 



function write_masked_operation_variables(model,config,l_asp,l_disp)
        a,b,Sa,c=masks(head(config),l_asp) # a,b,c are irrelevant for this function 
        a,b,c,Sd=masks(head(config),l_disp)
        return @variable(model, [1:Sa[1],1:Sa[2],1:Sd[1],1:Sd[2]],lower_bound= 0 ) # all of the positions of the mask for both the aspirate and dispense labware 
end 
 
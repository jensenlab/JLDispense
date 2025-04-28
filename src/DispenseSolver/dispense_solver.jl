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



function dispense_solver(sources::Vector{<:Well},targets::Vector{<:Well},configurations::Vector{<:Configuration};quiet::Bool=true,timelimit::Real=10,pad::Real=1.25,slack_tolerance::Real=0,overdraft_tolerance::Real=1e-8,require_nonzero::Bool=false,return_model::Bool=false,obj_tolerance=1e-2,obj_cutoff=1e-3,inoculation_quantity::Real=2, priority::Dict{JLIMS.Chemical,UInt64}=Dict{JLIMS.Chemical,UInt64}(),kwargs...) 

        all_labware = get_all_labware(sources,targets)

        all_chemicals =get_all_chemicals(stock.(sources),stock.(targets),priority) 

        all_organisms = get_all_organisms(stock.(sources),stock.(targets)) 

        src_stocks, tgt_stocks, s_enforced,t_enforced,tgt_caps = slot_stocks(sources,targets)

        src_quantities =ustrip.(map(x->uconvert(preferred_quantity_unit(x),JLIMS.quantity(x)),src_stocks))

        configs = vcat(configurations,null_robot) # add in the null robot to handle the "hidden" transfers that go directly from a well to itself
        
        # Create the concentration array for the sources and the target quantity array for the destinations
        source_concentrations=chemical_array(src_stocks,all_chemicals) 
        sc = Float64.(Matrix(Unitful.ustrip.(source_concentrations)))
        source_organisms = organism_array(src_stocks,all_organisms)
        so = Matrix(source_organisms)
        #create the target array for the destinations 
        target_quantities=chemical_array(tgt_stocks,all_chemicals;measure=quantity)
        tq = Float64.(Matrix(Unitful.ustrip.(target_quantities)))
        target_organisms = organism_array(tgt_stocks,all_organisms)
        to= Matrix(target_organisms)
 
        # check for missing source chemicals needed to complete the destinations
        missing_chemicals=filter(x-> sum(Unitful.ustrip.(source_concentrations[:,Symbol(x.name)]))==0,get_all_chemicals(stock.(targets)))       
        if length(missing_chemicals) > 0 
                throw(MissingChemicalError("No valid source of $(join(name.(missing_chemicals),", ")) available to complete the dispenses",missing_chemicals))
        end 

        # check for missing source organisms needed to create the targets 
        missing_organisms = filter(x-> sum(source_organisms[:,Symbol(JLIMS.name(x))])==0,get_all_organisms(stock.(targets)))
        if length(missing_organisms)>0 
                throw(MissingOrgansimError("No valid source of $(join(name.(missing_organisms),", ")) available to complete the dispenses",missing_organisms))
        end

        # update the chemical priority dict to account for the sources and targets 
        update_priority!(sources,targets,priority)

        W = length(src_stocks) # number of total wells in the model
        C = length(all_chemicals)
        O = length(all_organisms)
        L = length(all_labware)
        R = length(configs)


        Ma,Md,Sa,Sd = get_masks(configs,all_labware)
        lw_idx = vcat([fill(i,length(all_labware[i])) for i in 1:L]...)
        within_labware_index = vcat([1:length(all_labware[i]) for i in 1:L]...)

        MDict = Dict{Tuple{Int,Int,Int,Int},Int}()
        NDict = Dict{Tuple{Int,Int,Int,Int},Int}()
        m = 1 
        n=1

        Rm = []
        Rn = [] 
        mnozzles = Nozzle[]
        for r in 1:R 
                rmstarting= m
                rnstarting =n
                for l in 1:L 
                        Pa,Cha = Sa[r,l]
                        Pd,Chd = Sd[r,l]
                        for p in 1:Pa
                                for c in 1:Cha 
                                        MDict[(r,l,p,c)]=m 
                                        m += 1 
                                        push!(mnozzles,channels(head(configs[r]))[c])
                                end 
                        end 
                        for p in 1:Pd
                                for c in 1:Chd 
                                        NDict[(r,l,p,c)]= n 
                                        n += 1
                                end 
                        end 
                end 
                rmending = m-1 
                rnending = n-1
                push!(Rm,rmstarting:rmending)
                push!(Rn,rnstarting:rnending)


        end 
        M = length(keys(MDict))
        N = length(keys(NDict))


                        



        w_to_m = []
        w_to_n = []
        for w in 1:W 
                l = lw_idx[w]
                i = within_labware_index[w]
                wm_mapping = Int[]
                wn_mapping = Int[]
                for r in 1:R 
                        
                        ma = Ma[r,l]
                        md = Md[r,l]
                        Pa,Cha = Sa[r,l]
                        Pd,Chd = Sd[r,l]
                        for p in 1:Pa 
                                for c in 1:Cha  
                                        if ma(i,p,c) 
                                                push!(wm_mapping,MDict[(r,l,p,c)])
                                        end

                                end 
                        end
                        for p in 1:Pd 
                                for c in 1:Chd 
                                        if md(i,p,c)
                                                push!(wn_mapping,NDict[(r,l,p,c)])
                                        end 
                                end 
                        end  
                end 
                push!(w_to_m,wm_mapping)
                push!(w_to_n,wn_mapping)
        end


        ######
        # start defining model 
        ######
        model=Model(Gurobi.Optimizer) 
        if quiet 
            set_silent(model)
        end 
        JuMP.set_attribute(model,"TimeLimit",timelimit)

    
        # Define Model variables 
        @variable(model, Q[1:W,1:W]>=0) # q(i,j) = quantity of material transferred from well i to well j 
        @variable(model, V[1:M,1:N]>=0) # flow decision from m to n 
        @variable(model, Qa[1:W]>=0)
        @variable(model, Qd[1:W]>= 0 )
        @variable(model, PadLoss[1:M,1:R])
        @variable(model,chem_slacks[1:W,1:C])
        @variable(model,org_slacks[1:W,1:O])

        # constraints to ensure that we create the targets properly
        @constraint(model, Q'*sc .- chem_slacks .== tq ) # create the targets with the sources, allowing for some slack. This is a mass/volume balance. 
        @constraint(model, Q'*so .- org_slacks .== inoculation_quantity*to) # ensure that if a strain is dispensed, meaure the discrepancy with the strain slack
        #@constraint(model,org_slacks .== 0 ) # organisms must be hit exactly. They trump all other priorities. 
    


        # measure the physical masked operations 
        for i in 1:W 
                for j in 1:W 
                        idxs = Iterators.product(w_to_m[i],w_to_n[j])
                        idxs= CartesianIndex.(idxs)
                        @constraint(model, sum(V[idxs]) == Q[i,j] )
                end 
        end 

        for i in 1:W 
                if s_enforced[i]
                        ms = w_to_m[i]
                        nozzles = mnozzles[ms]
                        dead = map(x->ustrip(uconvert(u"µL",x.deadVol)),nozzles)
                        pad = map(x->x.deadVolFactor,nozzles)
                        for n in eachindex(nozzles) 
                                if nozzles[n] isa ContinuousNozzle
                                        pad[n] = pad[n] + dead[n]/ustrip(uconvert(u"µL",nozzles[n].maxAsp)) # the deadvolume should be a factor of the total aspiration for the flow for a continous nozzle
                                        dead[n]=0 # because its factored into pad 
                                
                                end
                        end
                        @constraint(model, sum( pad .* V[ms,:] .+ dead ) <= src_quantities[i]) # aspiration constraints 
                end
                if t_enforced[i]
                        ns = w_to_n[i]
                        @constraint(model, sum(V[:]) <= tgt_caps[i])
                end 
        end 


        for ra in 1:R 
                for rd in 1:R 
                        if ra != rd 
                                @constraint(model, V[Rm[ra],Rn[rd]].==0 ) # Flows cant go from one robot to another, they must stay within a particular robot 
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

        @objective(model, Min , sum(org_slacks.^2)) 
        optimize!(model) 

        current_slacks=abs.(JuMP.value.(org_slacks))
                 
        delta=slack_tolerance * inoculation_quantity*to[:,] # delta is the tolerance we give to updating the slack in higher priority levels, it is some fraction of the dispense target quantity for every destination. Wiggle room for the slack in future iterations 
        @constraint(model, org_slacks .>= -current_slacks .- delta)
        @constraint(model, org_slacks .<= current_slacks .+ delta)
        
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
                        if priority[all_chemicals[c]] <= level 
                                chem_weights[c]=true # activate the weight term for this ingredient on this pass
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
                                @constraint(model, chem_slacks[:,c] .<= current_slack .+ delta)
                        end 
                end 
        
                
        end 
        #reset the optimizer to remove the objective cutoff since we are switching objectives 
        set_optimizer(model,Gurobi.Optimizer)
        set_attribute(model,"TimeLimit",timelimit)
        optimize!(model) # resolve one last time with the final slack constraints -> we need to optimize before querying results for the secondary objectives 

        set_objective_sense(model, MOI.FEASIBILITY_SENSE)


        @objective(model, Min,sum(V)) # minimize total masked operations 
        optimize!(model)
        #=
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
        =#

        # check for solution optimality 

        c_slacks =JuMP.value.(chem_slacks)
        sol_quality=zeros(size(c_slacks))
        for (d,c) in CartesianIndices(c_slacks)
                slack= c_slacks[d,c]^2 
                target = tq[d,c]^2 
                if target == 0 && slack == 0 
                                sol_quality[d,c]= 0 
                elseif target == 0 && slack != 0 
                                error("stock target ($d), chemical $(all_chemicals[c]), unbounded error")
                else
                                sol_quality[d,c]= slack/target 
                end 
        end 
        if max(sol_quality) > obj_tolerance
                error("unacceptable solution: \n max error = $(max(sol_quality)*100)% \n tolerance limit = $(obj_tolerance*100)%")
        end 

        return model 

        dispense_ops = map(x->JumP.value.(x),V)
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
        start =0 
        if lw_index > 1 

                start = sum(map(x->length(x),all_labware[1:(lw_index-1)]))
        end
        return start+well_index
end 






function write_masked_operation_variables!(model,config,l_asp,l_disp)
        a,b,Sa,c=masks(head(config),l_asp) # a,b,c are irrelevant for this function 
        a,b,c,Sd=masks(head(config),l_disp)
        return @variable(model, [1:Sa[1],1:Sa[2],1:Sd[1],1:Sd[2]],lower_bound= 0 ) # all of the positions of the mask for both the aspirate and dispense labware 
end 
 

function get_masks(configs::Vector{<:Configuration},all_labware::Vector{Labware})
        
        M = map(I -> masks(I...),Iterators.product(head.(configs),all_labware))

        Ma,Md,Sa,Sd = map(x->getindex.(M,x),1:4)
        return Ma,Md,Sa,Sd
end 

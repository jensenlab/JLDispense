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



function dispense_solver(sources::Vector{<:Well},targets::Vector{<:Well},configs::Vector{<:Configuration},secondary_objectives...;robot_cost::Vector{<:Real}=ones(length(configs)),quiet::Bool=true, timelimit::Real=10,slack_tolerance::Real=1e-4,numerical_tolerance::Real=1e-8,require_nonzero::Bool=true,return_model::Bool=false,obj_tolerance=1e-3,inoculation_quantity::Real=2, priority::PriorityDict=Dict{JLIMS.Chemical,UInt64}(),round_digits::Int=1,kwargs...) 



        length(configs) == length(robot_cost) || ArgumentError("there must be a cost specified for each instrument configuration")

        all_labware = get_all_labware(sources,targets)

        all_chemicals =get_all_chemicals(stock.(sources),stock.(targets),priority) 

        all_organisms = get_all_organisms(stock.(sources),stock.(targets)) 

        src_stocks, tgt_stocks, s_enforced,t_enforced,tgt_caps = slot_stocks(sources,targets)
        src_quantities =ustrip.(map(x->uconvert(preferred_quantity_unit(x),JLIMS.quantity(x)),src_stocks))
 
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
        println(typeof(to))
        # check for missing source chemicals needed to complete the destinations
        missing_chemicals=filter(x-> sum(Unitful.ustrip.(source_concentrations[:,Symbol(x.name)]))==0,get_all_chemicals(stock.(targets)))       
        if length(missing_chemicals) > 0 
                throw(MissingChemicalError("No valid source of $(join(name.(missing_chemicals),", ")) available to complete the dispenses",missing_chemicals))
        end 

        # check for missing source organisms needed to create the targets 
        missing_organisms = filter(x-> sum(source_organisms[:,Symbol(JLIMS.name(x))])==0,get_all_organisms(stock.(targets)))
        if length(missing_organisms)>0 
                throw(MissingOrgansimError("No valid source of $(join(JLIMS.name.(missing_organisms),", ")) available to complete the dispenses",missing_organisms))
        end

        # update the chemical priority dict to account for the sources and targets 
        update_priority!(sources,targets,priority)

        W = length(src_stocks) # number of total wells in the model
        C = length(all_chemicals)
        O = length(all_organisms)
        L = length(all_labware)
        R = length(configs)

        for w in 1:W 
                if t_enforced[w]
                        tgt_caps[w] = ustrip(uconvert(u"µL",JLIMS.quantity(tgt_stocks[w])))
                end
        end
        asp,disp,plumb = get_masks(configs,all_labware) 
        lw_idx = vcat([fill(i,length(all_labware[i])) for i in 1:L]...)
        within_labware_index = vcat([1:length(all_labware[i]) for i in 1:L]...)
        Lw_ranges = [findall(x->x==i,lw_idx) for i in 1:L]
        MDict = Dict{Tuple{Int,Int,Int,Int},Int}()
        NDict = Dict{Tuple{Int,Int,Int,Int},Int}()

        m = 1 
        n=1

        Rm = []
        Rn = [] 
        for r in 1:R 
                rmstarting= m
                rnstarting =n
                for l in 1:L 
                        Pa,Cha = asp[r,l][2]
                        Pd,Chd = disp[r,l][2]
                        Mp,pistons = plumb[r,l]
                        for p in 1:Pa
                                for t in 1:pistons
                                        MDict[(r,l,p,t)]=m 
                                        m += 1 
                                end 
                        end 
                        for p in 1:Pd
                                for t in 1:pistons
                                        NDict[(r,l,p,t)]= n 
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
        MVec = collect(keys(MDict))
        NVec = collect(keys(NDict))
        println(Rm)
        println(Rn)
        # set up cost matrix for the cost of each flow 
        costs=ones(M,N)
        for a in eachindex(Rm)
                a_idxs=Rm[a]

                for b in eachindex(Rn)
                        b_idxs=Rn[b]
                        costs[a_idxs,b_idxs] .= sqrt(robot_cost[a]*robot_cost[b])
                end 
        end 



        w_to_m = []
        w_to_n = []
        w_to_m_nozzles=[]
        w_to_m_channelid = []
        w_to_n_channelid = []
        for w in 1:W 
                l = lw_idx[w]
                i = within_labware_index[w]
                wm_mapping = Int[]
                wm_nozzles = Nozzle[]
                wn_mapping = Int[]
                wm_channel = Int[]
                wn_channel = Int[]
                for r in 1:R 
                        
                        Ma = asp[r,l][1]
                        Md = disp[r,l][1]
                        Pa,Cha = asp[r,l][2]
                        Pd,Chd = disp[r,l][2]
                        Mp = plumb[r,l][1]
                        pistons = plumb[r,l][2]
                        for p in 1:Pa 
                                for c in 1:Cha
                                        for t in 1:pistons  
                                                if Ma(i,p,c) && Mp(t,c)
                                                        push!(wm_mapping,MDict[(r,l,p,t)])
                                                        push!(wm_nozzles,channels(head(configs[r]))[c])
                                                        push!(wm_channel,c)
                                                end
                                        end

                                end 
                        end
                        for p in 1:Pd 
                                for c in 1:Chd 
                                        for t in 1:pistons
                                                if Md(i,p,c) && Mp(t,c)
                                                        push!(wn_mapping,NDict[(r,l,p,t)])
                                                        push!(wn_channel,c)
                                                end 
                                        end 
                                end
                        end  
                end 
                push!(w_to_m,wm_mapping)
                push!(w_to_n,wn_mapping)
                push!(w_to_m_nozzles,wm_nozzles)
                push!(w_to_m_channelid,wm_channel)
                push!(w_to_n_channelid,wn_channel)
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
        @variable(model,chem_slacks[1:W,1:C])
        @variable(model,org_slacks[1:W,1:O])
        @variable(model,Qlw[1:L,1:L] >= 0)




        # constraints to ensure that we create the targets properly
        @constraint(model, Q'*sc .- chem_slacks .== tq ) # create the targets with the sources, allowing for some slack. This is a mass/volume balance. 

        if O >0 # only write constraint if there are organisms present  
        @constraint(model, Q'*so .- org_slacks .== inoculation_quantity*to) # ensure that if a strain is dispensed, meaure the discrepancy with the strain slack
        end
        # measure the physical masked operations 
        for i in 1:W 
                for j in 1:W 
                        idxs=Tuple{Int,Int}[]
                        ms = w_to_m[i]
                        ns = w_to_n[j] 
                        cha = w_to_m_channelid[i]
                        chd = w_to_n_channelid[j]
                        for f in eachindex(cha)
                                for g in eachindex(chd) 
                                        if cha[f] ==chd[g] # the pairing must use the same channel to execute the transfer 

                                                push!(idxs,(ms[f],ns[g]))
                                        end
                                end 
                        end 
                        idxs=CartesianIndex.(idxs)
                        if length(idxs) > 0 

                                if i ==j 
                                        @constraint(model,V[idxs] .== 0) # cannot use a 'real' robot to do operations that result in transfering from one well to itself
                                        #=
                                        if s_enforced[i] && t_enforced[j] 
                                                @constraint(model, Q[i,j] == ustrip(uconvert(u"µL",JLIMS.quantity(src_stocks[i]))) ) # fix the transfer amount 
                                        end 
                                        =#
                                else
                                        @constraint(model, sum(V[idxs]) == Q[i,j] ) # map the operations to well to well transfers 
                                end 
                        end 

                        

                end 
        end 

        for i in 1:L 
                for j in 1:L 
                        @constraint(model,Qlw[i,j] == sum(Q[Lw_ranges[i],Lw_ranges[j]]))
                        for r in 1:R 
                                if !can_slot_labware_pair(all_labware[i],all_labware[j],configs[r])
                                        lw_ms = union(w_to_m[Lw_ranges[i]]...)
                                        lw_ns = union(w_to_n[Lw_ranges[j]]...)
                                        ms = intersect(lw_ms,Rm[r])
                                        ns = intersect(lw_ns,Rn[r])

                                        @constraint(model, V[ms,ns] .== 0)
                                end 
                        end

                end 
        end 



        for i in 1:W 
                if s_enforced[i]
                        ms = w_to_m[i]
                        nozzles= w_to_m_nozzles[i]
                        dead::Vector{Real} = map(x->ustrip(uconvert(u"µL",x.deadVol)),nozzles)
                        pad::Vector{Real} = map(x->x.deadVolFactor,nozzles)
                        for n in eachindex(nozzles) 
                                if nozzles[n] isa ContinuousNozzle
                                        pad[n] += dead[n]/ustrip(uconvert(u"µL",nozzles[n].maxAsp)) # the deadvolume should be a factor of the total aspiration for the flow for a continous nozzle
                                
                                else 
                                        pad[n] += dead[n]/src_quantities[i]
                                end
                        end
                        @constraint(model, sum( pad .* V[ms,:]) <= src_quantities[i]) # aspiration constraints 
                        @constraint(model, sum(Q[i,:]) <= src_quantities[i]) # deals with the i to j transfers that aren't tracked
                else
                        @constraint(model, sum(Q[i,:]) == 0 )
                end
                if t_enforced[i]
                        #ns = w_to_n[i]
                        @constraint(model, sum(Q[:,i]) <= tgt_caps[i])
                else
                        @constraint(model, sum(Q[:,i])==0)
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


        if require_nonzero  # requiring nonzero flows for each needed destination well and chemical 
                for c in 1:C
                        sources_with_chemical=sc[:,c] .> 0 
                        for d in 1:W
                                if t_enforced[d]

                                        if tq[d,c] > 0 
                                                @constraint(model,sum( Q[:,d] .* sources_with_chemical) >=  numerical_tolerance ) # check that at least one transfer happens if an ingredient is needed in a destination, even if the optimial solution is to not dispense anything.  
                                        end 
                                end
                        end 
                end 
        end 

        ############################################################################################################################################################
        # SET UP OBJECTIVES

        ############################################################################################################################################################
        if O > 0
        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        org_weights=t_enforced
        T=findall(x->x==true,t_enforced)
        @objective(model, Min , sum(org_slacks[T,:].^2)) 
        optimize!(model) 

        current_slacks=abs.(JuMP.value.(org_slacks[T,:]))
                 
        delta=slack_tolerance * inoculation_quantity*to[T,:] # delta is the tolerance we give to updating the slack in higher priority levels, it is some fraction of the dispense target quantity for every destination. Wiggle room for the slack in future iterations 
        @constraint(model, org_slacks[T,:] .>= -current_slacks .- delta)
        @constraint(model, org_slacks[T,:] .<= current_slacks .+ delta)
        end
        set_objective_sense(model, MOI.FEASIBILITY_SENSE)

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

                end 
                current_slacks=abs.(JuMP.value.(chem_slacks))
                
                for c in 1:C
                        if chem_weights[c]
                                T = findall(x->x==true ,t_enforced)
                                delta=slack_tolerance * tq[T,c] # delta is the tolerance we give to updating the slack in higher priority levels, it is some fraction of the dispense target quantity for every destination. Wiggle room for the slack in future iterations 
                                current_slack=current_slacks[T,c]
                                @constraint(model, chem_slacks[T,c] .>= -current_slack .- delta)
                                @constraint(model, chem_slacks[T,c] .<= current_slack .+ delta)
                        end 
                end 
        
                
        end 
        #reset the optimizer to remove the objective cutoff since we are switching objectives 
        optimize!(model) # resolve one last time with the final slack constraints -> we need to optimize before querying results for the secondary objectives 
        current_slacks=abs.(JuMP.value.(chem_slacks))
        @constraint(model,chem_slacks .== current_slacks)
        #set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        @objective(model, Min,sum(costs .* V)) # minimize total masked operations 
        optimize!(model)
        Vval= JuMP.value.(V) 
        @constraint(model,sum(costs .* V) <= sum( costs .* Vval))

        optimize!(model) # re optimize so we can query, but further solutions will be constrained by the previous objective



        for obj! in secondary_objectives # run through secondary objectives found in objectives.jl 
                obj!(model;numerical_tolerance=numerical_tolerance)
        end 



        

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
        target_chems = get_all_chemicals(stock.(targets))
        c_slacks =JuMP.value.(chem_slacks)
        sol_quality=zeros(size(c_slacks))
        for (d,c)  in Tuple.(CartesianIndices(c_slacks))
                sol_quality[d,c]=0
                if all_chemicals[c] in target_chems 
                        if t_enforced[d]
                                slack= abs(c_slacks[d,c])
                                target = tq[d,c]
                                if target == 0 && isapprox(slack,0;atol=numerical_tolerance)
                                                sol_quality[d,c]= 0 
                                elseif target == 0 && !isapprox(slack,0;atol=numerical_tolerance) 
                                                println("stock target ($(d)), chemical: $(all_chemicals[c]), unbounded error, using max target value for $(all_chemicals[c]) to calculate error")
                                                sol_quality[c] = slack/ maximum(tq[:,c])
                                else
                                                sol_quality[c]= slack/target 
                                end 
                        end
                end
        end 
        max_error,max_idx = findmax(sol_quality)
        if max_error > obj_tolerance
                error("Unacceptable Solution: \n      max error = $(max_error*100)% \n        tolerance limit = $(obj_tolerance*100)% \n      stock_target: $(max_idx[1]) \n  chemical: $(all_chemicals[max_idx[2]])")
                 
        end 
        println("mean solution error: $(mean(sol_quality*100))")



        ops = JuMP.value.(V)

        config_dispenses = [[zeros(length(all_labware[a]),length(all_labware[d])) for a in 1:L , d in 1:L] for r in 1:R]
        slotting_indicators = [falses(L,L) for r in 1:R ]

        for i in 1:W 
                for j in 1:W  
                        la = lw_idx[i]
                        ld = lw_idx[j]
                        a = within_labware_index[i]
                        b = within_labware_index[j]
                        wm = w_to_m[i]
                        wn = w_to_n[j] 
                        wmch= w_to_m_channelid[i]
                        wnch = w_to_n_channelid[j]
                        for r in 1:R 
                                idxs=Tuple{Int,Int}[]
                                wmr_idxs = findall(x-> x in Rm[r],wm)
                                wnr_idxs = findall(x->x in Rn[r],wn)

                                wmr = wm[wmr_idxs]
                                wnr=wn[wnr_idxs]
                                wmrch=wmch[wmr_idxs]
                                wnrch=wnch[wnr_idxs]

                                for f in eachindex(wmrch)
                                        for g in eachindex(wnrch) 
                                                if wmrch[f] == wnrch[g] # the pairing must use the same channel to execute the transfer 
        
                                                        push!(idxs,(wmr[f],wnr[g]))
                                                end
                                        end 
                                end 
                                idxs=CartesianIndex.(idxs)
                                config_dispenses[r][la,ld][a,b] = round(sum(ops[idxs]),digits=round_digits)
                        end 
                end 
        end 
        for r in 1:R 
                slotting_indicators[r]= sum.(config_dispenses[r]) .> 0
        end

        if return_model 
                return config_dispenses,slotting_indicators,model, sol_quality
        else
                return config_dispenses,slotting_indicators
        end

                                

end










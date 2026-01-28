



function slot_stocks(sources::Vector{<:JLIMS.Well},targets::Vector{<:JLIMS.Well})
    all_labware=get_all_labware(sources,targets)
    N= sum(map(x->prod(JLIMS.shape(x)),all_labware))
    src_stocks= JLIMS.Stock[]
    src_enforced=Bool[]
    tgt_stocks=JLIMS.Stock[]
    tgt_enforced=Bool[]
    tgt_caps = Real[]
    src_ids = JLIMS.location_id.(sources)
    tgt_ids=JLIMS.location_id.(targets)
    for l in all_labware 
        for c in children(l)
            id = location_id(c)
            x=findfirst(t->t==id,src_ids)
            y = findfirst(t->t==id,tgt_ids)
            push!(tgt_caps, ustrip(uconvert(u"µL",JLIMS.wellcapacity(c))))

                
            
            if isnothing(x) 
                push!(src_stocks,JLIMS.Empty())
                push!(src_enforced,false)
            else
                push!(src_stocks,JLIMS.stock(sources[x]))
                push!(src_enforced,true)
            end 
            
            if isnothing(y)
                push!(tgt_stocks,JLIMS.Empty())
                push!(tgt_enforced,false)
            else
                push!(tgt_stocks,JLIMS.stock(targets[y]))
                push!(tgt_enforced,true)
            end 
            
        end 
    end 
    return src_stocks,tgt_stocks,src_enforced,tgt_enforced,tgt_caps
end


function can_slot_labware_pair(a::Labware,b::Labware,config::Configuration) 

    open_slots = Set{Tuple{DeckPosition,Int}}()

    for position in deck(config) 
        n_slots = prod(slots(position))
        for i in 1:min(n_slots,2) 
            push!(open_slots,(position,i))
        end 
    end 
    placed_a = false 
    placed_b = false 

    for s in open_slots 
        if can_place(a,s[1])
            placed_a=true 
            delete!(open_slots,s)
            break 
        end
        
    end 

    for s in open_slots 
        if can_place(b,s[1])
            placed_b=true 
            delete!(open_slots,s)
            break 
        end
    end 

    return placed_a && placed_b

end 



function update_priority(sources::Vector{<:JLIMS.Well},targets::Vector{<:JLIMS.Well},priority::PriorityDict)
    new_priority=PriorityDict()
    src_chems= get_all_chemicals(stock.(sources))
    tgt_chems = get_all_chemicals(stock.(targets)) 
        # Update priorities: priority increases with decreasing level

    # Level 0: All level 0 ingredeients are blocked from the design. They may not appear in the stock. 
    # Level 1+: All other ingredients are scheduled sequentially by priority. we allow slack for all nonzero priorities, where we try to minimize the slack for each ingredient and then constrain the slack for future priority levels. 
    # Level 2^(64)-1: Maximum allowable priority level, typically a solvent like water that we would use to back fill will be assigned maximum priority level. 
    for chem in keys(priority) 
        new_priority[chem] = priority[chem]
    end 


    for chem in tgt_chems
        if !in(chem,keys(new_priority))
            new_priority[chem]=UInt64(1)
        end 
    end 
    for chem in src_chems
        if !in(chem,keys(new_priority)) # if the user hasn't specified a source ingredient priority or included it in the targets (see above), assume it is priority 0. (These ingredients will be blocked, regardless of a discrete or continuous robot)
            new_priority[chem]=UInt64(0) 
        end 
    end 
    return new_priority 
end 


function pads(tf::TransferHead)
    pads=Real[]
    for n in channels(tf)
      push!(pads,n.deadVolFactor)
    end 
    return pads 
  end 
  
  function deadvols(tf::TransferHead)
    vols=Real[]
    for n in channels(tf)
        v=ustrip(uconvert(u"µL",n.deadVol))
        push!(vols,v)
    end
    return vols 
end 
    
  



function get_masks(h::TransferHead, l::Labware)
    Ma,Md,Sa,Sd =masks(h,l)
    plumb = plumbing_mask(h) 

    asp = (Ma,Sa)
    disp = (Md,Sd)

    return asp,disp,plumb 
end 


function get_masks(configs::Vector{<:Configuration},all_labware::Vector{<:Labware})
        
    M = map(I -> get_masks(I...),Iterators.product(head.(configs),all_labware))

    asp,disp,plumb = map(x->getindex.(M,x),1:3)
    return asp,disp,plumb
end 





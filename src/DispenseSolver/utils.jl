preferred_quantity_unit(ing::JLIMS.Solid) =u"mg"
preferred_quantity_unit(ing::JLIMS.Liquid) = u"µL"


preferred_quantity_unit(stock::JLIMS.Stock)= u"µL"
preferred_quantity_unit(stock::JLIMS.Mixture) =u"mg"




const chemical_access_dict=Dict(
    JLIMS.Solid => JLIMS.solids ,
    JLIMS.Liquid => JLIMS.liquids
)

PriorityDict = Dict{JLIMS.Chemical,UInt64}

const emptywell= Well{1}(0,"placeholder",nothing,Empty())

""" 
    concentration(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the concentration of an ingredient in a stock using the preferred units for that ingredient and stock 
    

"""
function concentration(stock::JLIMS.Stock,ingredient::JLIMS.Chemical) 
    if ingredient in stock 
        return uconvert(preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient]/JLIMS.quantity(stock))
    else
        return 0*preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock)
    end 
end 


""" 
    quantity(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the quantity of an ingredient in a stock using the preferred units for that ingredient
    

"""
function quantity(stock::JLIMS.Stock,ingredient::JLIMS.Chemical)
    if ingredient in stock
        return uconvert(preferred_quantity_unit(ingredient),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient])
    else
        return 0*preferred_quantity_unit(ingredient)
    end 
end 


function chemical_array(stocks::Vector{<:JLIMS.Stock},ingredients::Vector{<:JLIMS.Chemical};measure::Function=concentration) # can return concentration or quantity 
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

        


function organism_array(stocks::Vector{<:JLIMS.Stock},orgs::Vector{JLIMS.Organism})
    pairs=Iterators.product(orgs,stocks) |> collect 
    out = Base.splat(in).(pairs)

    return DataFrame(out',JLIMS.name.(orgs)) # output expected bo be cultures by strains instead of strains by cultures
end  


function labware_id_and_type(lw::JLIMS.Labware)
    return (JLIMS.location_id(lw),typeof(lw))
end 

function get_all_labware(sources::Vector{<:JLIMS.Well},targets::Vector{<:JLIMS.Well})

    all(map(x-> !isnothing(x),JLIMS.parent.(sources))) || ArgumentError("all source wells must have a searchable parent")
    all(map(x-> !isnothing(x),JLIMS.parent.(targets))) || ArgumentError("all target wells must have a searchable parent")

    all_labware = unique(labware_id_and_type,vcat(JLIMS.parent.(sources),JLIMS.parent.(targets)))

    return all_labware 
end 


function get_all_chemicals(source::JLIMS.Stock)
         # Gather all ingredients contained in the sources, destinations, and priority list 
         solids = chemicals(JLIMS.solids(source))
         liqs = chemicals(JLIMS.liquids(source))
         return collect(union(solids,liqs))
end 

function get_all_chemicals(sources::Vector{<:JLIMS.Stock})
    return collect(union(get_all_chemicals.(sources)...))
end 

function get_all_chemicals(priority::PriorityDict)
    return collect((keys(priority)))
end 

function get_all_chemicals(sources::Vector{<:JLIMS.Stock},targets::Vector{<:JLIMS.Stock},priority::PriorityDict)
    return union(get_all_chemicals(sources),get_all_chemicals(targets),get_all_chemicals(priority))
end 


get_all_organisms(x::JLIMS.Stock) = collect(JLIMS.organisms(x))

function get_all_organisms(sources::Vector{<:JLIMS.Stock})
    return collect(union(get_all_organisms.(sources)...))
end 

function get_all_organisms(sources::Vector{<:JLIMS.Stock},targets::Vector{<:JLIMS.Stock})
    return collect(union(get_all_organisms(sources),get_all_organisms(targets)))
end 



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



function update_priority!(sources::Vector{<:JLIMS.Well},targets::Vector{<:JLIMS.Well},priority::PriorityDict)
    src_chems= get_all_chemicals(stock.(sources))
    tgt_chems = get_all_chemicals(stock.(targets)) 
        # Update priorities: priority increases with decreasing level

    # Level 0: All level 0 ingredeients are blocked from the design. They may not appear in the stock. 
    # Level 1+: All other ingredients are scheduled sequentially by priority. we allow slack for all nonzero priorities, where we try to minimize the slack for each ingredient and then constrain the slack for future priority levels. 
    # Level 2^(64)-1: Maximum allowable priority level, typically a solvent like water that we would use to back fill will be assigned maximum priority level. 

    for chem in tgt_chems
        if !in(chem,keys(priority))
            priority[chem]=UInt64(1)
        end 
    end 
    for chem in src_chems
        if !in(chem,keys(priority)) # if the user hasn't specified a source ingredient priority or included it in the targets (see above), assume it is priority 0. (These ingredients will be blocked, regardless of a discrete or continuous robot)
            priority[chem]=UInt64(0) 
        end 
    end 
    return priority 
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


function get_masks(configs::Vector{<:Configuration},all_labware::Vector{Labware})
        
    M = map(I -> get_masks(I...),Iterators.product(head.(configs),all_labware))

    asp,disp,plumb = map(x->getindex.(M,x),1:3)
    return asp,disp,plumb
end 





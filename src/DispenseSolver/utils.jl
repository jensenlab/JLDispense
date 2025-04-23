preferred_quantity_unit(ing::JLIMS.Solid) =u"mg"
preferred_quantit_unit(ing::JLIMS.Liquid) = u"µL"


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
        return uconvert(preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient]/quantity(stock))
    else
        return 0*preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock)
    end 
end 


""" 
    quantity(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the quantity of an ingredient in a stock using the preferred units for that ingredient
    

"""
function quantity(stock::JLIMS.Stock,ingredient::JLIMS.Chemical)
    if ingredient in ingredients(stock.composition)
        return uconvert(preferred_quantity_unit(ingredient),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient])
    else
        return 0*preferred_quantity_unit(ingredient)
    end 
end 


function chemical_array(stocks::Vector{<:JLIMS.Stock},ingredients::Set{<:JLIMS.Chemical};measure::Function=concentration) # can return concentration or quantity 
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
    return out' # output expected bo be cultures by strains instead of strains by cultures
end  




function get_all_labware(sources::Vector{JLIMS.Well},targets::Vector{JLIMS.Well})

    all(map(x-> !isnothing(x),JLIMS.parent(sources))) || ArgumentError("all source wells must have a searchable parent")
    all(map(x-> !isnothing(x),JLIMS.parent(targets))) || ArgumentError("all target wells must have a searchable parent")

    all_labware = unique(vcat(JLIMS.parent.(sources),JLIMS.parent.(targets)))

    return all_labware 
end 


function get_all_chemicals(source::JLIMS.Stock)
         # Gather all ingredients contained in the sources, destinations, and priority list 
         solids = chemicals(solids(source))
         liqs = chemicals(liquids(source))
         return Set(union(solids,liqs))
end 

function get_all_chemicals(sources::Vector{<:JLIMS.Stock})
    return union(get_all_chemicals.(sources)...)
end 

function get_all_chemicals(priority::PriorityDict)
    return Set(keys(priority))
end 

function get_all_chemicals(sources::Vector{<:JLIMS.Stock},targets::Vector{<:JLIMS.Stock},priority::PriorityDict)
    return union(get_all_chemicals(sources),get_all_chemicals(targets),get_all_chemicals(priority))
end 


get_all_organisms(x::JLIMS.Stock) = JLIMS.organisms(x) 

function get_all_organisms(sources::Vector{<:JLIMS.Stock})
    return union(get_all_organisms.(sources)...)
end 

function get_all_organisms(sources::Vector{<:JLIMS.Stock},targets::Vector{<:JLIMS.Stock})
    return union(get_all_organisms(sources),get_all_organisms(targets))
end 



function slot_stocks(sources::Vector{JLIMS.Well},targets::Vector{JLIMS.Well})
    all_labware=get_all_labware(sources,targets)
    N= sum([prod(shape(l)...) for l in all_labware])
    src_stocks= JLIMS.Stock[]
    src_enforced=Bool[]
    tgt_stocks=JLIMS.Stock[]
    tgt_enforced=Bool[]
    src_ids = JLIMS.location_id.(sources)
    tgt_ids=JLIMS.location_id.(targets)
    for l in all_labware 
        for c in children(l)
            id = location_id(c)
            x=findfirst(y->y==id,src_ids)
            if isnothing(x) 
                push!(src_stocks,JLIMS.Empty())
                push!(src_enforced,false)
            else
                push!(src_stocks,JLIMS.stock(sources[x]))
                push!(src_enforced,true)
            end 
            x = findfirst(y->y==id,tgt_ids)
            if isnothing(x)
                push!(tgt_stocks,JLIMS.Empty())
                push!(tgt_enforced,false)
            else
                push!(tgt_stocks,JLIMS.stock(targets[x]))
                push!(tgt_enforced,true)
            end 
        end 
    end 
    return src_stocks,tgt_stocks,src_enforced,tgt_enforced
end



function update_priority!(sources::Vector{JLIMS.Well},targets::Vector{JLIMS.Well},priority::PriorityDict)
    src_chems= all_chemicals(sources)
    tgt_chems = all_chemicals(targets) 
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


function pad(tf::TransferHead)
    pads::Vector{Real}[]
    for n in nozzles(tf)
      push!(pads,n.deadVolFactor)
    end 
    return pads 
  end 
  
  function deadvols(tf::TransferHead)
    vols::Vector{Real}[]
    for n in nozzles(tf)
        v=ustrip(uconvert(u"µL"),n.deadVol)
        push!(vols,v)
    end
    return vols 
end 
    
  







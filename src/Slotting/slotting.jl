
# slotting functions find a layout for labware on a robot with a particular configuration




function slotting_greedy(labware::Vector{<:Labware},config::Configuration)

    for lw in labware 
        can_place(lw,deck(config)) || error("labware type $(typeof(lw)) cannot be placed on a $(name(config)) deck")
    end 

    slotting = SlottingDict()
    open_slots = Set{Tuple{DeckPosition,Int}}()

    for position in deck(config) 
        n_slots = prod(slots(position))
        for i in 1:n_slots 
            push!(open_slots,(position,i))
        end 
    end 
    not_placed=Labware[]
    for lw in labware 
        placed=false
        for s in open_slots 
            if can_place(lw,s[1])
                slotting[lw]= s 
                delete!(open_slots,s)
                placed=true
                break 
            end 
        end 

        if !placed
            push!(not_placed,lw)
        end 
    end 
    return slotting
end 



function packing_greedy(pairings::Vector{Tuple{Labware,Labware}},config::Configuration;kwargs...)

    slotting_dicts=SlottingDict[]
    remaining= deepcopy(pairings)
    all_labware=union(remaining...)
    while length(remaining) > 0 
        slotting=slotting_greedy(all_labware,config)
        pairings_remaining = filter(x-> any(map(y->!in(y,keys(slotting)),x)),remaining)
        if pairings_remaining == remaining 
            error("unsolvable packing arrangement, $(pairings_remaining)")
        end 
        push!(slotting_dicts,slotting)
        remaining = deepcopy(pairings_remaining)
        if length(remaining)==0 
            break 
        end 
        all_labware = union(remaining...)
    end 
    return slotting_dicts 
end 

function packing_greedy(pairings::Vector{Tuple{Labware,Labware}},config::CobraConfiguration;kwargs...)

    # cobra only has two slots, one dedicated to source and one to destination. A greedy packing arrangement just loops through the pairings
    slotting_dicts=SlottingDict[]
    remaining= deepcopy(pairings)
    for pairing in pairings
        d=SlottingDict(
            pairing[1] => (deck(config)[1],1),
            pairing[2] => (deck(config)[2],1)
        )
        push!(slotting_dicts,d)
    end 
    return slotting_dicts 
end





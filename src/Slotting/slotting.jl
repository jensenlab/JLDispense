
# slotting functions find a layout for labware on a robot with a particular configuration

function slotting_greedy(labware::Vector{<:Labware},config::Configuration)

    slotting = SlottingDict()
    all_slots = Set{Tuple{DeckPosition,Int}}()

    for position in deck(config) 
        n_slots = prod(slots(position))
        for i in 1:n_slots 
            push!(all_slots,(position,i))
        end 
    end 

    for lw in labware 
        for s in all_slots 
            if can_place(lw,s[1])
                slotting[lw]= s 
                delete!(all_slots,s)
                break 
            end 
        end 
        if !in(lw,keys(slotting))
            error("cannot find an open slot for $lw, use differnt labware or change the configuration")
        end 
    end 
    return slotting
end 






"""
    protocol_compilier(dispenses,sources,destinations,config,slotting) 

Given a set of dispense plans `dispenses` from DispenseSolver, find labware slotting configurations for the labware involved, and compile the protocol using the each instruments `dispense` function. 

This function writes the protocols to the current working directory by defualt. The directory can be changed with the keyword argument `directory`

see also: [`JLDispense.dispense`](@ref) 

"""
function protocol_compiler(dispenses::AbstractArray{Matrix{Float64}},sources::Vector{<:Well},destinations::Vector{<:Well},config::Configuration,slotting::BitMatrix;directory=pwd(),packing_method=packing_greedy,kwargs...)



    all_labware=JLDispense.get_all_labware(sources,destinations) 
    all_names = [vec(["$(JLIMS.name(l))_$(JLIMS.name(ch))" for ch in children(l)]) for l in all_labware]

    slotting_pairs = Tuple{Labware,Labware}[]

    for idx in CartesianIndices(slotting)
        i,j=Tuple(idx)
        if slotting[idx]
            push!(slotting_pairs,(all_labware[i],all_labware[j]))
        end 
    end 

    slotting_dicts = packing_method(slotting_pairs,config)
    n_protocols = length(slotting_dicts)
    protocol_names=random_protocol_name(n_protocols)

    for s in eachindex(slotting_dicts)
        labware=collect(keys(slotting_dicts[s]))
        println(labware)
        println(labware[1] in all_labware)
        idxs=map(x->findfirst(y->labware_id_and_type(y) == labware_id_and_type(x),all_labware),labware)
        println(idxs)
        disps=dispenses[idxs,idxs]
        d = hvcat(length(labware),permutedims(disps)...)' # transpose so matrix is dxs instead of sxd. Each source is named as a column in a dataframe
        d=d .* u"µL"
        design = DataFrame(d,vcat(all_names[idxs]...))
        dispense(config,design,directory,protocol_names[s],labware,slotting_dicts[s];kwargs...)
    end 

    

    return nothing 
end 


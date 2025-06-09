

function scheduler(directory::String,sources::Vector{<:Well},targets::Vector{<:Well},configs::Vector{<:Configuration},secondary_objectives...;kwargs...)



    dispenses ,slotting = dispense_solver(sources,targets,configs,secondary_objectives...; return_model=false,kwargs...)

    for r in eachindex(configs)
        instructor(dispenses[r],sources,targets,configs[r],slotting[r];directory=directory,kwargs...)
    end 

end 
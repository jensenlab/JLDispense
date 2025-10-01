
"""
    scheduler(directory,sources,targets,configs,secondary_objectives...;kwargs...)

A wrapper for running `dispense_solver` followed by `protocol_compiler` for each instrument 

see also: [`JLDispense.dispense_solver`](@ref), [`JLDispense.protocol_compiler`](@ref)

"""
function scheduler(directory::String,sources::Vector{<:Well},targets::Vector{<:Well},configs::Vector{<:Configuration},secondary_objectives...;kwargs...)



    dispenses ,slotting = dispense_solver(sources,targets,configs,secondary_objectives...; return_model=false,kwargs...)

    for r in eachindex(configs)
        protocol_compiler(dispenses[r],sources,targets,configs[r],slotting[r];directory=directory,kwargs...)
    end 

end 
"""
    dispense(config::MantisConfiguration, design::DataFrame,directory::AbstractString,protocol_name::AbstractString,destination::Type{<:Labware})

Create Mantis dipsense instructions for SBS plate dispensing

## Arguments 
* `config`: A Configuration object defining the Mantis
* `design`: a (# of sources) x (# of destinations) dataframe containing the volume of each transfer in µL.
* `directory`: output file directory
* `protocol_name`: protocol name  
* `destinations`: A JLIMS.Labware destination subtype  
"""
function dispense(config::Configuration,design::DataFrame, directory::AbstractString,protocol_name::AbstractString,labware::Type{<:Labware};slotting::Function = slotting_greedy, kwargs...)
    lw = map( (x,y)-> x(1,"labware$(i)"),(labware,1:length(labware)))
    layout = slotting(lw,config)
    return dispense(config,design,directory,protocol_name,lw,layout;kwargs...)
end 
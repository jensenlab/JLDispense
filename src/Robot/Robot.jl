abstract type Robot end 


abstract type RobotProperties end 
abstract type RobotConfiguration end 


struct DeckPosition 
    name::AbstractString
    slots::Integer
    compatible_containers::Vector{JLIMS.Container}
end 


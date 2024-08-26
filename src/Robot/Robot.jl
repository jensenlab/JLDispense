abstract type Robot end 


abstract type RobotProperties end 
abstract type RobotConfiguration end 


struct DeckPosition 
    name::AbstractString
    is_source::Bool
    slots::Integer
    compatible_containers::Union{Vector{JLIMS.Container},Missing}
end 


dwp96_2ml=JLIMS.Container("deep_well_plate_96_2_ml",2.0u"mL",(8,12))
dwp96_1ml=JLIMS.Container("deep_well_plate_96_1_ml",1.0u"mL",(8,12))
wp96=JLIMS.Container("plate_96",200.0u"µL",(8,12))
wp384=JLIMS.Container("plate_384",80.0u"µL",(16,24))
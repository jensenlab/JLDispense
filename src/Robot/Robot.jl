abstract type Robot end 


abstract type RobotProperties end 
abstract type RobotConfiguration end 


struct DeckPosition 
    name::AbstractString
    can_aspirate::Bool
    can_dispense::Bool
    slots::Integer
    compatible_containers::Union{Vector{JLIMS.Container},Missing}
end 


dwp96_2ml=JLIMS.Container("deep_well_plate_96_2_ml",2.0u"mL",(8,12))
dwp96_1ml=JLIMS.Container("deep_well_plate_96_1_ml",1.0u"mL",(8,12))
wp96=JLIMS.Container("plate_96",200u"µL",(8,12))
wp384=JLIMS.Container("plate_384",80u"µL",(16,24))
conical_15=JLIMS.Container("conical_15ml",15u"mL",(1,1))
conical_50=JLIMS.Container("conical_50ml",50u"mL",(1,1))




#### import a proxy database of lab objects 
chemicals=JLIMS.parse_chemical_csv("test_ingredients.csv")
strains=JLIMS.parse_strain_csv("test_strains.csv")
water=chemicals[1]
iron_nitrate=chemicals[2]
magnesium_sulfate=chemicals[4]
SMU_UA159=strains[1]
conical_50=Container("conical_50ml",50u"mL",(1,1))
WP96=Container("plate_96",200u"µL",(8,12))
dwp96_2ml=JLIMS.Container("deep_well_plate_96_2_ml",2.0u"mL",(8,12))

function ing(name;ings=chemicals)
    return ings[findfirst(x->x.name==name,ings)]
end 

function str(name; strs=strains)
    return strs[findfirst(x->x.name==name,strs)]
end 
# define stocks 
water=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent"
    )),
    50u"mL",
    JLIMS.Well(1,1,1,conical_50)
)
alanine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("alanine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(2,2,1,conical_50)
)

arginine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("arginine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(3,3,1,conical_50)
)

asparagine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("asparagine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(4,4,1,conical_50)
)

aspartic_acid=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("aspartic_acid")=>5u"g/L",
        ing("hydrochloric_acid")=>0.2u"M"
    )),
    50u"mL",
    JLIMS.Well(5,5,1,conical_50)
)

glutamine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("glutamine")=>20u"g/L"
    )),
    50u"mL",
    JLIMS.Well(6,6,1,conical_50)
)

glutamic_acid=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("glutamic_acid")=>5u"g/L",
        ing("hydrochloric_acid")=>0.2u"M"
    )),
    50u"mL",
    JLIMS.Well(7,7,1,conical_50)
)

glycine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("glycine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(8,8,1,conical_50)
)

histidine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("histidine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(9,9,1,conical_50)
)

isoleucine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("isoleucine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(10,10,1,conical_50)
)

leucine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("leucine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(11,11,1,conical_50)
)

lysine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("lysine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(12,12,1,conical_50)
)

methionine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("methionine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(13,13,1,conical_50)
)

phenylalanine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("phenylalanine")=>10u"g/L",
        ing("hydrochloric_acid")=>0.2u"M"
    )),
    50u"mL",
    JLIMS.Well(14,14,1,conical_50)
)
proline=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("proline")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(15,15,1,conical_50)
)
serine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("serine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(16,16,1,conical_50)
)

threonine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("threonine")=>20u"g/L"
    )),
    50u"mL",
    JLIMS.Well(17,17,1,conical_50)
)
tryptophan=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("tryptophan")=>10u"g/L",
        ing("hydrochloric_acid")=>0.2u"M"
    )),
    50u"mL",
    JLIMS.Well(18,18,1,conical_50)
)

tyrosine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("tyrosine")=>10u"g/L",
        ing("hydrochloric_acid")=>0.2u"M"
    )),
    50u"mL",
    JLIMS.Well(19,19,1,conical_50)
)

valine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("valine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(20,20,1,conical_50)
)

cysteine=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=> 100u"percent",
        ing("l_cysteine")=>10u"g/L"
    )),
    50u"mL",
    JLIMS.Well(21,21,1,conical_50)
)

cdm_2x_glucose1=JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("adenine")=>0.04u"g/L",
        ing("biotin")=> 0.0004u"g/L",
        ing("glucose")=>20u"g/L",
        ing("folic_acid")=>0.0016u"g/L",
        ing("guanine")=>0.04u"g/L",
        ing("iron_sulfate")=>0.01u"g/L",
        ing("iron_nitrate")=>0.002u"g/L",
        ing("magnesium_sulfate")=>1.4u"g/L",
        ing("manganese_sulfate")=>0.01u"g/L",
        ing("b_nadph")=>0.005u"g/L",
        ing("niacinamide")=>0.002u"g/L",
        ing("paba")=>0.0004u"g/L",
        ing("pantothenate")=>0.004u"g/L",
        ing("pyridoxal_hydrochloride")=>0.002u"g/L",
        ing("pyridoxamine")=>0.002u"g/L",
        ing("riboflavin")=>0.004u"g/L",
        ing("sodium_acetate")=>9u"g/L",
        ing("thiamine")=>0.002u"g/L",
        ing("uracil")=>0.04u"g/L",
        ing("vitamin_b12")=>0.0002u"g/L",
        ing("water")=>100u"percent"
    )),
    50u"mL",
    JLIMS.Well(22,22,1,conical_50)
)



smu= Culture(
    str("SMU_UA159"),
    JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=>100u"percent"
    )),
    50u"mL",
    JLIMS.Well(23,23,1,conical_50)
))

ssa= Culture(
    str("SSA_SK36"),
    JLIMS.Stock(
    JLIMS.Composition(Dict(
        ing("water")=>100u"percent"
    )),
    50u"mL",
    JLIMS.Well(24,24,1,conical_50)
))


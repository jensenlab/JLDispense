db_file = "test_jldispense.db"

seed = 52486181

JLIMS.create_db(db_file)

JLIMS.@connect_SQLite db_file 

reservior= generate_location(JLConstants.DeepReservior,"plate_reservior")

deep_well= generate_location(JLConstants.DeepWP96,"deepwell")

wp96=generate_location(JLConstants.WP96,"96-well")
wp384=generate_location(JLConstants.WP384,"384-well")
dwp96=generate_location(JLConstants.DeepWP96,"deep_well")
conical50_1=generate_location(JLConstants.Conical50,"conical_1")
conical50_2=generate_location(JLConstants.Conical50,"conical_2")
conical50_3=generate_location(JLConstants.Conical50,"conical_3")
conical50_4=generate_location(JLConstants.Conical50,"conical_4")
conical50_5=generate_location(JLConstants.Conical50,"conical_5")
conical50_6=generate_location(JLConstants.Conical50,"conical_6")
conical50_7=generate_location(JLConstants.Conical50,"conical_7")
conical50_8=generate_location(JLConstants.Conical50,"conical_8")
conical50_9=generate_location(JLConstants.Conical50,"conical_9")
conical50_10=generate_location(JLConstants.Conical50,"conical_10")
conical50_11=generate_location(JLConstants.Conical50,"conical_11")
conical50_12=generate_location(JLConstants.Conical50,"conical_12")

conicals= [conical50_1,conical50_2,conical50_3,conical50_4,conical50_5,conical50_6,conical50_7,conical50_8,conical50_9,conical50_10,conical50_11,dwp96]


rm(db_file)

chems  = Dict(
    chem"ampicillin" => 0.1u"mg",
    chem"sucrose" => 0.25u"mg",
    chem"chloramphenicol" => 0.3u"mg",
    chem"acid_red_1"=> 1u"µg",
    chem"alanine"=> 50u"µg"
)

orgs = [
    org"SMU_UA159",
    org"SSA_SK36"
]


cdm_2x = 2/5 * cdm_glucose_500mL - 100u"mL" * chem"water"  # 2x cdm with glucose 

for child in children(reservior)
    child.stock = cdm_2x # 2x cdm with glucose 
end 


for child in children(conical50_1)
    child.stock = 50u"mL" * chem"water" 
end 

for child in children(conical50_2)
    child.stock = 50u"mL" * chem"water" 
end 

for child in children(conical50_3)
    child.stock = 0.5u"g" * chem"ampicillin" + 40u"mL"* chem"water" + 10u"mL" * chem"ethanol" 
end

for child in children(conical50_4)
    child.stock = 3u"g" * chem"sucrose" + 10u"mL" * chem"water" 
end 

for child in children(conical50_5)
    child.stock = 0.2u"g" * chem"chloramphenicol" + 20u"mL" * chem"water"
end 

for child in children(conical50_6) 
    child.stock = 2u"g" * chem"chloramphenicol" + 30u"mL" * chem"ethanol" + 20u"mL" * chem"water" 
end 

for child in children(conical50_7)
    child.stock = 0.1u"g" * chem"acid_red_1" + 25u"mL" * chem"water" 
end 

for child in children(conical50_8)
    child.stock = 0.5u"g" * chem"alanine" + 50u"mL" * chem"water" 
end 

for child in children(conical50_9)
    child.stock = org"SMU_UA159" + 25u"mL"* chem"water" 
end 

for child in children(conical50_10)
    child.stock = org"SSA_SK36" + 25u"mL" * chem"water" 
end 

for child in children(conical50_11)
    child.stock= 30u"mL" * chem"water"
end 

conical_11_target=deepcopy(conical50_11)

for child in children(wp96)
    child.stock = 1/1000 * cdm_2x  + 100u"µL" * chem"water"
    n = rand(1:4)
    for _ in 1:n
        chemical = rand(keys(chems))
        val=rand(0:0.1:1)
        child.stock += val*chems[chemical]*chemical
    end
    porg =0.8 
    if rand() < porg 
        org = rand(orgs)
        child.stock += org 
    end 
end


for child in children(dwp96)
    child.stock = 2u"mL" *chem"water" 
end 


    





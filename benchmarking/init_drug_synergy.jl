

db_file = "./benchmarking/benchmark_jldispense_drug_synergy.db"

seed = 52486181

JLIMS.create_db(db_file)

JLIMS.@connect_SQLite db_file 

source_plate=generate_location(JLConstants.DeepWP96,"deepwell")



plt = generate_location(JLConstants.WP96,"96well")
d1reservior=generate_location(JLConstants.DeepReservior,"drug1 reservior")
d2reservior=generate_location(JLConstants.DeepReservior,"drug2 reservior")
media_reservior=generate_location(JLConstants.DeepReservior,"media reservir")

rm(db_file)


s1 = 100u"g" * chem"ampicillin" + 100u"mL" * chem"water"
s2 = 100u"g" * chem"spectinomycin" + 100u"mL" * chem"water" 

m = 1/5 * JLConstants.cdm_glucose_500mL 

children(d1reservior)[1].stock = s1
children(d2reservior)[1].stock = s2
children(media_reservior)[1].stock = m


maxdrug = 4
wellvol = 200
for i in 1:8 
    for j in 1:8 
        children(plt)[i,j].stock = ((maxdrug)*(1-(i-1)/8)u"µL" / quantity(s1) |> NoUnits) * s1  + ((maxdrug)*(1-(j-1)/8)u"µL" / quantity(s2) |> NoUnits )* s2 + ((192u"µL" / quantity(m) |> NoUnits) *m)
    end 
end 





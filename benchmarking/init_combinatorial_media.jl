

db_file = "./benchmarking/comb_media.db"

seed = 52486181

JLIMS.create_db(db_file)

JLIMS.@connect_SQLite db_file 

source_plate=generate_location(JLConstants.DeepWP96,"deepwell")

total_plates = 5

dest_plates = Location[]

for i in 1:total_plates 
    plt = generate_location(JLConstants.WP96,"96well$i")
    push!(dest_plates,plt)
end 
n_reserviors= 3
reserviors= Location[]
for i in 1:n_reserviors
    reservior=generate_location(JLConstants.DeepReservior,"reservior")
    for ch in children(reservior)
        ch.stock = 100u"mL" *chem"water"  
    end
    push!(reserviors,reservior)
end

n_bottles = 3 
bottles = Location[]
for i in 1:n_bottles
    bottle =generate_location(JLConstants.Bottle250mL, "water_bottle$i")
    for ch in children(bottle)
        ch.stock = 100u"mL" * chem"water" 
    end 
    push!(bottles,bottle)
end 


rm(db_file)


include("benchmarking_chems.jl")



chems = collect(keys(chem_dict))[1:48]

for k in eachindex(chems) 
    w = children(source_plate)[k] 
    w.stock = chem_dict[chems[k]] * chems[k] + 2u"mL" * chem"water" 
end 

for k in length(chems)+1:length(chems)+4
    w = children(source_plate)[k]
    w.stock = 2u"mL" * chem"water"
end 














    







db_file = "./benchmarking/benchmark_jldispense.db"

seed = 52486181

JLIMS.create_db(db_file)

JLIMS.@connect_SQLite db_file 

source_plate=generate_location(JLConstants.DeepWP96,"deepwell")

total_plates = 50

dest_plates = Location[]

for i in 1:total_plates 
    plt = generate_location(JLConstants.WP96,"96well$i")
    push!(dest_plates,plt)
end 
n_reserviors= 12 
reserviors= Location[]
for i in 1:n_reserviors
    reservior=generate_location(JLConstants.DeepReservior,"reservior")
    for ch in children(reservior)
        ch.stock = 100u"mL" *chem"water" 
    end
    push!(reserviors,reservior)
end



rm(db_file)


include("benchmarking_chems.jl")



chems = collect(keys(chem_dict))[1:48]

for k in eachindex(chems) 
    w = children(source_plate)[k] 
    x = children(source_plate)[k+48]
    w.stock = chem_dict[chems[k]] * chems[k] + 2u"mL" * chem"water" 
    x.stock= chem_dict[chems[k]] * chems[k] + 2u"mL" * chem"water" 
end 

for n in (length(chems)+1):96
    w = children(source_plate)[n]
    w.stock = 2u"mL" *chem"water" 
end 












    





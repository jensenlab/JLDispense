

using JLIMS,JLConstants,JLDispense, SQLite,CSV,DataFrames, Unitful,StatsBase,Random,Dates,JuMP


Random.seed!(48207531)
date=string(Dates.today())

db_file = "./benchmarking/comb_media.db"
outfile="./benchmarking/benchmark_$date.csv"

JLIMS.create_db(db_file)

JLIMS.@connect_SQLite db_file 


source1 = generate_location(Conical50,"glucose_tube")
source2= generate_location(Conical50,"erm_water")
source3 = generate_location(Conical50,"erm_ethanol")
target = generate_location(Conical50,"target")

w1 = JLIMS.children(source1)[1]
w2 = JLIMS.children(source2)[1]
w3 = JLIMS.children(source3)[1] 
t = JLIMS.children(target)[1]


w1.stock = 16u"mg" * chem"glucose" + 1u"mL" * chem"water" 
w2.stock = 2u"mg" * chem"erythromycin" + 1u"mL" *chem"water" 
w3.stock = 20u"mg" *chem"erythromycin" + 1u"mL" *chem"ethanol" 


t.stock = 10u"mg" * chem"glucose" + 5u"mg" * chem"erythromycin" + 1u"mL" * chem"water" 


sources=[w1,w2,w3]
targets =[t] 

a,b,mdl1 = JLDispense.dispense_solver(sources,targets,[JLDispense.p1000];obj_tolerance=1,return_model=true)

v = JuMP.value.(mdl1[:Q])


vols1 = v[1:3,4] /1000

tprime1 = sum(stock.(sources) .* vols1 )
slacks1= JuMP.value.(mdl1[:chem_slacks])

a,b,mdl2 = JLDispense.dispense_solver(sources,targets,[JLDispense.p1000];obj_tolerance=1,return_model=true,priority=JLDispense.PriorityDict(chem"erythromycin" => UInt64(2), chem"ethanol" => UInt64(0),chem"water"=>typemax(UInt64)))

v = JuMP.value.(mdl2[:Q])


vols2 = v[1:3,4] /1000

tprime2 = sum(stock.(sources) .* vols2 )
slacks2= JuMP.value.(mdl2[:chem_slacks])
a,b,mdl3 = JLDispense.dispense_solver(sources,targets,[JLDispense.p1000];obj_tolerance=1,return_model=true,priority=JLDispense.PriorityDict(chem"ethanol" => UInt64(2),chem"water"=>typemax(UInt64)))

v = JuMP.value.(mdl3[:Q])


vols3 = v[1:3,4] /1000

tprime3 = sum(stock.(sources) .* vols3 )
slacks3= JuMP.value.(mdl3[:chem_slacks])

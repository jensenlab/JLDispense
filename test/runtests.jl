using JLDispense, JLIMS, JLConstants , Unitful, Plots , SQLite,DataFrames
init_file= "init_tests.jl"

include(init_file)





tgt_wells= vec(children(wp96))

src_wells = vcat(map(x->vec(children(x)),conicals)...)

src_wells=vcat(src_wells,children(reservior)[1])

instruments =[JLDispense.p20,JLDispense.p2,JLDispense.platemaster]
ins= [JLDispense.p2]
out = JLDispense.dispense_solver(src_wells,tgt_wells,instruments)

show(out)
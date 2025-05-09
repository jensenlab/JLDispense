using JLDispense, JLIMS, JLConstants , Unitful, Plots , SQLite,DataFrames, StatsBase, JuMP,CSV
init_file= "init_tests.jl"

include(init_file)





tgt_wells= vcat(map(x->vec(children(x)),[wp96,conical_11_target])...)

src_wells = vcat(map(x->vec(children(x)),conicals)...)

src_wells=vcat(src_wells,children(reservior)[1])

instruments =[JLDispense.p20,JLDispense.platemaster,JLDispense.multichannel_p100,JLDispense.mantis_hv,JLDispense.mantis_lv,JLDispense.cobra,JLDispense.nimbus,JLDispense.tempest_lv]
instruments=[JLDispense.nimbus,JLDispense.cobra,JLDispense.platemaster]
inst_cost = [100,1,12,3,5,1,2,3]
ins= [JLDispense.p2]
priority=JLDispense.PriorityDict(
    chem"ethanol" => UInt(2),
    chem"water" => typemax(UInt)
)

sec_objectives = (min_operations!,min_sources!,min_labware_crossover!)
sec_objectives = (min_operations!,) 
sec_objectives=()

disp,slotting,out = JLDispense.dispense_solver(src_wells,tgt_wells,instruments,sec_objectives...;priority=priority,obj_tolerance=1e-2,return_model=true)

timetest = @timed  JLDispense.dispense_solver(src_wells,tgt_wells,instruments,sec_objectives...;priority=priority,obj_tolerance=1e-2,return_model=true)
optimize!(out)
Q=out[:Q] 
V= JuMP.value.(out[:V])

write_files = true 

if write_files 

    cons = all_constraints(out,include_variable_in_set_constraints=true)
    confile= "/Users/BDavid/Desktop/constraints.txt"
    if isfile(confile)
        rm(confile) 
    end 

    open(confile,"w") do file 
        for con in cons 
            println(file,con)
        end 
    end 

    qfile ="/Users/BDavid/Desktop/quantities.csv"
    if isfile(qfile)
        rm(qfile)
    end 
    vfile="/Users/BDavid/Desktop/flows.csv"
    if isfile(vfile)
        rm(vfile)
    end 
    CSV.write(qfile,DataFrame(JuMP.value.(Q),:auto))

    CSV.write(vfile,DataFrame(V,:auto))
end 





println("slotting")
println(sum.(slotting))



println("sum dispenses = $(sum(JuMP.value.(Q))) µL")
println("sum operations = $(sum(V)) µL")
println("time to solve =  $(timetest.time) s" )
disp_directory="./dispense_lists"

if !isdir(disp_directory)
    mkdir(disp_directory)
end

JLDispense.scheduler(disp_directory,src_wells,tgt_wells,instruments,sec_objectives...;priority=priority,obj_tolerance=1e-2)
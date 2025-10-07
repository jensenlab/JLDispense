using JLIMS,JLConstants,JLDispense, SQLite,CSV,DataFrames, Unitful,StatsBase,Random,Dates, JuMP




include("./benchmarking/init_drug_synergy.jl")


robots = [JLDispense.p200,JLDispense.platemaster,JLDispense.multichannel_p100_v,JLDispense.multichannel_p100_h]
sources = [d1reservior,d2reservior,media_reservior]

srcs= vcat(map(x->vec(children(x)),sources)...)
targets= vec(children(plt)[1:8,1:8])

priority=JLDispense.PriorityDict(
     chem"water"=>typemax(UInt64)
)

out = @timed JLDispense.dispense_solver(srcs,targets,robots,min_operations!;priority=priority,numerical_tolerance = 1e-8,return_model=true,timelimit=60)

vals = out.value[1] 
model = out.value[3]
solution_summary(model)
#=
cons = all_constraints(model;include_variable_in_set_constraints=false)

open("cons.txt","w") do f 
    for c in cons 
        println(f,c)
    end 
=#

n = 50 

tm = []

for i in 1:n

    out = @timed JLDispense.dispense_solver(srcs,targets,robots,min_operations!;priority=priority,numerical_tolerance = 1e-8,return_model=true,timelimit=60)
    push!(tm,out.time)
end 
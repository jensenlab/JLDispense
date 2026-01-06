using JLIMS,JLConstants,JLDispense, SQLite,CSV,DataFrames, Unitful,StatsBase,Random,Dates,JuMP

include("./benchmarking/init_combinatorial_media.jl")

Random.seed!(48207531)
date=string(Dates.today())
outfile="./benchmarking/benchmark_$date.csv"

instruments = [JLDispense.p200,JLDispense.multichannel_p1000_h,JLDispense.multichannel_p100_v, JLDispense.nimbus, JLDispense.platemaster,JLDispense.cobra,JLDispense.tempest_hv ]

L=5
p=0.7
for l in 1:L
    plt = dest_plates[l]
    for well in children(plt)
        well.stock = 200u"µL" * chem"water"
        for ch in chems 
            if rand() < p 
                well.stock += 0.1u"mg" * ch
            end 
        end 
    end 
end 

sources = vec(children(source_plate))
rvs= vcat(map(x->vec(children(x)),reserviors[1:1])...)

bts = vcat(map(x->vec(children(x)),bottles[1:1])...)

sources1=vcat(sources,rvs)
sources2=vcat(sources,bts)
targets = vcat(map(x->vec(children(x)),dest_plates)...)

priority=JLDispense.PriorityDict(
    chem"water"=>typemax(UInt64)
)

costs = [2,1,1,1,2,1,1,1]
out1 = @timed JLDispense.dispense_solver(sources1,targets,instruments;robot_cost=costs,priority=priority,obj_tolerance=1e-3,return_model=true,timelimit=60)

out2 = @timed JLDispense.dispense_solver(sources2,targets,instruments;robot_cost=costs, priority=priority,obj_tolerance=1e-3,return_model=true,timelimit=60)

a=map(x-> sum(sum.(x)),out1.value[1])
b= map(x-> sum(sum.(x)),out2.value[1])
inst = ["single channel", "multichannel horizontal","multichannel vertical","nimbus","platemaster","cobra" ,"tempest"]

out_df = DataFrame(instrument=inst, expt1=a,expt2=b)

CSV.write("./combinatorial_media_by_instrument.csv",out_df)

mdl1 = out1.value[3]
mdl2= out2.value[3]
Q1 = JuMP.value.(mdl1[:Q])
V1 = JuMP.value.(mdl1[:V])
Q2 = JuMP.value.(mdl2[:Q])
V2= JuMP.value.(mdl2[:V])

expts = ["experiment1","experiment2"]
sumQ = [sum(Q1),sum(Q2)]
sumV = [sum(V1),sum(V2)]
nQ = [sum(Q1 .>0),sum(Q2 .> 0)]
nV = [sum(V1 .>0),sum(V2 .> 0)]
summary_df = DataFrame(Experiment=expts,TotalVolume=sumQ,TotalFlow=sumV,Transfers=nQ,Flows=nV)
CSV.write("./combinatorial_media_summary.csv",summary_df)
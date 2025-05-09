using JLIMS,JLConstants,JLDispense, SQLite,CSV,DataFrames, Unitful,StatsBase,Random

include("init_benchmark.jl")

Random.seed!(48207531)

robots = [JLDispense.p200,JLDispense.platemaster,JLDispense.multichannel_p100,JLDispense.nimbus,JLDispense.cobra]

KK = [12,24,36,48]
LL = 2:6
RR = 1:5
p = 0.5 
n_reps = 5



expts = length(KK) *length(LL) * length(RR) * n_reps 

KK = 12 
LL = 5 
RR = 1 
n_reps =1

Settings = Iterators.product(KK,LL,RR)

data =DataFrame(chemicals=Int[],wells=Int[],robots=Int[],instruments=String[],error=Real[],time=Real[])
ex = 1 
for n in 1:n_reps 
    for S in Settings 
        k=S[1]
        L=S[2]
        R=S[3]
        chemicals = sample(chems,k;replace=false)
        for l in 1:L
            plt = dest_plates[l]
            for well in children(plt)
                well.stock = 200u"µL" * chem"water"
                for ch in chemicals 
                    if rand() < p 
                        well.stock += 0.1u"mg" * ch
                    end 
                end 
            end 
        end 

        sources = vec(children(source_plate))
        rvs= vcat(map(x->vec(children(x)),reserviors)...)
        sources=vcat(sources,rvs)
        targets = vcat(map(x->vec(children(x)),dest_plates[1:L])...)

        instruments = vcat(robots[1], sample(robots[2:end],R-1,replace=false))

        priority=JLDispense.PriorityDict(
            chem"water"=>typemax(UInt64)
        )

        out = @timed JLDispense.dispense_solver(sources,targets,instruments;priority=priority,obj_tolerance=1e-3,return_model=true,timelimit=300)
        ins = join(map(x->x.name,instruments),",")
        disp,slotting,model,quality =out.value 
        tm = out.time 
        n_wells = L*96+96
        q=maximum(quality)
        push!(data,[k,n_wells,R,ins,q,tm])
        CSV.write("./benchmarking/benchmark.csv",data)
        println("status: $ex / $expts" )
        global ex+=1 

    end 
end 



KK = [12,24,36,48]
LL = 1:5
RR = 1:5
p = 0.5 
n_reps = 5



expts = length(KK) *length(LL) * length(RR) * n_reps 

Settings = Iterators.product(KK,LL,RR)

data =DataFrame(chemicals=Int[],wells=Int[],robots=Int[],instruments=String[],error=Real[],time=Real[])
ex = 1 
for n in 1:n_reps 
    for S in Settings 
        k=S[1]
        L=S[2]
        R=S[3]
        chemicals = sample(chems,k;replace=false)
        for l in 1:L
            plt = dest_plates[l]
            for well in children(plt)
                well.stock = 200u"µL" * chem"water"
                for ch in chemicals 
                    if rand() < p 
                        well.stock += 0.1u"mg" * ch
                    end 
                end 
            end 
        end 

        sources = vec(children(source_plate))
        rvs= vcat(map(x->vec(children(x)),reserviors)...)
        sources=vcat(sources,rvs)
        targets = vcat(map(x->vec(children(x)),dest_plates[1:L])...)

        instruments = vcat(robots[1], sample(robots[2:end],R-1,replace=false))

        priority=JLDispense.PriorityDict(
            chem"water"=>typemax(UInt64)
        )

        out = @timed JLDispense.dispense_solver(sources,targets,instruments;priority=priority,obj_tolerance=1e-3,return_model=true,timelimit=300)

        ins = join(map(x->x.name,instruments),",")
        disp,slotting,model,quality =out.value 
        tm = out.time 
        n_wells = L*96+96
        q=maximum(quality)
        push!(data,[k,n_wells,R,ins,q,tm])
        CSV.write("./benchmarking/benchmark.csv",data)
        println("status: $ex / $expts" )
        global ex+=1 

    end 
end 

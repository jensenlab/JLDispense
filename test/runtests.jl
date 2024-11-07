using JLIMS, JLDispense, Unitful, CSV, DataFrames 

include("init_tests.jl")


srcs=vcat(dwp_stocks,cdm_2x_glucose2,smu)
robots=[cobra_default,tempest_default]

source_compatibility=falses(length(srcs),length(robots))
source_compatibility[1:96,1].= true 
source_compatibility[97:end,2].=true

t,m=@time dispense_solver(srcs,cultures,[cobra_default,tempest_default],minimize_overdrafts!,minimize_robots!,minimize_sources!,minimize_transfers!;quiet=true,return_model=true) 

print(t)


pairs=Iterators.product(srcs,robots) |> collect;

out=Base.splat(is_compatible_source).(pairs)

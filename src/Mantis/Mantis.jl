
struct MantisProperties <: RobotProperties 
    isdiscrete::Bool
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    positions::Vector{DeckPosition}
    compatible_stocks::Vector{DataType}
end 
  
  mutable struct MantisConfiguration <: RobotConfiguration
    type::AbstractString #tube or pipette tip
    enable_pause::Bool
    prime_vol::Unitful.Volume
    predispense_vol::Unitful.Volume
    revovery_vol::Unitful.Volume
  end
  
  struct Mantis <:Robot
    name::AbstractString 
    properties::MantisProperties 
    configuration::MantisConfiguration
  
  end
  
  const mantis_names=Dict{JLIMS.Container,String}(
  wp96=>"PT3-96-Assay.pd.txt",
  wp384=>"PT9-384-Assay.pd.txt")

  
  const default_mantis_hv=DeckPosition("High Volume",true,1,missing)
  const default_mantis_lc3=DeckPosition("LC3",true,1,missing)
  const default_mantis_destination=DeckPosition("Destination",false,1,[wp96,wp384])
  
  mantis_default = Mantis("Default Mantis",
  MantisProperties(true,0.1u"µL",100u"mL",vcat([default_mantis_hv for _ in 1:4],[default_mantis_lc3 for _ in 1:24],[default_mantis_destination]),[JLIMS.LiquidStock]),
  MantisConfiguration("pipette_tip",false,6u"µL",0.5u"µL",6u"µL")
  )
#vcat([default_mantis_hv for _ in 1:4],[default_mantis_lc3 for _ in 1:24]),


"""
    mantis(directory::AbstractString,protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Mantis) where {T <: JLIMS.Stock,U <:JLIMS.Stock}

Create Mantis dipsense instructions for SBS plate dispensing

  ## Arguments 
  * `directory`: output file directory
  * `protocol_name`: protocol name  
  * `design`: a (# of sources) x (# of destinations) dataframe containing the volume of each transfer in µL.
  * `sources`: the source stocks
  * `destinations`: the destination stocks 
"""
function mantis(directory::AbstractString,protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Mantis) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
        s_slots=sum(map(x->x.slots,filter(x->x.is_source,robot.properties.positions)))
        s_labware=unique(map(x->x.well.labwareid,sources))
        s_slots >= length(s_labware) ? nothing : error("Number of source labware ($(length(s_labware))) exceeds capacity for $(robot.name) ($s_slots)") 

        allequal(map(x->x.well.labwareid,destinations)) ? nothing : error("All destination stocks must be on the same labware")

        dest=unique(map(x->x.well.container,destinations))[1]
        R,C=dest.shape
        source_names=map(x->"WellID $(x.well.id)",sources)
        dispenses=DataFrame(zeros(R*C,length(sources)),source_names)

        for d in eachindex(destinations) 
          idx=destinations[d].well.wellindex
          dispenses[idx,:] .= ustrip.(uconvert.((u"µL",),design[:,d]))
        end 
        full_dir=joinpath(directory,protocol_name)
        if ~isdir(full_dir)
            mkdir(full_dir)
        end 
        dispenses=dispenses[:,findall(x->sum(dispenses[:,x]) > 0,names(dispenses))]
        filename=joinpath(full_dir,protocol_name*".dl.txt")
        n_stocks=ncol(dispenses)
        delay_header=vcat(n_stocks,repeat([0,""],n_stocks))
        outfile=open(filename,"w")
        platefilename=mantis_names[dest]
        print(outfile,join(["[ Version: 5 ]"],'\t'),"\r\n")
        print(outfile,platefilename,"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")
        print(outfile,join([1],'\t'),"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")

        for i in 1:n_stocks 
            vols=Vector(dispenses[:,i])
            vols=round.(reshape(vols,R,C),digits=1)
            print(outfile,join([names(dispenses)[i],"","Normal"],'\t'),"\r\n")
            print(outfile,join(["Well",1],'\t'),"\r\n")
            for r in 1:R
                print(outfile,join(vols[r,:],'\t'),"\r\n")
            end 
        end 
        close(outfile)
        write(joinpath(full_dir,"config.json"),JSON.json(robot))
        return make_transfer_table(sources,destinations,design)
end 


function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Mantis,kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}

  design=dispense_solver(sources,destinations,robot,minimize_overdrafts!,minimize_sources!,minimize_transfers!;pad=1.1,kwargs...)

  srcs_needed=filter(x->ustrip(sum(design[x,:])) > 0,eachindex(sources))

  srcs=sources[srcs_needed]

  des=design[srcs_needed,:]
  all_labware=unique(map(x->x.well.labwareid,destinations))
  dest_idxs=[findall(x->x.well.labwareid==i,destinations) for i in all_labware]
  tt=DataFrame[]
  pn=AbstractString[]
  for i in eachindex(all_labware)
    protocol_name=random_protocol_name()
    push!(tt,mantis(directory,protocol_name,des[:,dest_idxs[i]],srcs,destinations[dest_idxs[i]],robot))
    push!(pn,protocol_name)
  end 
  return pn,tt
end 

#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]

protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,mantis_default)

=# 
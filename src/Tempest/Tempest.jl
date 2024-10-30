struct TempestProperties <: RobotProperties 
    isdiscrete::Bool
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    positions::Vector{DeckPosition}
    compatible_stocks::Vector{DataType}
end 
  
  mutable struct TempestConfiguration <: RobotConfiguration
    type::AbstractString #tube or pipette tip
    enable_pause::Bool
    prime_vol::Unitful.Volume
    predispense_vol::Unitful.Volume
    revovery_vol::Unitful.Volume
  end
  
  struct Tempest <:Robot
    name::AbstractString 
    properties::TempestProperties 
    configuration::TempestConfiguration
  
  end
  
  const tempest_names=Dict{JLIMS.Container,String}(
  wp96=>"PT3-96-Assay.pd.txt",
  wp384=>"PT9-384-Assay.pd.txt")

  
  const default_tempest_source=DeckPosition("Front Rack",true,false,6,missing)
  const default_tempest_destination=DeckPosition("Destination",false,true,1,[wp96,wp384])
  
  tempest_default = Tempest("Default Tempest",
  TempestProperties(true,0.2u"µL",1000u"mL",vcat([default_tempest_source],[default_tempest_destination for _ in 1:10]),[JLIMS.LiquidStock]),
  TempestConfiguration("tube",false,300u"µL",0.5u"µL",300u"µL")
  )


#= Adam Dama Python Code 
        with open(path, "w", newline="", encoding="utf-8") as f:
            # Header
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["Version            :", 6])
            writer.writerow(
                ["Plate type name    :", plate_type.formulatrix_plate_template_name]
            )
            writer.writerow(["Priority Delays    :", 0])

            # Pipette volumes
            for idx, r in enumerate(reagents):
                writer.writerow(["Reagent Name    :", r])
                writer.writerow(["Barcode            :"])
                writer.writerow(["Priority           :", 1])
                vols = worklist.iloc[:, idx].values
                vols = vols.reshape(plate_type.shape)
                writer.writerows(vols)
=# 
"""
    tempest(design::DataFrame,filepath::String,destination::ContainerName)

Create Tempest dipsense instructions for SBS plate dispensing

  ## Arguments 
  * `design`: a (# of runs) x (# of stocks) dataframe containing the volume of each stock for each run in µL.
  * `filepath`: the ouput path of the dispense file. The function automatically strips any incorrect file extension and adds the appropriate one.
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
"""
function tempest(directory::AbstractString,protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Tempest) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    s_slots=sum(map(x->x.slots,filter(x->x.is_source,robot.properties.positions)))
    s_labware=unique(map(x->x.well.labwareid,sources))
    s_slots >= length(s_labware) ? nothing : error("Number of source labware ($(length(s_labware))) exceeds capacity for $(robot.name) ($s_slots)") 

    allequal(map(x->x.well.labwareid,destinations)) ? nothing : error("All destination stocks must be on the same labware")
    destination=destinations[1].well.container
    R,C=destination.shape
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
    filename=joinpath(full_dir,protocol_name*".dl.txt")
    n_stocks=ncol(dispenses)
    #delay_header=vcat(n_stocks,repeat([0,""],n_stocks))
    outfile=open(filename,"w")
    platefilename=tempest_names[destination]
    if platefilename==""
        error(ArgumentError("Selected plate type $(destination.name) is not supported by Tempest."))
    end 
    print(outfile,join(["Version            :", 6],'\t'),"\r\n")
    print(outfile,join(["Plate type name    :", platefilename],'\t'),"\r\n")
    print(outfile,join(["Priority Delays    :", 0],'\t'),"\r\n")


    for i in 1:n_stocks 
        vols=Vector(dispenses[:,i])
        vols=reshape(vols,R,C)
        print(outfile,join(["Reagent Name    :", source_names[i]],'\t'),"\r\n")
        print(outfile,join(["Barcode            :"],'\t'),"\r\n")
        print(outfile,join(["Priority           :", 1],'\t'),"\r\n")
        for r in 1:R
            print(outfile,join(vols[r,:],'\t'),"\r\n")
        end 
    end 
    close(outfile)
    write(joinpath(full_dir,"config.json"),JSON.json(robot))
    return make_transfer_table(sources,destinations,design)
end 

"""
    multi_tempest(designs::Vector{DataFrame},dlnames::Vector{String},directory::String,mdlname::String,destination::ContainerName)

Create Tempest multidipsense instructions for SBS plate dispensing

  ## Arguments 
  * `design`: a vector of k dataframes where each dataframe is a (# of runs) x (# of stocks) dataframe containing the volume of each stock for each run in µL.
  * `dlnames`" The name of each dispense list. There should be k names, one for each design.
  * `directory`: the ouput directory of the dispense file. The function automatically strips any incorrect file extension and adds the appropriate one.
  * `mdlname`: The name of the multidispense list file. The funcion outputs the multidispense list in directory
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
"""
function multi_tempest(directory::AbstractString,protocol_name::AbstractString,design::DataFrame,sources::Vector{T},destinations::Vector{U},robot::Tempest) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    s_slots=sum(map(x->x.slots,filter(x->x.is_source,robot.properties.positions)))
    s_labware=unique(map(x->x.well.labwareid,sources))
    s_slots >= length(s_labware) ? nothing : error("Number of source labware ($(length(s_labware))) exceeds capacity for $(robot.name) ($s_slots)") 

    d_slots=sum(map(x->x.slots,filter(x->!x.is_source,robot.properties.positions)))
    d_labware=unique(map(x->x.well.labwareid,destinations))
    d_slots >= length(d_labware) ? nothing : error("Number of destination labware ($(length(d_labware))) exceeds capacity for $(robot.name) ($d_slots)") 

    d_idxs=[findall(x->x.well.labwareid==i,destinations) for i in d_labware]

    dlnames=[protocol_name*"_$i" for i in eachindex(d_idxs)]

    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
        mkdir(full_dir)
    end 




    for i in eachindex(d_idxs)
        dests=destinations[d_idxs[i]]
        d = design[:,d_idxs[i]]
        if sum(Matrix(d)) > 0u"µL"
            tempest(full_dir,dlnames[i],d,sources,dests,robot)
        end 
    end 

    mdlbasename=protocol_name*".mdl.txt"
    outpath=joinpath(full_dir,mdlbasename)

    

    outfile=open(outpath,"w")

    print(outfile,join(["LP:$(destination.tempest).pd.txt"],'\t'),"\r\n")
    for i in 1:n
        print(outfile,join(["P:$i","TP:1","DL:$(basenames[i])"],'\t'),"\r\n")
    end
    close(outfile)
    write(joinpath(full_dir,"config.json"),JSON.json(robot))
    return make_transfer_table(sources,destinations,design)
end 



function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Tempest,kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    design=dispense_solver(sources,destinations,robot,minimize_overdrafts!,minimize_sources!,minimize_transfers!;pad=1.1,kwargs...)


    srcs_needed=filter(x->ustrip(sum(design[x,:])) > 0,eachindex(sources))

    srcs=sources[srcs_needed]

    des=design[srcs_needed,:]

    s_slots=sum(map(x->x.slots,filter(x->x.is_source,robot.properties.positions)))
    s_labware=unique(map(x->x.well.labwareid,srcs))
    S=length(s_labware)

    d_slots=sum(map(x->x.slots,filter(x->!x.is_source,robot.properties.positions)))
    d_labware=unique(map(x->x.well.labwareid,destinations))
    D=length(d_labware)
    outer_runs=cld(S,s_slots)
    inner_runs=cld(D,d_slots)

    sl_set=[(s_slots*(i-1) + 1):min(s_slots*(i),S) for i in 1:outer_runs]
    dl_set=[(d_slots*(i-1) + 1):min(d_slots*(i),D) for i in 1:inner_runs]

    s_idxs=[findall(x->in(x.well.labwareid,s_labware[sl_set[i]]),srcs) for i in 1:outer_runs]
    d_idxs=[findall(x->in(x.well.labwareid,d_labware[dl_set[i]]),destinations) for i in 1:inner_runs]
    tt=DataFrame[]
    pn=AbstractString[]
    for i in 1:outer_runs 
        for j in 1:inner_runs 
            protocol_name=random_protocol_name()
            if length(dl_set[j])>1 
                push!(tt,multi_tempest(directory,protocol_name,des[s_idxs[i],d_idxs[j]],srcs[s_idxs[i]],destinations[d_idxs[j]],robot))
            else
                push!(tt,tempest(directory,protocol_name,des[s_idxs[i],d_idxs[j]],srcs[s_idxs[i]],destinations[d_idxs[j]],robot))
            end 
            push!(pn,protocol_name)
        end 
    end 
    return pn,tt

end 

#= 

using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,tempest_default)
=# 
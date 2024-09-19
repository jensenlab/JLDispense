
struct NimbusProperties <: RobotProperties 
    isdiscrete::Bool
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxASP::Unitful.Volume
    positions::Vector{DeckPosition}
    compatible_stocks::Vector{DataType}
  end 
  
  mutable struct NimbusConfiguration <: RobotConfiguration
    max_tip_use::Integer
    comment::AbstractString
  end
  
  struct Nimbus <:Robot
    name::AbstractString 
    properties::NimbusProperties 
    configuration::NimbusConfiguration
  end









# Define our available racks
#15 mL tube
TubeRack15ML_0001=DeckPosition("TubeRack15ML_0001",true,24,[conical_15])
# 50 mL tube
TubeRack50ML_0001=DeckPosition("TubeRack50ML_0001",true,6,[conical_50])
TubeRack50ML_0002=DeckPosition("TubeRack50ML_0002",true,6,[conical_50])
TubeRack50ML_0003=DeckPosition("TubeRack50ML_0003",true,6,[conical_50])
TubeRack50ML_0004=DeckPosition("TubeRack50ML_0004",true,6,[conical_50])
TubeRack50ML_0005=DeckPosition("TubeRack50ML_0005",true,6,[conical_50])
TubeRack50ML_0006=DeckPosition("TubeRack50ML_0006",true,6,[conical_50])
# 2 mL deep well plate
Cos_96_DW_2mL_0001=DeckPosition("Cos_96_DW_2mL_0001",false,1,[dwp96_2ml,wp96])
Cos_96_DW_2mL_0002=DeckPosition("Cos_96_DW_2mL_0002",false,1,[dwp96_2ml,wp96])


default_nimbus_config=[
    TubeRack50ML_0001,
    TubeRack50ML_0002,
     TubeRack50ML_0003,
     TubeRack50ML_0004,
     TubeRack50ML_0005,
     TubeRack50ML_0006,
    Cos_96_DW_2mL_0001
]

nimbus_default=Nimbus("Nimbus Default",
NimbusProperties(false,50u"µL",1u"mL",1000u"µL",default_nimbus_config,[JLIMS.LiquidStock]),
NimbusConfiguration(25,"default")
)

function rack_codes(n) 
    codes=["" for _ in 1:n]
    alphabet=collect('A':'Z')
    k=length(alphabet)
    i=1
    alphacounter=0
    for i = 1:n
        codes[i]=repeat(alphabet[mod(i-1,k)+1],cld(i,k))
    end 
    return codes 
end 


function circle(x,y,r)

    θ=LinRange(0,2*π,500)
    x .+ r*sin.(θ), y .+ r*cos.(θ)
end 

function visualize_rack(name::AbstractString,setup::Vector{AbstractString})
    n=length(setup)
    size=(1,1)
    fontsize=14
    wrapwidth=20
    if n == 6 
        size=(2,3)
        wrapwidth=30

    elseif n==24 
        size=(4,6)
        fontsize=10
        wrapwidth=20
    else
        ArgumentError("Only 50mL and 15mL racks are supported")
    end 
    r,c=size
    vals=reshape(setup,r,c)
    plt=plot(grid=false,size=(1200,800),yflip=true,legend=false,dpi=300,xticks=1:c,yticks=(1:r,rack_codes(r)),xmirror=true,tickdirection=:none,tickfontsize=18)

    for x in 1:c
        for y in 1:r 
            annotate!(x,y,text(TextWrap.wrap(vals[y,x],width=wrapwidth),:center,fontsize))
            plot!(circle(x,y,0.5),color="black")
        end 
    end 


    plot!(ylims=(0.5,r+0.5),xlims=(0.5,c+0.5))
    plot!(title=name,titlefontsize=20)
    return plt 




end 

function convert_nimbus_design(design::DataFrame,sources::Vector{T},destinations::Vector{U},source_deck::DeckPosition,destination_deck::DeckPosition,robot::Nimbus) where {T <: JLIMS.Stock,U <:JLIMS.Stock}


    if !in(source_deck,robot.properties.positions) || !in(destination_deck,robot.properties.positions)
        ArgumentError("Provide a valid Source Rack and Destination Rack Code")
    end 

    if ncol(design) != length(destinations)
        ArgumentError("Number of columns in design must match the size of the destination")
    end 


    if nrow(design) != length(sources)
        ArgumentError("Number of rows in design must match the size of the source")
    end 

    source_id=String[]
    source_position=Int64[]
    volume=Unitful.Volume[]
    destination_id=String[]
    destination_position=String[]
    alphabet=collect('A':'Z')
    R,C=destinations[1].well.container.shape
    for row in 1:nrow(design)
        for col in 1:ncol(design)
            push!(source_id,source_deck.name)
            push!(source_position,row)
            push!(volume,design[row,col])
            push!(destination_id,destination_deck.name)

            c=cld(col,R)
            r=col-R*(c-1)
            pos=string(alphabet[r],c)
            push!(destination_position,pos)
        end 
    end 

    out=DataFrame("Source Labware ID"=>source_id,
        "Source Position ID"=> source_position,
        "Volume (uL)"=> volume,
        "Destination Labware ID"=> destination_id,
        "Destination Position ID"=> destination_position,
    )

    return out
end

"""
    nimbus(dispense_list::NimbusDispenseList,filepath::String)

Create Hamilton Nimbus dipsense instructions

  ## Arguments 
  * `dispense_list`: a NimbusDispenseList object encoding the design source and destination.
  * `filepath`: the ouput path of the dispense file in a .csv format

  ## Keyword Arguments
  * `nimbus_config`: A vector of NimbusRack objects that specify the configuration of the numbus. the defual configuration is 5 50ml conical racks and a well plate rack. 
"""
function nimbus(directory::AbstractString, protocol_name::AbstractString, design::DataFrame, sources::Vector{T}, destinations::Vector{U}, robot::Nimbus;kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    source_decks=filter(x->x.is_source,robot.properties.positions)
    s_slots=map(x->x.slots,source_decks)
    s_labware=unique(map(x->x.well.labwareid,sources))

    S=length(s_labware)
    sum(s_slots) >= S ? nothing : error("source labware ($S) exceed the number of available Nimbus source slots $(s_slots)")

    destination_decks=filter(x->!x.is_source,robot.properties.positions)
    d_slots=map(x->x.slots,destination_decks)
    d_labware=unique(map(x->x.well.labwareid,destinations))
    D=length(d_labware)
    sum(d_slots) >= D ? nothing : error("destination labware ($D) exceed the number of available Nimbus destination slots $(d_slots)")
    

    s_idxs=[findall(x->x.well.labwareid==i,sources) for i in s_labware]
    d_idxs=[findall(x->x.well.labwareid==i,destinations) for i in d_labware]
    cum_s_slots=vcat(0,cumsum(s_slots))
    cum_d_slots=vcat(0,cumsum(d_slots))
    s_runs=findfirst(x->x>=S,cumsum(s_slots))
    d_runs=findfirst(x->x>=D,cumsum(d_slots))
    sl_set=[(cum_s_slots[i]+ 1):min(cum_s_slots[i+1],S) for i in 1:s_runs]
    dl_set=[(cum_d_slots[i] + 1):min(cum_d_slots[i+1],D) for i in 1:d_runs]

    source_loading=DataFrame("Rack"=>AbstractString[],"Position"=>Integer[],"LabwareID"=>AbstractString[])
    destination_loading=DataFrame("Rack"=>AbstractString[],"Position"=>Integer[],"LabwareID"=>AbstractString[])

    dispenses=DataFrame("Source Labware ID"=>AbstractString[],"Source Position ID"=>Integer[],"Volume (uL)"=>Real[],"Destination Labware ID"=>AbstractString[],"Destination Position ID"=>AbstractString[])
    for d in 1:d_runs
        for s in 1:s_runs 
            ss=vcat(s_idxs[sl_set[s]]...)
            dd=vcat(d_idxs[dl_set[d]]...)
            df=convert_nimbus_design(design[ss,dd],sources[ss],destinations[dd],source_decks[s],destination_decks[d],robot)
            src_load=DataFrame("Rack"=>[source_decks[s].name for _ in 1:s_slots[s]],"Position"=> collect(1:s_slots[s]),"LabwareID"=>vcat(map(x->x.well.labwareid,sources[ss]),["" for _ in 1:(s_slots[s]-length(sources[ss]))]))
            dest_load=DataFrame("Rack"=>[destination_decks[d].name for _ in 1:d_slots[d]], "Position"=>collect(1:d_slots[d]),"LabwareID"=>vcat(unique(map(x->x.well.labwareid,destinations[dd])),["" for _ in 1:(d_slots[d]-length(unique(map(x->x.well.labwareid,destinations[dd]))))]))
            append!(source_loading,src_load)
            append!(destination_loading,dest_load)
            n=nrow(df)
            dispense_df=DataFrame([[],[],[],[],[],],names(df))
            for i = 1:n 
                vol=df[i,"Volume (uL)"]
                while vol > 1e-6*u"µL" 
                    shotvol=min(robot.properties.maxVol,vol) # maximum shot volume of 1 ml 
                    push!(dispense_df,(df[i,"Source Labware ID"],df[i,"Source Position ID"],ustrip(uconvert(u"µL",shotvol)),df[i,"Destination Labware ID"],df[i,"Destination Position ID"]))
                    vol-=shotvol
                end 
            end 
            append!(dispenses,dispense_df)
        end 
    end 

    n=nrow(dispenses)
    
    change_tip=zeros(Int64,n)
    for i in 2:n
        if dispenses[i,"Source Position ID"] != dispenses[i-1,"Source Position ID"]
            change_tip[i]=1
        end 
    end 

    windowsize=robot.configuration.max_tip_use
    for i in 1:n-windowsize+1   
        window=change_tip[i:i+windowsize-1]
        if sum(window)==0
            change_tip[i+windowsize-1]=1
        end 
    end 
    dispenses[!,"Change Tip Before"].= change_tip
    
    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
        mkdir(full_dir)
    end 

    source_loading=unique(source_loading)
    destination_loading=unique(destination_loading)

    src_racks=unique(source_loading[:,"Rack"])
    for rack in src_racks 
        entries=subset(source_loading, :Rack => x->x.==rack)
        sort!(entries,[:Position])
        plt=visualize_rack(rack,entries[:,:LabwareID])
        png(plt,joinpath(full_dir,rack*".png"))
    end 
    CSV.write(joinpath(full_dir,"destination_loading.csv"),destination_loading)


    
    CSV.write(joinpath(full_dir,protocol_name*".csv"),dispenses)
    write(joinpath(full_dir,"config.json"),JSON.json(robot))
    return make_transfer_table(sources,destinations,design)
end 


function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Nimbus,kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
    design=dispense_solver(sources,destinations,robot,minimize_overdrafts!,minimize_sources!,minimize_transfers!;pad=1.25,kwargs...)


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
            push!(tt,nimbus(directory,protocol_name,des[s_idxs[i],d_idxs[j]],srcs[s_idxs[i]],destinations[d_idxs[j]],robot))
            push!(pn,protocol_name)
        end 
    end 
    return pn,tt

end 

#= 
using JLDispense,JLD2,JLIMS

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
all_actives=JLD2.load("./src/Mixer/all_actives.jld2")["all_actives"]
#t,m=dispense_solver(all_actives,sources,nimbus_default;return_model=true)
protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",all_actives,sources,nimbus_default)

=#


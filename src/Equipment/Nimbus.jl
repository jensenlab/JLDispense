

#=
struct Nozzle
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    is_discrete::Bool
    multidispense::Bool
end=#

const nimbus_nozzle= ContinuousNozzle(50u"µL",1000u"µL",1000u"µL",25u"µL",1,false)


abstract type NimbusHead <: TransferHead end 

struct NimbusSingleChannelHead <: NimbusHead 
    nozzles::Nozzle
    NimbusSingleChannelHead()=new(nimbus_nozzle)
end 



# Define our available racks
#15 mL tube
tuberack15mL_0001=SBSPosition("TubeRack15ML_0001",Set([JLConstants.Conical15]),(4,6),true,true,false,false,circle)
# 50 mL tube
tuberack50mL_0001=SBSPosition("TubeRack50ML_0001",Set([JLConstants.Conical50]),(2,3),true,true,false,false,circle)
tuberack50mL_0002=SBSPosition("TubeRack50ML_0002",Set([JLConstants.Conical50]),(2,3),true,true,false,false,circle)
tuberack50mL_0003=SBSPosition("TubeRack50ML_0003",Set([JLConstants.Conical50]),(2,3),true,true,false,false,circle)
tuberack50mL_0004=SBSPosition("TubeRack50ML_0004",Set([JLConstants.Conical50]),(2,3),true,true,false,false,circle)
tuberack50mL_0005=SBSPosition("TubeRack50ML_0005",Set([JLConstants.Conical50]),(2,3),true,true,false,false,circle)
tuberack50mL_0006=SBSPosition("TubeRack50ML_0006",Set([JLConstants.Conical50]),(2,3),true,true,false,false,circle)

# 2 mL deep well plate
Cos_96_DW_2mL_0001=SBSPosition("Cos_96_DW_2mL_0001",Set([JLConstants.DeepWP96,JLConstants.WP96]),(1,1),true,true,false,false,rectangle)
Cos_96_DW_2mL_0002=SBSPosition("Cos_96_DW_2mL_0002",Set([JLConstants.DeepWP96,JLConstants.WP96]),(1,1),true,true,false,false,rectangle)


nimbus_deck = [Cos_96_DW_2mL_0001 tuberack50mL_0001 tuberack50mL_0002 EmptyPosition(); tuberack50mL_0003 tuberack50mL_0004 tuberack50mL_0005 tuberack50mL_0006]

struct NimbusSettings <: InstrumentSettings 
    max_tip_use::Integer
end 

NimbusConfiguration = Configuration{NimbusHead,Deck,NimbusSettings}

const nimbus = NimbusConfiguration("Nimbus",NimbusSingleChannelHead(),nimbus_deck,NimbusSettings(25))




function masks(h::NimbusSingleChannelHead,l::Labware)
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    Md=Ma 
    return Ma,Md
end 


function masks(h::NimbusSingleChannelHead,l::JLConstants.WP384)
    function Ma(w::Integer,p::Integer,c::Integer) # Nimbus cannot access 384 well plates
        return false
    end 
    Md=Ma
    return Ma,Md 
end 


function convert_design(design::DataFrame,labware::Vector{<:Labware}, slotting::SlottingDict,config::NimbusConfiguration) 
    # helper function that converts the design and slotting scheme into operations for the nimbus 
    
    all(map(x-> x[1] in deck(config),values(slotting))) || ArgumentError("All deck positions in the SlottingDict must be present on the deck")


    source_id=String[]
    source_position=Union{String,Integer}[]
    volume=Unitful.Volume[]
    destination_id=String[]
    destination_position=Union{String,Integer}[]
    alphabet=collect('A':'Z')
    for col in 1:ncol(design) #sources
        source = get_labware(labware,col)
        s_id , s_pos = slotting[source]
        for row in 1:nrow(design) #destinations
            if design[row,col] == 0 
                continue # skip if no volume transferred 
            end 
            push!(source_id,s_id.name)
            if length(source) ==1 
                push!(source_position,s_pos)
            else
                r,c = cartesian(falses(shape(source)...),col)
                pos=string(alphabet[r],c)
                push!(source_position,pos)
            end 
            push!(volume,design[row,col])
            destination =get_labware(labware,row)
            d_id,d_pos = slotting[destination]
            push!(destination_id,d_id)
            if length(destination) == 1 
                push!(destination_position,d_pos)
            else
                r,c = cartesian(falses(shape(destination)...),row)
                pos=string(alphabet[r],c)
                push!(destination_position,pos)
            end 
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


function slottingdict2dataframe(slotting::SlottingDict)
    labware_ids = Integer[]
    labware_names = String[]
    deck_position = String[]
    slot = Union{String,Integer}[]

    for s in keys(slotting)
        push!(labware_ids,location_id(s))
        push!(labware_names,name(s))
        push!(deck_position,slotting[s][1])
        push!(slot,slotting[s][2])
    end 
    
    return DataFrame(Rack=deck_position,Position=slot,LabwareID=labware_ids,Name=labware_names)
end 


"""
    nimbus(dispense_list::NimbusDispenseList,filepath::String)

Create Hamilton Nimbus dipsense instructions

  ## Arguments 
  * `dispense_list`: a NimbusDispenseList object encoding the design source and destination.
  * `filepath`: the ouput path of the dispense file in a .csv format

  ## Keyword Arguments
  * `nimbus_config`: A vector of NimbusRack objects that specify the configuration of the numbus. the default configuration is 5 50ml conical racks and a well plate rack. 
"""
function dispense(config::NimbusConfiguration,design::DataFrame, directory::AbstractString,protocol_name::AbstractString,labware::Vector{<:Labware},slotting::SlottingDict = slotting_greedy(labware,config); render_loading=true)
    # input error handling 
    W= nrow(design) 
    lws = length.(labware) 
    W == sum(lw) || ArgumentError("Dimension mismatch between design ($W) and number of wells in labware ($(sum(lw)))")
    all(map(x-> x in keys(slotting),labware) || ArgumentError("All labware must be slotted"))
    allunique(values(slotting)) || ArgumentError("Only one labware can be assigned to a given slot")
    is_square(design) || ArgumentError("design must be a square dataframe")

    df=convert_design(design,labware,slotting,config)
    n=nrow(df)
    dispense_df=DataFrame([[],[],[],[],[],],names(df))
    for i = 1:n 
        vol=df[i,"Volume (uL)"]
        while vol > 1e-6*u"µL" 
            shotvol=min(nozzles(head(config)).maxVol,vol) # maximum shot volume of 1 ml 
            push!(dispense_df,(df[i,"Source Labware ID"],df[i,"Source Position ID"],ustrip(uconvert(u"µL",shotvol)),df[i,"Destination Labware ID"],df[i,"Destination Position ID"]))
            vol-=shotvol
        end 
    end 
    append!(dispenses,dispense_df)

    n=nrow(dispenses)
    
    change_tip=zeros(Int64,n)
    for i in 2:n
        if dispenses[i,"Source Position ID"] != dispenses[i-1,"Source Position ID"]
            change_tip[i]=1
        end 
    end 

    windowsize=settings(config).max_tip_use
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


    if render_loading 
        loading = slottingdict2dataframe(slotting) 

        plt = plot(slotting, config)
        png(joinpath(full_dir,"loading.png"))

        CSV.write(joinpath(full_dir,"loading.csv"),loading)
    end 

    
    CSV.write(joinpath(full_dir,protocol_name*".csv"),dispenses)
    write(joinpath(full_dir,"config.json"),JSON.json(config))
    return nothing 
end 

#=
function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Nimbus;kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}
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
            push!(tt,nimbus(directory,protocol_name,des[s_idxs[i],d_idxs[j]],srcs[s_idxs[i]],destinations[d_idxs[j]],robot;kwargs...))
            push!(pn,protocol_name)
        end 
    end 
    return pn,tt

end 
=#
#= 
using JLDispense,JLD2,JLIMS

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
all_actives=JLD2.load("./src/Mixer/all_actives.jld2")["all_actives"]
#t,m=dispense_solver(all_actives,sources,nimbus_default;return_model=true)
protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",all_actives,sources,nimbus_default)
=#


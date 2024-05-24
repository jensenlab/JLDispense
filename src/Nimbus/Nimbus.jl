

struct NimbusDispenseList
    design::DataFrame
    source::String
    destination::String
end 




struct NimbusRack
    name::String
    size::Tuple{Int64,Int64}
    labware::ContainerName
    is_source::Bool
end 



# Define our available racks
#15 mL tube
TubeRack15ML_0001=NimbusRack("TubeRack15ML_0001",(4,6),containers[:conical_15ml],true)
# 50 mL tube
TubeRack50ML_0001=NimbusRack("TubeRack50ML_0001",(2,3),containers[:conical_50ml],true)
TubeRack50ML_0002=NimbusRack("TubeRack50ML_0002",(2,3),containers[:conical_50ml],true)
TubeRack50ML_0003=NimbusRack("TubeRack50ML_0003",(2,3),containers[:conical_50ml],true)
TubeRack50ML_0004=NimbusRack("TubeRack50ML_0004",(2,3),containers[:conical_50ml],true)
TubeRack50ML_0005=NimbusRack("TubeRack50ML_0005",(2,3),containers[:conical_50ml],true)
TubeRack50ML_0006=NimbusRack("TubeRack50ML_0006",(2,3),containers[:conical_50ml],true)
# 2 mL deep well plate
Cos_96_DW_2mL_0001=NimbusRack("Cos_96_DW_2mL_0001",(1,1),containers[:dWP96_2ml],false)
Cos_96_DW_2mL_0002=NimbusRack("Cos_96_DW_2mL_0002",(1,1),containers[:dWP96_2ml],false)


default_nimbus_config=[
    TubeRack50ML_0001,
    TubeRack50ML_0002,
     TubeRack50ML_0003,
     TubeRack50ML_0004,
     TubeRack50ML_0005,
     TubeRack50ML_0006,
    Cos_96_DW_2mL_0001
]

function convert_nimbus_design(design::DataFrame, source::String,destination::String;nimbus_config::Vector{NimbusRack}=default_nimbus_config,kwargs...)

    rack_codes=map(x->x.name,nimbus_config)


    if !in(source,rack_codes) || !in(destination,rack_codes)
        ArgumentError("Provide a valid Source Rack and Destination Rack Code")
    end 
    nimbus_config_dict=Dict(rack_codes .=> nimbus_config)

    dw=prod(nimbus_config_dict[destination].size)*prod(nimbus_config_dict[destination].labware.shape)
    if nrow(design) != dw
        ArgumentError("Number of rows in design must match the size of the destination")
    end 

    sw = prod(nimbus_config_dict[source].size)*prod(nimbus_config_dict[source].labware.shape)
    if ncol(design) != sw
        ArgumentError("Number of columns in design must match the size of the source")
    end 

    source_id=String[]
    source_position=Int64[]
    volume=Float64[]
    destination_id=String[]
    destination_position=String[]
    alphabet=collect('A':'Z')
    R,C=nimbus_config_dict[destination].labware.shape
    for col in 1:ncol(design)
        for row in 1:nrow(design)
            push!(source_id,source)
            push!(source_position,col)
            push!(volume,design[row,col])
            push!(destination_id,destination)

            c=cld(row,R)
            r=row-R*(c-1)
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
function nimbus(dispense_list::NimbusDispenseList, filepath::String;kwargs...)

    df=convert_nimbus_design(dispense_list.design,dispense_list.source,dispense_list.destination;kwargs...)
    n=nrow(df)
    dispense_df=DataFrame([[],[],[],[],[],],names(df))
    for i = 1:n 
        vol=df[i,"Volume (uL)"]
        while vol >0 
            dfrow=deepcopy(df[i,:])
            shotvol=min(1000,vol) # maximum shot volume of 1 ml 
            dfrow["Volume (uL)"]=shotvol
            push!(dispense_df,dfrow)
            vol-=shotvol
        end 
    end 

    n=nrow(dispense_df)
    
    change_tip=zeros(Int64,n)
    for i in 2:n
        if dispense_df[i,"Source Position ID"] != dispense_df[i-1,"Source Position ID"]
            change_tip[i]=1
        end 
    end 

    windowsize=25 
    for i in 1:n-windowsize+1
        window=change_tip[i:i+windowsize-1]
        if sum(window)==0
            change_tip[i+windowsize-1]=1
        end 
    end 
    dispense_df[!,"Change Tip Before"].= change_tip
    

    
    CSV.write(filepath,dispense_df)
end 


function nimbus(dispense_list::Vector{NimbusDispenseList}, filepath::String;kwargs...)
    designs=map(x->x.design,dispense_list)
    sources=map(x->x.source,dispense_list)
    destinations=map(x->x.destination,dispense_list)
    dfs=convert_nimbus_design.(designs,sources,destinations;kwargs...)
    df=vcat(dfs...)
    n=nrow(df)
    dispense_df=DataFrame([[],[],[],[],[],],names(df))
    for i = 1:n 
        vol=df[i,"Volume (uL)"]
        while vol >0 
            dfrow=deepcopy(df[i,:])
            shotvol=min(1000,vol)
            dfrow["Volume (uL)"]=shotvol
            push!(dispense_df,dfrow)
            vol-=shotvol
        end 
    end 

    n=nrow(dispense_df)
    
    change_tip=zeros(Int64,n)
    for i in 2:n
        if dispense_df[i,"Source Position ID"] != dispense_df[i-1,"Source Position ID"]
            change_tip[i]=1
        end 
    end 

    windowsize=25 
    for i in 1:n-windowsize+1
        window=change_tip[i:i+windowsize-1]
        if sum(window)==0
            change_tip[i+windowsize-1]=1
        end 
    end 
    dispense_df[!,"Change Tip Before"].= change_tip
    

    
    CSV.write(filepath,dispense_df)
end 


#= 
using DataFrames, JLDispense

design=100*ones(96,6)

dl=NimbusDispenseList(DataFrame(design,:auto),"TubeRack50ML_0001","Cos_96_DW_2mL_0001")
nimbus(dl,"/Users/BDavid/Desktop/test_nimbus.csv")

=#



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

const tempest_lv_nozzle = DiscreteNozzle(0.2u"µL",10u"µL",1.1)
const tempest_hv_nozzle = DiscreteNozzle(1u"µL",25u"µL",1.1)
abstract type TempestHead <: TransferHead end 

struct TempestLVHead <: TempestHead 
    channels::AbstractArray{Nozzle}
    TempestLVHead() = new(fill(tempest_lv_nozzle,8))
end 

struct TempestHVHead <: TempestHead 
    channels::AbstractArray{Nozzle}
    TempestHVHead() = new(fill(tempest_hv_nozzle,8))
end 

abstract type TempestDeckPosition <: DeckPosition end 

struct TempestInput <: TempestDeckPosition
    labware::Set{Type{<:Labware}}
    TempestInput() = new(Set([JLConstants.Bottle,JLConstants.Conical]))
end 

struct TempestMagazine <: TempestDeckPosition
    labware::Set{Type{<:Labware}}
    TempestMagazine() = new(Set([JLConstants.MicroPlate]))
end  

struct TempestSettings <: InstrumentSettings 
end 

tempest_deck= vcat(fill(TempestInput(),6),fill(TempestMagazine(),24))

TempestConfiguration=Configuration{TempestHead,Deck,TempestSettings}

const tempest_lv = TempestConfiguration("Tempest Low Volume",TempestLVHead(),tempest_deck,TempestSettings())
const tempest_hv = TempestConfiguration("Tempest High Volume",TempestHVHead(),tempest_deck,TempestSettings())



function masks(h::TempestHead,l::JLConstants.WP96) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 8
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,12
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # tempest cannot aspirate from well plates 
        return false
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == c && wj == pn 
    end 
    return Ma,Md
end 

function masks(h::TempestHead,l::JLConstants.WP384) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 8
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 2,24
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # tempest cannot aspirate from well plates 
        return false
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi = 2*(c-1)+pm && wj == pn 
    end 
    return Ma,Md,(0,0),(Pi*Pj,C)
end 


function masks(h::TempestHead,l::Union{JLConstants.Bottle,JLConstants.Tube}) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # tempest cannot aspirate from well plates 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        return false 
    end 
    return Ma,Md,(Pi*Pj,C),(0,0)
end 





  
    const tempest_names=Dict{Type{<:Labware},String}(
    WP96=>"PT3-96-Assay.pd.txt",
    WP384=>"PT9-384-Assay.pd.txt")


    function convert_design(design::DataFrame,labware::Vector{<:Labware},slotting::SlottingDict,config::TempestConfiguration)

        all(map(x-> x[1] in deck(config),values(slotting))) || ArgumentError("All deck positions in the SlottingDict must be present on the deck")
        source_cols=[]
        dispense_rows=[]
        source_names= []
        for col in 1:ncol(design)
            lw=get_labware(labware,col)
            if can_aspirate(head(config),slotting[lw],lw) 
                push!(source_cols,col)
                push!(source_names,JLIMS.name(lw))
            end 
        end 
        for row in 1:nrow(design) 
            lw=get_labware(labware,row)
            if can_dispense(head(config),slotting[lw],lw)
                push!(dispense_rows,row)
            end 
        end 
    
        out_design = design[dispense_rows,source_cols]
    
        dest_labware= filter(x-> can_dispense(head(config),slotting[x],x),labware)
        allequal(map(x->typeof(x),dest_labware)) || error("all destination labware must be the same type for a multidispense list tempest run")
        if length(dest_labware) > 1 
            dest_lw = get_labware.((dest_labware,),nrow(out_design))
            outs= DataFrame[]

            for l in dest_labware
                idxs=findall(x->x==l,dest_lw)
                push!(outs,out_design[idxs,:])
            end
            return outs ,dest_labware, fill(source_names,length(dest_labware))
        else
            return out_design, dest_labware[1],source_names
        end
    end 
    
function write_tempest_dl(dispenses::DataFrame,labware::Labware, filename::String,source_names::Vector{String})
    n_stocks=ncol(dispenses)
    outfile=open(filename,"w")
    platefilename=tempest_names[typeof(labware)]
    if platefilename==""
        error(ArgumentError("Selected plate type $(typeof(labware)) is not supported by Tempest."))
    end 
    print(outfile,join(["Version            :", 6],'\t'),"\r\n")
    print(outfile,join(["Plate type name    :", platefilename],'\t'),"\r\n")
    print(outfile,join(["Priority Delays    :", 0],'\t'),"\r\n")
    R,C=JLIMS.shape(labware)

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
end

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
    dispense(config::TempestConfiguration,design::DataFrame, directory::AbstractString,protocol_name::AbstractString,destinations::Labware)

Create Tempest multidipsense instructions for SBS plate dispensing

  ## Arguments 
  * `config`: A Configuration object defining the tempest 
  * `design`: A (# of total runs) x (# of stocks) DataFrame containing the volume of each stock for each run in µL. 
  * `directory`: the ouput directory of the dispense file. The function automatically strips any incorrect file extension and adds the appropriate one.
  * `protocol_name`: The name of the multidispense list file. The funcion outputs the multidispense list in directory
  * `destination`: A  destination labware. Can either be a JLIMS.Labware object or A JLIMS.Labware subtype 
"""
function dispense(config::TempestConfiguration,design::DataFrame,directory::AbstractString,protocol_name::AbstractString,labware::Vector{<:Labware},slotting::SlottingDict=slotting_greedy(labware,config))

    dispenses, destinations, source_names =convert_design(design,labware,slotting,config)
    if ~isdir(full_dir)
        mkdir(full_dir)
    end
    N= length(dispenses)
    full_dir=joinpath(directory,protocol_name)
    if N == 1  # single dispense file 
        filename=joinpath(full_dir,protocol_name*".dl.txt")
        write_tempest_dl(dispenses,destinations,filename,source_names)
    else  # multidispense 

        basenames = [protocol_name*i for i in 1:N]
        for i in 1:N 
            filename=joinpath(full_dir,basenames[i]*".dl.txt")
            write_tempest_dl(dispenses[i],destinations[i],filename,source_names[i])
        end 
        mdlbasename=protocol_name*".mdl.txt"
        outpath=joinpath(full_dir,mdlbasename)
        outfile=open(outpath,"w")
        print(outfile,join(["LP:$(tempest_names[destinations[1]]).pd.txt"],'\t'),"\r\n") # we previously check that all destinations are the same type, so we can safely take the first one
        for i in 1:n
            print(outfile,join(["P:$i","TP:1","DL:$(basenames[i])"],'\t'),"\r\n")
        end
        close(outfile)
    end
    
    #delay_header=vcat(n_stocks,repeat([0,""],n_stocks))

    write(joinpath(full_dir,"config.json"),JSON.json(config))
    return nothing 
end 



#=
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
=#
#= 

using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,tempest_default)
=# 
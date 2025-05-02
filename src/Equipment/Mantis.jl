





#=
struct DiscreteNozzle <: Nozzle
    minVol::Unitful.Volume
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
end 
=#
const mantis_lv_nozzle = DiscreteNozzle(0.1u"µL",10u"µL",1.1)
const mantis_hv_nozzle = DiscreteNozzle(1u"µL",25u"µL",1.1)

abstract type MantisHead <: TransferHead end 

struct MantisLVHead <: MantisHead 
    channels::AbstractArray{Nozzle}
    MantisLVHead() = new([mantis_lv_nozzle])
end 

struct MantisHVHead <: MantisHead 
    channels::AbstractArray{Nozzle}
    MantisHVHead() =new([mantis_hv_nozzle])
end 

struct MantisSettings <:InstrumentSettings
    MantisSettings() = new()
end 

    
const lc3 = StackPosition("LC3",Set([JLConstants.TipReservior,JLConstants.Conical15]),(24,1),true,false,false,false,circle)
const mantis_hv_pos = StackPosition("Mantis High Volume" , Set([JLConstants.Conical]),(4,1),true,false,false,false,circle)
const mantis_main = SBSPosition("Mantis Center Position",Set([JLConstants.MicroPlate,JLConstants.BreakawayPCRPlate]),(1,1),false,true,false,false,rectangle)

const mantis_deck = vcat(mantis_hv_pos,mantis_main,lc3)


MantisConfiguration = Configuration{MantisHead,Deck,MantisSettings}

const mantis_lv = MantisConfiguration("Mantis Low Volume",MantisLVHead(),mantis_deck,MantisSettings())
const mantis_hv = MantisConfiguration("Mantis High Volume",MantisHVHead(),mantis_deck,MantisSettings())

### define deck access functions 

function plumbing_mask(h::MantisHead)
    pistons = 1
    channels = 1
    function Mp(p::Integer,c::Integer)
        1 <= p <= pistons || return false 
        1 <= c <= channels || return false 
      return p == c 
    end
    return Mp ,pistons
  end 
  


function masks(h::MantisHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = Wi,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # Mantis cannot aspirate from well plates 
        return false
    end 
    function Md(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == pm && wj == pn 
    end 
    return Ma,Md, (0,0),(Pi*Pj,C)
end 

function masks(h::MantisLVHead,l::Union{JLConstants.Conical15,JLConstants.TipReservior}) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = JLIMS.shape(l)
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # Mantis cannot aspirate from well plates 
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

function masks(h::MantisHVHead,l::JLConstants.Conical50) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) # Mantis cannot aspirate from well plates 
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



function convert_design(design::DataFrame,labware::Vector{<:Labware},slotting::SlottingDict,config::MantisConfiguration)

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

    if length(dest_labware) > 1 
        error("multiple destination labware detected. Mantis only supports a single destination labware")
    end 

    return out_design, dest_labware[1],source_names
end 

        




"""
    dispense(config::MantisConfiguration, design::DataFrame,directory::AbstractString,protocol_name::AbstractString,destination::Labware)

Create Mantis dipsense instructions for SBS plate dispensing

## Arguments 
* `config`: A Configuration object defining the Mantis
* `design`: a (# of sources) x (# of destinations) dataframe containing the volume of each transfer in µL.
* `directory`: output file directory
* `protocol_name`: protocol name  
* `destinations`: A JLIMS.Labware destination object  
"""
function dispense(config::MantisConfiguration, design::DataFrame,directory::AbstractString,protocol_name::AbstractString,labware::Vector{<:Labware},slotting::SlottingDict = slotting_greedy(labware,config);render_loading=true)
        dispenses,destination,source_names= convert_design(design,labware,slotting,config)
        full_dir=joinpath(directory,protocol_name)
        if ~isdir(full_dir)
            mkdir(full_dir)
        end 
        filename=joinpath(full_dir,protocol_name*".dl.txt")
        n_stocks=ncol(dispenses)
        delay_header=vcat(n_stocks,repeat([0,""],n_stocks))
        outfile=open(filename,"w")
        platefilename=mantis_names[typeof(destination)]
        print(outfile,join(["[ Version: 5 ]"],'\t'),"\r\n")
        print(outfile,platefilename,"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")
        print(outfile,join([1],'\t'),"\r\n")
        print(outfile,join(delay_header,'\t'),"\r\n")

        for i in 1:n_stocks 
            vols=ustrip.(uconvert((u"µL",),design[:,i]))
            vols=round.(reshape(vols,R,C),digits=1)
            print(outfile,join([source_names[i],"","Normal"],'\t'),"\r\n")
            print(outfile,join(["Well",1],'\t'),"\r\n")
            for r in 1:R
                print(outfile,join(vols[r,:],'\t'),"\r\n")
            end 
        end 
        close(outfile)
        write(joinpath(full_dir,"config.json"),JSON.json(robot))
        if render_loading 
            plt = plot(slotting,config)
            png(joinpath(full_dir,"loading.png"))
        end 
        return nothing 
end 






#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]

protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,mantis_default)

=# 


const platemaster_nozzle = ContinuousNozzle(2u"µL",200u"µL",200u"µL",0u"µL",1,false)



struct PlateMasterHead<: FixedTransferHead
    channels::AbstractArray{Nozzle}
    PlateMasterHead() = new(fill(platemaster_nozzle,8,12))
end 

plate_master_compat_labware=Set([JLConstants.WellPlate])

pm1=SBSPosition("Slot 1", plate_master_compat_labware,(1,1),false,true,false,false,"rectangle")
pm2=SBSPosition("Slot 2", plate_master_compat_labware,(1,1),false,true,false,false,"rectangle")
pm3=SBSPosition("Slot 3", plate_master_compat_labware,(1,1),false,true,false,false,"rectangle")
pm4=SBSPosition("Slot 4", plate_master_compat_labware,(1,1),false,true,false,false,"rectangle")


const platemaster_deck= [pm1 pm2; pm3 pm4]
struct PlateMasterSettings <: InstrumentSettings 
    setting::String
    PlateMasterSettings()=new("placeholder")
end 


PlateMasterConfiguration = Configuration{PlateMasterHead,Deck,PlateMasterSettings}

const platemaster =PlateMasterConfiguration("PlateMaster",PlateMasterHead(),platemaster_deck,PlateMasterSettings())



function plumbing_mask(h::PlateMasterHead)
    pistons = 1
    channels = 96 
    function Mp(p::Integer,c::Integer)
        1 <= p <= pistons || return false 
        1 <= c <= channels || return false 
      return true
    end
    return Mp ,pistons
  end 


function masks(h::PlateMasterHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    Ci,Cj=8,12
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    C=falses(Ci,Cj)
    S = (Pi*Pj,Ci*Cj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= Ci*Cj || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        ci,cj=cartesian(C,c)
        return wi == ci && wj == cj 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 
function masks(h::PlateMasterHead,l::JLConstants.DeepReservior) # for generic 96 well plates, we will define a separate method for 384 well plates 
    Ci,Cj=8,12
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    C=falses(Ci,Cj)
    S = (Pi*Pj,Ci*Cj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= Ci*Cj || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        ci,cj=cartesian(C,c)
        return true 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 



function masks(h::PlateMasterHead,l::JLConstants.WP384) 
    Ci,Cj=8,12
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 2,2
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    C=falses(Ci,Cj)
    S=(Pi*Pj,Ci*Cj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels 
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= Ci*Cj || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        ci,cj=cartesian(C,c)
        return wi == 2*(ci-1)+pm  && wj  == 2*(cj-1)+pn
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 




function dispense(config::PlateMasterConfiguration, design::DataFrame, directory::AbstractString,protocol_name::AbstractString,labware::Vector{<:Labware},slotting::SlottingDict=slotting_greedy(labware,config);render_loading=true,kwargs...) 
    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
      mkdir(full_dir)
    end 
    CSV.write(joinpath(full_dir,"dispenses.csv"),design)
    write(joinpath(full_dir,"config.json"),JSON.json(config))
    #print("entire $(experiment_name) folder must be moved to Dropbox -> JensenLab -> Cobra")
    if render_loading 
      plt = plot(slotting,config)
      png(joinpath(full_dir,"loading.png"))
    end 
    return nothing 
end


    


    

    



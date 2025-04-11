

const platemaster_nozzle = ContinuousNozzle(2u"µL",200u"µL",200u"µL",0u"µL",1,false)



struct PlateMasterHead<: FixedTransferHead
    nozzles::AbstractArray{Nozzle}
    PlateMasterHead() = new(fill(platemaster_nozzle,8,12))
end 

struct PlateMasterDeckPosition <: DeckPosition 
    labware::Set{Type{<:Labware}}
    PlateMasterDeckPosition() = new(Set([JLConstants.WellPlate]))
end


const platemaster_deck= fill(PlateMasterDeckPosition(),4)
struct PlateMasterSettings <: InstrumentSettings 
    PlateMasterSettings()= new()
end 


PlateMasterConfiguration = Configuration{PlateMasterHead,Deck{PlateMasterDeckPosition},PlateMasterSettings}

const platemaster =PlateMasterConfiguration(PlateMasterHead(),platemaster_deck,PlateMasterSettings())

function can_aspirate(h::PlateMasterlHead, d::PlateMasterDeckPosition,l::Labware) 
    return can_place(l,d)
end
  
function can_dispense(h::PlateMasterHead,d::PlateMasterDeckPosition,l::Labware) 
    return can_place(l,d)
end




function masks(h::PlateMasterHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    Ci,Cj=1,1
    Wi,Wj=shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    C=falses(Ci,Cj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= Ci*Cj || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        ci,cj=cartesian(C,c)
        return wi == ci  && wj == cj
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 



function masks(h::PlateMasterHead,l::JLConstants.WP384) 
    Ci,Cj=1,1
    Wi,Wj=shape(l)
    Pi,Pj = 2,2
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    C=falses(Ci,Cj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels 
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= Ci*Cj || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        ci,cj=cartesian(C,c)
        return wi % 2 == 2-pm  && wj % 2  == 2-pn
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 





    


    

    



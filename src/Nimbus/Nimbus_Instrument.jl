
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
abstract type NimbusHead <: TransferHead end 

struct NimbusSingleChannelHead <: NimbusHead 
    nozzles::Nozzle
end 


const nimbus_nozzle= Nozzle(50u"µL",1000u"µL",1000u"µL",25u"µL",1,false,false)


const nimbus_single_head = NimbusSingleChannelHead(nimbus_nozzle)



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




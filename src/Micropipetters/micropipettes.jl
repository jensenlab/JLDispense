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
const p1000_nozzle = Nozzle(100u"µL",1000u"µL",1000u"µL",0u"µL",1,false,false)
const p200_nozzle = Nozzle(20u"µL",200u"µL",200u"µL",0u"µL",1,false,false)
const p20_nozzle = Nozzle(2u"µL",20u"µL",20u"µL",0u"µL",1,false,false)
const p2_nozzle = Nozzle(0.2u"µL",2u"µL",2u"µL",0u"µL",1,false,false)
const p100_nozzle = Nozzle(10u"µL",100u"µL",100u"µL",0u"µL",1,false,false)
const p10_nozzle = Nozzle(1u"µL",10u"µL",10u"µL",0u"µL",1,false,false)



struct SingleChannelHead <: FixedTransferHead
    nozzles::Nozzle
end 

const p1000_head =SingleChannelHead(p1000_nozzle)
const p200_head =SingleChannelHead(p200_nozzle)
const p20_head =SingleChannelHead(p20_nozzle)
const p2_head =SingleChannelHead(p2_nozzle)

function get_masks(h::SingleChannelHead,l::Labware)
    C= 1
    Wi,Wj=shape(l)
    Pi,Pj =shape(l)
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
    Md = deepcopy(Ma) 
    return Ma,Md
end 





struct EightChannelHead <: FixedTransferHead 
    nozzles::Vector{Nozzle}
end 

const p1000_multichannel_head =EightChannelHead( fill(p1000_nozzle,8))
const p100_multichannel_head = EightChannelHead(fill(p100_nozzle,8))
const p10_multichannel_head = EightChannelHead(fill(p10_nozzle,8))



function masks(h::EightChannelHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= length(nozzles(h))
    Wi,Wj=shape(l)
    Pi,Pj = 1,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == c && wj == pn 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 

function masks(h::EightChannelHead,l::JLConstants.WP384) 
    C= length(nozzles(h))
    Wi,Wj=shape(l)
    Pi,Pj = 2,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == 2*(c-1)+pm && wj == pn 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 






    
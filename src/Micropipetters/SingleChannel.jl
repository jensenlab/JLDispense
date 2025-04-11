#=
struct ContinuousNozzle
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    multidispense::Bool
end=#
const p1000_nozzle = ContinuousNozzle(100u"µL",1000u"µL",1000u"µL",0u"µL",1,false)
const p200_nozzle = ContinuousNozzle(20u"µL",200u"µL",200u"µL",0u"µL",1,false)
const p20_nozzle = ContinuousNozzle(2u"µL",20u"µL",20u"µL",0u"µL",1,false)
const p2_nozzle = ContinuousNozzle(0.2u"µL",2u"µL",2u"µL",0u"µL",1,false)
const p100_nozzle = ContinuousNozzle(10u"µL",100u"µL",100u"µL",0u"µL",1,false)
const p10_nozzle = ContinuousNozzle(1u"µL",10u"µL",10u"µL",0u"µL",1,false)



struct SingleChannelHead <: FixedTransferHead
    nozzles::Nozzle
end 

const p1000_head =SingleChannelHead(p1000_nozzle)
const p200_head =SingleChannelHead(p200_nozzle)
const p20_head =SingleChannelHead(p20_nozzle)
const p2_head =SingleChannelHead(p2_nozzle)



struct SingleChannelSettings <: InstrumentSettings 
end 

const single_channel_settings= SingleChannelSettings()



const single_channel_deck = Deck[UnconstrainedDeckPosition(),UnconstrainedDeckPosition()]

SingleChannelConfiguration = Configuration{SingleChannelHead,Deck{UnconstrainedDeckPosition},SingleChannelSettings}


const p1000 = SingleChannelConfiguration(p1000_head,single_channel_deck,single_channel_settings)
const p200 = SingleChannelConfiguration(p200_head,single_channel_deck,single_channel_settings)
const p20 = SingleChannelConfiguration(p20_head,single_channel_deck,single_channel_settings)
const p2 = SingleChannelConfiguration(p2_head,single_channel_deck,single_channel_settings)

### define deck access functions 


function can_aspirate(h::SingleChannelHead, d::UnconstrainedDeckPosition,l::Labware) 
    return can_place(l,d)
end
function can_dispense(h::SingleChannelHead,d::UnconstrainedDeckPosition,l::Labware) 
    return can_place(l,d)
end



function masks(h::SingleChannelHead,l::JLIMS.Labware)
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



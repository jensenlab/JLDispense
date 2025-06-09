
const p100_nozzle = ContinuousNozzle(10u"µL",100u"µL",100u"µL",0u"µL",1,false)
const p10_nozzle = ContinuousNozzle(1u"µL",10u"µL",10u"µL",0u"µL",1,false)

struct EightChannelHead <: FixedTransferHead 
    channels::AbstractArray{Nozzle}
end 

const p1000_multichannel_head =EightChannelHead(fill(p1000_nozzle,8))
const p100_multichannel_head = EightChannelHead(fill(p100_nozzle,8))
const p10_multichannel_head = EightChannelHead(fill(p10_nozzle,8))

struct EightChannelSettings <: InstrumentSettings 
end 


const eight_channel_settings= EightChannelSettings()

const eight_channel_labware = Set([JLConstants.WellPlate])
const eight_channel_position = SBSPosition("Eight_Channel_Position", eight_channel_labware,(4,6),true,true,false,false,"rectangle")

const eight_channel_deck = [eight_channel_position]

EightChannelConfiguration = Configuration{EightChannelHead,Deck{SBSPosition},EightChannelSettings}


const multichannel_p1000 = EightChannelConfiguration("P-1000 Multichannel",p1000_multichannel_head,eight_channel_deck,eight_channel_settings)
const multichannel_p100 = EightChannelConfiguration("P-100 Multichannel",p100_multichannel_head,eight_channel_deck,eight_channel_settings)
const multichannel_p10 = EightChannelConfiguration("P-10 Multichannel",p10_multichannel_head,eight_channel_deck,eight_channel_settings)
# deck access functions
function can_aspirate(h::EightChannelHead, d::DeckPosition,l::Labware) 
    return can_place(l,d)
end 
function can_dispense(h::EightChannelHead,d::DeckPosition,l::Labware) 
    return can_place(l,d)
end


function plumbing_mask(h::EightChannelHead)
    pistons = 1
    channels = 8
    function Mp(p::Integer,c::Integer)
        1 <= p <= pistons || return false 
        1 <= c <= channels || return false 
      return true
    end
    return Mp,pistons
  end 


function masks(h::EightChannelHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 8
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    S = (Pi*Pj,C)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return  wj == pn && wi == c
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 

function masks(h::EightChannelHead,l::JLConstants.WP384) 
    C= 8
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 2,Wj
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    S = (Pi*Pj,C)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == 2*(c-1) + pm && wj == pn 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 



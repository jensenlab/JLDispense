
const p100_nozzle = ContinuousNozzle(10u"µL",100u"µL",100u"µL",0u"µL",1,false)
const p10_nozzle = ContinuousNozzle(1u"µL",10u"µL",10u"µL",0u"µL",1,false)

abstract type EightChannelHead <: FixedTransferHead end 


    # Eight channel with head placed in a vertical orientation (each channel corresponds to a row)
struct EightChannelHeadVertical <: EightChannelHead
    channels::AbstractArray{Nozzle}
end 

const p1000_multichannel_head_v =EightChannelHeadVertical(fill(p1000_nozzle,8))
const p100_multichannel_head_v = EightChannelHeadVertical(fill(p100_nozzle,8))
const p10_multichannel_head_v = EightChannelHeadVertical(fill(p10_nozzle,8))


 # Eight channel with head placed in a horizontal orientation (each channel corresponds to a column )
struct EightChannelHeadHorizontal <: EightChannelHead
    channels::AbstractArray{Nozzle}
end 

const p1000_multichannel_head_h =EightChannelHeadHorizontal(fill(p1000_nozzle,8))
const p100_multichannel_head_h = EightChannelHeadHorizontal(fill(p100_nozzle,8))
const p10_multichannel_head_h = EightChannelHeadHorizontal(fill(p10_nozzle,8))

struct EightChannelSettings <: InstrumentSettings 
end 


const eight_channel_settings= EightChannelSettings()

const eight_channel_labware = Set([JLConstants.WellPlate])
const eight_channel_position = SBSPosition("Eight_Channel_Position", eight_channel_labware,(4,6),true,true,false,false,"rectangle")

const eight_channel_deck = [eight_channel_position]

EightChannelConfiguration = Configuration{EightChannelHead,Deck{SBSPosition},EightChannelSettings}


const multichannel_p1000_v = EightChannelConfiguration("P-1000 Multichannel (Vertical)",p1000_multichannel_head_v,eight_channel_deck,eight_channel_settings)
const multichannel_p100_v = EightChannelConfiguration("P-100 Multichannel (Vertical)",p100_multichannel_head_v,eight_channel_deck,eight_channel_settings)
const multichannel_p10_v = EightChannelConfiguration("P-10 Multichannel (Vertical)",p10_multichannel_head_v,eight_channel_deck,eight_channel_settings)
const multichannel_p1000_h = EightChannelConfiguration("P-1000 Multichannel (Horizontal)",p1000_multichannel_head_h,eight_channel_deck,eight_channel_settings)
const multichannel_p100_h = EightChannelConfiguration("P-100 Multichannel (Horizontal)",p100_multichannel_head_h,eight_channel_deck,eight_channel_settings)
const multichannel_p10_h = EightChannelConfiguration("P-10 Multichannel (Horizontal)",p10_multichannel_head_h,eight_channel_deck,eight_channel_settings)
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


function masks(h::EightChannelHeadVertical,l::JLConstants.WP96) # for generic 96 well plates, we will define a separate method for 384 well plates 
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

function masks(h::EightChannelHeadVertical,l::JLConstants.WP384) 
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


function masks(h::EightChannelHeadHorizontal,l::JLConstants.WP96) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 8
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = Wi, 5
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
        return  pm == wi && wj == pn +c -1 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 

function masks(h::EightChannelHeadHorizontal,l::JLConstants.WP384) 
    C= 8
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = Wi,10
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
        return wi ==  pm  && wj == pn + 2*(c-1) 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end


function masks(h::EightChannelHead,l::JLConstants.DeepReservior) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C = 8 
    Wi,Wj=JLIMS.shape(l)
    Pi,Pj = 1,1
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    S = (Pi*Pj,C)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C|| return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return true 
    end 
    Md = deepcopy(Ma) 
    return Ma,Md,S,S
end 

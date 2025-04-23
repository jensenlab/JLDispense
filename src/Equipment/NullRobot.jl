#=
struct ContinuousNozzle
    minVol::Unitful.Volume
    maxVol::Unitful.Volume
    maxAsp::Unitful.Volume 
    deadVol::Unitful.Volume 
    deadVolFactor::Real 
    multidispense::Bool
end=#
const unconstrained_nozzle = ContinuousNozzle(0u"µL",Inf*u"µL",Inf*u"µL",0u"µL",1,false)


const unconstrained_head =SingleChannelHead(unconstrained_nozzle)


const null_robot = SingleChannelConfiguration("Null Robot",unconstrained_head,single_channel_deck,single_channel_settings)











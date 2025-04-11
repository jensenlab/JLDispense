
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
abstract type TempestHead <: TransferHead end 

struct TempestLVHead <: TempestHead 
    nozzles::Vector{Nozzle}
end 

struct TempestHVHead <: TempestHead 
    nozzles::Vector{Nozzle}
end 

const tempest_lv_nozzle = Nozzle(0.2u"µL",10u"µL",1.1)
const tempest_hv_nozzle = Nozzle(1u"µL",25u"µL",1.1)


const tempest_lv_head = TempestLVHead(tempest_lv_nozzle)
const tempest_hv_head = TempestHVHead(tempest_hv_nozzle)




function masks(h::TempestHead,l::JLConstants.WP96) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 8
    Wi,Wj=shape(l)
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
    Wi,Wj=shape(l)
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
    return Ma,Md
end 


function masks(h::TempestHead,l::Union{JLConstants.Bottle,JLConstants.Tube}) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= 1
    Wi,Wj=shape(l)
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
    return Ma,Md
end 

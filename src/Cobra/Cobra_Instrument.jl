

struct CobraHead <: TransferHead
    nozzles::AbstractArray{Nozzle}
end 

const cobra_nozzle = Nozzle(0.3u"µL",40u"µL",800u"µL",20u"µL",1.1,false,true)


const cobra_head = CobraHead(fill(cobra_nozzle,4,1))


function masks(h::CobraHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    C= length(nozzles(h))
    Wi,Wj=shape(l)
    Pi,Pj = 5,12
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Ma(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == c+pm-1 && wj == pn 
    end 
    C= length(nozzles(h))
    Wi,Wj=shape(l)
    Pi,Pj = 11,12
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Md(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == c+pm-4 && wj == pn 
    end 
    return Ma,Md
end 

function masks(h::CobraHead,l::JLConstants.WP384) 
    C= length(nozzles(h))
    Wi,Wj=shape(l)
    Pi,Pj = 10,12
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
    C= length(nozzles(h))
    Wi,Wj=shape(l)
    Pi,Pj = 22,12
    W = falses(Wi,Wj)
    P=falses(Pi,Pj)
    function Md(w::Integer,p::Integer,c::Integer) 
        # w=wells, p=positions, c=channels
        1 <= w <= Wi*Wj || return false 
        1 <= p <= Pi*Pj || return false 
        1 <= c <= C || return false 
        wi,wj=cartesian(W,w)
        pm,pn=cartesian(P,p)
        return wi == 2*(c-1)+pm-6 && wj == pn 
    end 
    return Ma,Md
end 





    


    

    











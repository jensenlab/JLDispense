

const platemaster_nozzle = Nozzle(2u"µL",200u"µL",200u"µL",0u"µL",1,false,false)



struct PlateMasterHead<: FixedTransferHead
    nozzles::AbstractArray{Nozzle}
    PlateMasterHead() = new(fill(platemaster_nozzle,8,12))
end 


function masks(h::PlateMasterHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
    Ci,Cj=8,12
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
    Ci,Cj=8,12
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
        return wi == 2*(ci-1)+pm  && wj == 2*(cj-1)+pn
    end 
    Md = deepcopy(Ma) 
    return Ma,Md
end 





    


    

    



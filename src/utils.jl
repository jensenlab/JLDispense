

function cartesian(x::AbstractArray,i)
    return Tuple(CartesianIndices(x)[i])
end 

function linear(x::AbstractArray,idxs...)
    return LinearIndeices(x)[idxs...]
end 


function transfer_table(source::Labware,destination::Labware,design::DataFrame)
    transfer_table=DataFrame(Source=Integer[],Destination=Integer[],Quantity=Real[],Unit=AbstractString[])
    r=nrow(design)
    c=ncol(design)
    for col in 1:c
        for row in 1:r 
            val=design[row,col]
            quantity=ustrip(val)
            if quantity==0 
                continue 
            else 
                source=JLIMS.location_id(source[cartesian(source,row)])
                destination=JLIMS.location_id(destination[cartesian(destination,col)])
                un=string(unit(val))
                push!(transfer_table,(source,destination,quantity,un))
            end 
        end 
    end 
    return transfer_table
end 

function is_square(df::DataFrame)
    r,c=size(df)
    return r == c 
end 


function get_labware(labware::Vector{Labware},index::Integer) 

    return vcat(map(x->fill(x,length(x)),labware)...)[index] 
end 

function rack_codes(n) 
    codes=["" for _ in 1:n]
    alphabet=collect('A':'Z')
    k=length(alphabet)
    i=1
    alphacounter=0
    for i = 1:n
        codes[i]=repeat(alphabet[mod(i-1,k)+1],cld(i,k))
    end 
    return codes 
end 


function circle(x,y,d)
    r=d/2
    θ=LinRange(0,2*π,500)
    x .+ r*sin.(θ), y .+ r*cos.(θ)
end 

function rectangle(x,y,w,h= w )
    a = w/2
    b = h/2 
    up = LinRange(y-b,y+b,500)
    over = LinRange(x-a,x+a,500)
    i = length(up)
    j= length(over)
    arg1 = vcat(over,fill(x+a,i),reverse(over),fill(x-a,i))
    arg2=vcat(fill(y-b,j),up,fill(y+b,j),reverse(up))
    return arg1,arg2
    
end 


plotting_shape(::Labware) = rectangle 
plotting_shape(::JLConstants.Tube) =circle


function slotting_greedy(labware::Vector{<:Labware},config::Configuration)

    slotting = SlottingDict()
    all_slots = Set{Tuple{DeckPosition,Int}}()

    for position in deck(config) 
        n_slots = prod(position.slots)
        for i in 1:n_slots 
            push!(all_slots,(position,i))
        end 
    end 

    for lw in labware 
        for s in all_slots 
            if can_place(lw,s[1])
                slotting[lw]= s 
                delete!(all_slots,s)
                break 
            end 
        end 
        if !in(lw,keys(slotting))
            error("cannot find an open slot for $lw, use differnt labware or change the configuration")
        end 
    end 
    return slotting
end 






function plot(slotting::SlottingDict,config::Configuration;wrapwidth::Integer=20,fontsize::Integer=14)
    positions = unique(map(x->x[1],values(slotting)))
    n=length(positions)

    

    plot_shape=size(deck(config))

    plts = [] 
    for pos in positions 
        r,c=slots(pos)
        k = r*c 
        lw= filter(x-> slotting[x][1]==pos,keys(slotting))
        vals = fill("",r,c)
        plotfun = Matrix{Function}(fill(circle,r,c))
        for l in lw 
            idx = slotting[l][2]
            vals[idx]= JLIMS.name(l)
            plotfun[idx]=plotting_shape(l)
        end 
        plt=plot(grid=false,size=(1200,800),yflip=true,legend=false,dpi=300,xticks=1:c,yticks=(1:r,rack_codes(r)),xmirror=true,tickdirection=:none,tickfontsize=18)

        for x in 1:c
            for y in 1:r 
                annotate!(x,y,text(TextWrap.wrap(vals[y,x],width=wrapwidth),:center,fontsize))
                plot!(plotfun[y,x](x,y,1),color="black")
            end 
        end 

        plot!(ylims=(0.5,r+0.5),xlims=(0.5,c+0.5))
        plot!(title=pos.name,titlefontsize=20)
        push!(plts,plt) 
    end





    plot(plts...,layout=plot_shape,legend=false,title=typeof(config))

end 

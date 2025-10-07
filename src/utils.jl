

function cartesian(x::AbstractArray,i)
    return Tuple(CartesianIndices(x)[i])
end 

function linear(x::AbstractArray,idxs...)
    return LinearIndices(x)[idxs...]
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
                src_well = children(source)[cartesian(children(source),row)...]
                src=JLIMS.location_id(src_well)
                dst_well = children(destination)[cartesian(children(destination),col)...]
                dst=JLIMS.location_id(dst_well)
                un=string(unit(val))
                push!(transfer_table,(src,dst,quantity,un))
            end 
        end 
    end 
    return transfer_table
end 

function is_square(df::DataFrame)
    r,c=size(df)
    return r == c 
end 


function get_labware(labware::Vector{<:Labware},index::Integer) 

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

function rectangle(x,y,w)
    h=w
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


#plotting_shape(::Labware) = rectangle 
#plotting_shape(::JLConstants.Tube) = circle




plotting_fun=Dict(
    "rectangle"=>rectangle,
    "circle"=>circle
)



function plot(slotting::SlottingDict,config::Configuration;wrapwidth::Integer=20,titlefontsize::Integer=18,fontsize::Integer=14,plotsize=(1200,800))

    plot_shape=(0,0)
    if deck(config) isa Vector
        plot_shape=(length(deck(config)),1)
    else

        plot_shape=size(deck(config))
    end 
    plts = [] 
    for row in 1:plot_shape[1]
        for col in 1:plot_shape[2]

            pos = deck(config)[row,col]
            plt=plot(pos,titlefontsize=titlefontsize,tickfontsize=titlefontsize)
            r,c = 1,1
            if !isa(pos,EmptyPosition)
                r,c=slots(pos)
                k = r*c 
                lw= filter(x-> slotting[x][1]==pos,keys(slotting))
                vals = fill("",r,c)
                plotfun = Matrix{Function}(fill(circle,r,c))
                for l in lw 
                    idx = slotting[l][2]
                    vals[idx]= JLIMS.name(l)
                end 
                

                for x in 1:c
                    for y in 1:r 
                        annotate!(x,y,text(TextWrap.wrap(vals[y,x],width=wrapwidth),:center,fontsize))
                        plot!((plotting_fun[pos.plotting_shape])(x,y,1),color="black")
                    end 
                end 
            end 
            

            plot!(ylims=(0.5,r+0.5),xlims=(0.5,c+0.5))
            plot!(title=TextWrap.wrap(name(pos),width=wrapwidth),titlefontsize=titlefontsize)
            push!(plts,plt) 
        end
    end 
    finalsize = plotsize .* plot_shape
    plot(plts...,size=finalsize,layout=plot_shape,legend=false,plot_title=TextWrap.wrap(name(config),width=wrapwidth),margins=(10,:mm),titlefontsize=titlefontsize)

end 


function plot(pos::SBSPosition;titlefontsize=18,tickfontsize=18) 
    r,c = slots(pos)
    plot(grid=false,size=(1200,800),yflip=true,legend=false,dpi=300,xticks=1:c,yticks=(1:r,rack_codes(r)),xmirror=true,tickdirection=:none,tickfontsize=tickfontsize)
end 


function plot(pos::StackPosition;titlefontsize=18,tickfontsize=18)
    r,c = slots(pos) 
    plot(grid=false,size=(1200* r/4,800),yflip=false,legend=false,dpi=300,xticks=1:c,yticks=(1:r),xmirror=true,tickdirection=:none,tickfontsize=tickfontsize)
end 

function plot(pos::EmptyPosition;titlefontsize=18,tickfontsize=18)
    r,c=(1,1)
    plot(showaxis=false,grid=false,size=(1200,800),yflip=true,legend=false,dpi=300,ticks=:none)
end 

function plot(pos::UnconstrainedPosition;titlefontsize=18,tickfontsize=18)
    r,c= slots(pos)
    plot(grid=false,size=(1200* r/4,800*c/4),yflip=false,legend=false,dpi=300,xticks=1:c,yticks=(1:r),xmirror=true,tickdirection=:none,tickfontsize=tickfontsize)
end 
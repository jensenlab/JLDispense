
function preferred_ingredient_quantity(ing::JLIMS.Ingredient)
    pref_ing_quant=Dict(:solid=>u"g",:liquid=>u"µL",:organism=>u"OD*µL")
    return pref_ing_quant[ing.class]
end 

function preferred_stock_quantity(stock::JLIMS.Stock)
    pref_stock_quant=Dict(JLIMS.SolidStock=> u"g",JLIMS.LiquidStock => u"µL")
    return pref_stock_quant[typeof(stock)]
end 

function stock_concentration_array(stocks::Vector{T}) where T<: JLIMS.Stock
    ings=map(x->ingredients(x.composition),stocks)
    all_ings=unique(vcat(ings...))
    S=length(stocks)
    concentrations=DataFrame()
        for i in all_ings 
            concs=Any[]
            for s in 1:S
                if i in ings[s]
                    push!(concs,uconvert(preferred_ingredient_quantity(i)/preferred_stock_quantity(stocks[s]),stocks[s].composition.ingredients[i]))
                else 
                    push!(concs,0*preferred_ingredient_quantity(i)/preferred_stock_quantity(stocks[s]))
                end
            end 
            concentrations[:,Symbol(i.name)]=concs
        end 

    return concentrations
end 

            
function stock_quantity_array(stocks::Vector{T}) where T<:JLIMS.Stock
    ings=map(x->ingredients(x.composition),stocks)
    all_ings=unique(vcat(ings...))
    S=length(stocks)
    quantities=DataFrame()
    for i in all_ings 
        
        quants=Any[]
        for s in 1:S
            squant=stocks[s].quantity
            if i in ings[s]
                push!(quants,uconvert(preferred_ingredient_quantity(i),stocks[s].composition.ingredients[i]*squant))
            else 
                push!(quants,0*preferred_ingredient_quantity(i))
            end
        end 
        quantities[:,Symbol(i.name)]=quants
    end 
    return quantities 
end 


function stock_quantity_array(stock::JLIMS.Stock)
    return stock_quantity_array([stock])
end 

function stock_concentration_array(stock::JLIMS.Stock)
    return stock_concentration_array([stock])
end 





function mixer(sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Robot;directory=pwd())

    robot_type=typeof(robot) 
    tt=DataFrame() 
    protocol_name=""
    if robot_type == Human 
        tt,protocol_name=human(directory,sources,destinations,robot) 
    else
        error("$robot_type is unsupported by the mixer")
    end 
     
    return tt, protocol_name
end 


const string_unit_substitution= Dict(
    "%" => "percent",

    "" => "NoUnits"

)

const unit_string_substitution = Dict( 
    "percent" => "%",
    "NoUnits" => "" 
)


preferred_quantity_unit(ing::JLIMS.Solid) =u"mg"
preferred_quantity_unit(ing::JLIMS.Liquid) = u"µL"


preferred_quantity_unit(stock::JLIMS.Stock)= u"µL"
preferred_quantity_unit(stock::JLIMS.Mixture) =u"mg"




const chemical_access_dict=Dict(
    JLIMS.Solid => JLIMS.solids ,
    JLIMS.Liquid => JLIMS.liquids
)

PriorityDict = Dict{JLIMS.Chemical,UInt64}

const emptywell= Well{1}(0,"placeholder",nothing,Empty())

""" 
    concentration(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the concentration of an ingredient in a stock using the preferred units for that ingredient and stock 
    

"""
function concentration(stock::JLIMS.Stock,ingredient::JLIMS.Chemical) 
    if ingredient in stock 
        return uconvert(preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient]/JLIMS.quantity(stock))
    else
        return 0*preferred_quantity_unit(ingredient)/preferred_quantity_unit(stock)
    end 
end 


""" 
    quantity(stock::JLIMS.Stock,ingredient::JLIMS.Ingredient)

Return the quantity of an ingredient in a stock using the preferred units for that ingredient
    

"""
function quantity(stock::JLIMS.Stock,ingredient::JLIMS.Chemical)
    if ingredient in stock
        return uconvert(preferred_quantity_unit(ingredient),(chemical_access_dict[typeof(ingredient)])(stock)[ingredient])
    else
        return 0*preferred_quantity_unit(ingredient)
    end 
end 





function string_to_unit(str::String)

    if str in keys(string_unit_substitution) 
        return Unitful.uparse(string_unit_substitution[str])
    else
        return Unitful.uparse(str)
    end 
end 

function unit_to_string(unit::Unitful.Units) 
    ustr = string(unit) 
    if ustr in keys(unit_string_substitution)
        return unit_string_substitution[ustr]
    else 
        return ustr 
    end 
end 
    


function string_to_reagent(str::String,unit::Unitful.Units;chem_context::Vector{Module}=[JLIMS],kwargs...) 
    
    try 
        return chemparse(str;chem_context=chem_context)
    catch 
    end 

    try 
        return orgparse(str;chem_context=chem_context)
    catch 
    end 

    @warn("reagent $str not registered. parsing $str assuming it is a chemical. No chemical properties known.")
    if unit isa Unitful.DensityUnits || unit isa Unitful.MassUnits || unit isa Unitful.AmountUnits || unit isa Unitful.MolarityUnits 
        # a solid mass concentration in vc format or a mass in q format 
        return Solid(str,missing,missing,missing)
    elseif unit isa Unitful.DimensionlessUnits || unit isa Unitful.Volume # a %v/v concentration in vc format or a volume in q format 
        return Liquid(str,missing,missing,missing) 
    end 

end 


function reagent_to_string(chem::JLIMS.Chemical; chem_context::Vector{Module}=[JLIMS],kwargs...)
        all_chem_symbols = vcat(map(x->filter(y-> JLIMS.chemstr_check_bool(getfield(x,y)),names(x;all=true)),chem_context)...) # all chemical symbols 
        chem_strs = String.(all_chem_symbols)
        chem_dict = Dict(chemparse.(chem_strs;chem_context=chem_context)  .=> String.(chem_strs))
        if chem in keys(chem_dict) 
            return chem_dict[chem] # JLIMS registered chemicals appear in the chem dict. Their name is a display name and not necessarily the key used to find the chemical
        else
            return JLIMS.name(chem) # Chemicals that are not registered (generated on the fly) can be re-generated with just the name 
        end 
end 


function get_vc_reagents(vc::DataFrame,units::DataFrame;kwargs...)
    colnames =names(vc) 
    # in the vc format, each column other than the "volume" column is a reagent 
    reagent_names = setdiff(colnames,["volume"])
    # take units from first row to supply to the reagent parser 
    un = units[1,reagent_names]

    reagents = string_to_reagent.(reagent_names,string_to_unit.(collect(un));kwargs...) 

    return Dict(reagent_names .=> reagents)
end 

function get_quant_reagents(quant::DataFrame,units::DataFrame;kwargs...)
    reagent_names  =names(quant) 
    # take units from first row to supply to the reagent parser 
    un = units[1,reagent_names]

    reagents = string_to_reagent.(reagent_names,string_to_unit.(collect(un));kwargs...) 

    return Dict(reagent_names .=> reagents)
end 
    



function check_csv_inputs(df::DataFrame,units::DataFrame) 
    # check that units contains either 1 or N rows 
    N = nrow(df) 
    U = nrow(units) 
    if U != N && U != 1     # check that units contains either 1 or N rows 
        error("units file must either have a single row or a corresponding row for each stock")
    end 
    if names(df) != names(units)     # check that the column names of each table match 
        error("column names must be consistent for each table")
    end 
end 

function is_single_unit(units::DataFrame) 
    return nrow(units) == 1 
end 


function chemical_array(stocks::Vector{<:JLIMS.Stock},ingredients::Vector{<:JLIMS.Chemical};measure::Function=concentration,kwargs...) # can return concentration or quantity 
    out=DataFrame()
        for i in ingredients 
            vals=Any[]
            for s in stocks
                push!(vals,measure(s,i))
            end 
            out[:,reagent_to_string(i;kwargs...)]=vals
        end 

    return out
end 



## User formats 
# volconc : a volume x concentration table format for stocks 
# quant : a quantity based table format for defining stocks 
# stock: a JLIMS stock representation 

## Parsers for different dispensesolver inputs 

function volconc_to_stock(vc::DataFrame,units::DataFrame;kwargs...)
    check_csv_inputs(vc,units) 

    single_unit = is_single_unit(units) 

    reagent_dict = get_vc_reagents(vc,units;kwargs...)
    N = nrow(vc)
    stocks = Stock[]
    urow = 1 # assume single unit by default 
    for row in 1:N
        if !single_unit   # set the units file to follow the vc file, rather than use the single row 
            urow = row 
        end 

        st = Empty() 
        vol = vc[row,"volume"] * string_to_unit(units[urow,"volume"]) # parse the volume
        for r in keys(reagent_dict) 
            quant = ((vc[row,r] * string_to_unit(units[urow,r])) *vol)
            st +=  quant * reagent_dict[r] 
        end 
        push!(stocks,st)
    end 
    return stocks 
end 



function stock_to_volconc(stocks::Vector{<:Stock}; kwargs...)
    chems = get_all_chemicals(stocks) 
    concs = chemical_array(stocks,chems;measure=concentration,kwargs...)
    vols = JLIMS.quantity.(stocks) 

    vc = DataFrame()
    un = DataFrame() 
    vc.volume = ustrip.(vols) 
    un.volume = unit_to_string.(unit.(vols))

    vc= hcat(vc,ustrip.(concs))
    un = hcat(un,unit_to_string.(unit.(concs)))

    return vc, un 
end 


function quant_to_stock(quant::DataFrame,units::DataFrame;kwargs...)
    check_csv_inputs(quant,units)
    single_unit = is_single_unit(units) 

    reagent_dict = get_quant_reagents(quant,units;kwargs...)
    N = nrow(quant)
    stocks = Stock[]
    urow = 1 # assume single unit by default 
    for row in 1:N
        if !single_unit   # set the units file to follow the vc file, rather than use the single row 
            urow = row 
        end 

        st = Empty() 
        for r in keys(reagent_dict) 
            q = quant[row,r] * string_to_unit(units[urow,r])
            st +=  q * reagent_dict[r] 
        end 
        push!(stocks,st)
    end 
    return stocks 
end 

function stock_to_quant(stocks::Vector{<:Stock};kwargs...)
    chems = get_all_chemicals(stocks) 
    quants = chemical_array(stocks,chems;measure=quantity,kwargs...)

    q = ustrip.(quants) 
    un = unit_to_string.(unit.(quants) )
    return q, un 
end 


function volconc_to_quant(vc::DataFrame,units::DataFrame;kwargs...)
    s = volconc_to_stock(vc,units;kwargs...)
    return stock_to_quant(s;kwargs...)
end 

function quant_to_volconc(quant::DataFrame,units::DataFrame;kwargs...)
    s= quant_to_stock(quant,units;kwargs...)
    return stock_to_volconc(s;kwargs...) 
end 










        


function organism_array(stocks::Vector{<:JLIMS.Stock},orgs::Vector{JLIMS.Organism})
    pairs=Iterators.product(orgs,stocks) |> collect 
    out = Base.splat(in).(pairs)

    return DataFrame(out',JLIMS.name.(orgs)) # output expected bo be cultures by strains instead of strains by cultures
end  


function labware_id_and_type(lw::JLIMS.Labware)
    return (JLIMS.location_id(lw),typeof(lw))
end 

function get_all_labware(sources::Vector{<:JLIMS.Well},targets::Vector{<:JLIMS.Well})

    all(map(x-> !isnothing(x),JLIMS.parent.(sources))) || ArgumentError("all source wells must have a searchable parent")
    all(map(x-> !isnothing(x),JLIMS.parent.(targets))) || ArgumentError("all target wells must have a searchable parent")

    all_labware = unique(labware_id_and_type,vcat(JLIMS.parent.(sources),JLIMS.parent.(targets)))

    return all_labware 
end 


function get_all_chemicals(source::JLIMS.Stock)
         # Gather all ingredients contained in the sources, destinations, and priority list 
         solids = chemicals(JLIMS.solids(source))
         liqs = chemicals(JLIMS.liquids(source))
         return collect(union(solids,liqs))
end 

function get_all_chemicals(sources::Vector{<:JLIMS.Stock})
    return collect(union(get_all_chemicals.(sources)...))
end 

function get_all_chemicals(priority::PriorityDict)
    return collect((keys(priority)))
end 

function get_all_chemicals(sources::Vector{<:JLIMS.Stock},targets::Vector{<:JLIMS.Stock},priority::PriorityDict)
    return union(get_all_chemicals(sources),get_all_chemicals(targets),get_all_chemicals(priority))
end 


get_all_organisms(x::JLIMS.Stock) = collect(JLIMS.organisms(x))

function get_all_organisms(sources::Vector{<:JLIMS.Stock})
    return collect(union(get_all_organisms.(sources)...))
end 

function get_all_organisms(sources::Vector{<:JLIMS.Stock},targets::Vector{<:JLIMS.Stock})
    return collect(union(get_all_organisms(sources),get_all_organisms(targets)))
end 
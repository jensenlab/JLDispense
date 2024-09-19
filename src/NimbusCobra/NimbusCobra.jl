
function generate_temporary_labware(container::JLIMS.Container)

    n_wells=prod(container.shape)
    lw_id=string(UUIDs.uuid1())
    return JLIMS.Well.(1:n_wells,(lw_id,),1:n_wells,(container,))
end 



function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot_nimbus::Nimbus,robot_cobra::Cobra;lw_gen=generate_temporary_labware,kwargs...) where {T <: JLIMS.Stock,U <:JLIMS.Stock}

    compatible_containers=unique(vcat(map(x->x.compatible_containers,filter(x->x.is_source,robot_nimbus.properties.positions))...))

    srcs=filter(x->in(x.well.container,compatible_containers),sources)

    full_design=dispense_solver(srcs,destinations,human_default,minimize_overdrafts!,minimize_sources!,minimize_transfers!;pad=1.25,kwargs...)
    source_quantities=sum.(eachrow(full_design))

    pad_factor=1.3

    source_quantities=source_quantities*pad_factor
    deadvol=300u"µL"
    cap=dwp96_2ml.capacity-deadvol

    n_wells=prod(dwp96_2ml.shape)

    wells_needed=cld.(source_quantities,(cap,))

    plates_needed=cld(sum(wells_needed),n_wells)

    wells=vcat(lw_gen.([dwp96_2ml for _ in 1:plates_needed])...)
    wellidx=1
    intermediates=JLIMS.Stock[]
    for i in eachindex(source_quantities)
        filled=0u"µL"
        
        while filled < source_quantities[i]
            src=srcs[i]
            fillvol=min(wells[wellidx].container.capacity-deadvol,source_quantities[i])
            filled+= fillvol
            new_stock=JLIMS.Stock(src.composition,fillvol+deadvol,wells[wellidx])
            push!(intermediates,new_stock)
            wellidx+=1
        end 
    end 
    pn_nimbus,tt_nimbus=mixer(directory,srcs,intermediates,robot_nimbus;kwargs...)
    pn_cobra,tt_cobra=mixer(directory,intermediates,destinations,robot_cobra;kwargs...)

    return vcat(pn_nimbus,pn_cobra),vcat(tt_nimbus,tt_cobra)

end 
#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]
all_actives=JLD2.load("./src/Mixer/all_actives.jld2")["all_actives"]
intermediates=JLD2.load("./src/Mixer/intermediates.jdl2")["intermediates"]
intermediates=mixer("/Users/BDavid/Desktop/",all_actives,destinations,nimbus_default,cobra_default)

=#
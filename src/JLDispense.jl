module JLDispense

using 
    CSV,
    DataFrames

include("./container_names.jl")
include("./Cobra/Cobra.jl")
include("./Mantis/Mantis.jl")
include("./Tempest/Tempest.jl")
include("./Nimbus/Nimbus.jl")

export ContainerName,
containers

export CobraSettings,
SoftLinxSettings,
fill_protocol_template,
fill_softlinx_template,
snake_order,
CobraCSV,
fill_design,
design2protocols,
protocols2softlinx,
cobra


export nimbus,
convert_nimbus_design,
default_nimbus_config,
NimbusDispenseList,
NimbusRack

export mantis
export tempest, multi_tempest


end # module JLDispense
